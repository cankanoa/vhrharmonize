"""Atmospheric correction step module."""

from __future__ import annotations

import os
import re
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Protocol, Tuple

import numpy as np
import rasterio
from shapely.geometry.base import BaseGeometry
from vhrharmonize.io.progress import current_callback, progress, raster_windows, reports_progress
from vhrharmonize.io.progress_transport import local_worker_progress

from vhrharmonize.io.geospatial import get_image_percentile_value
from vhrharmonize.io.workflow_utils import remove_output_files
from .base import FunctionPlugin
from vhrharmonize.io.logging import _log, _logged_operation


@dataclass(frozen=True)
class Py6SRunResult:
    """Summary of a shared Py6S run."""

    output_raster: str
    effective_params: Dict[str, Any]
    auto_atmos_estimate: Optional[Dict[str, Any]]


@dataclass(frozen=True)
class FLAASHRunResult:
    """Summary of a shared FLAASH run."""

    output_raster: str
    params: Dict[str, Any]
    params_output_path: str


# ---------------------------------------------------------------------------
# Py6S
# ---------------------------------------------------------------------------


class AtmosphericCorrector(Protocol):
    """Atmospheric correction adapter interface."""

    def run(self, input_raster: str, output_raster: str, **kwargs: Any) -> str:
        """Run correction and return output path."""
        ...


def _py6s_atmosphere_profile(AtmosProfile: Any, value: str) -> Any:
    """Resolve a Py6S atmosphere profile.
    Args:
        AtmosProfile: Py6S AtmosProfile class.
        value: Atmosphere profile name.
    Returns:
        Py6S predefined atmosphere profile.
    """
    mapping = {
        "tropical": AtmosProfile.Tropical,
        "midlatitude_summer": AtmosProfile.MidlatitudeSummer,
        "midlatitude_winter": AtmosProfile.MidlatitudeWinter,
        "subarctic_summer": AtmosProfile.SubarcticSummer,
        "subarctic_winter": AtmosProfile.SubarcticWinter,
        "us_standard_1962": AtmosProfile.USStandard1962,
    }
    key = value.strip().lower()
    if key not in mapping:
        raise ValueError(f"Unsupported Py6S atmosphere profile: {value}")
    return AtmosProfile.PredefinedType(mapping[key])


def _py6s_aerosol_profile(AeroProfile: Any, value: str) -> Any:
    """Resolve a Py6S aerosol profile.
    Args:
        AeroProfile: Py6S AeroProfile class.
        value: Aerosol profile name.
    Returns:
        Py6S predefined aerosol profile.
    """
    mapping = {
        "continental": AeroProfile.Continental,
        "maritime": AeroProfile.Maritime,
        "urban": AeroProfile.Urban,
        "desert": AeroProfile.Desert,
        "biomass_burning": AeroProfile.BiomassBurning,
        "stratospheric": AeroProfile.Stratospheric,
    }
    key = value.strip().lower()
    if key not in mapping:
        raise ValueError(f"Unsupported Py6S aerosol profile: {value}")
    return AeroProfile.PredefinedType(mapping[key])


@_logged_operation("atmospheric_correction", inputs=("input_raster",), outputs=("output_raster",), worker_progress=True)
def run_py6s(
    input_raster: str,
    output_raster: str,
    *,
    variables: dict | None = None,
    solar_zenith: float | None = None,
    solar_azimuth: float | None = None,
    view_zenith: float | None = None,
    view_azimuth: float | None = None,
    day: int | None = None,
    month: int | None = None,
    band_wavelengths_um: list[float] | None = None,
    dn_to_radiance_factors: list[float] | None = None,
    dn_to_radiance_offsets: list[float] | None = None,
    dem_file_path: str | None = None,
    dem_ground_percentile: float = 50.0,
    footprint_geometry: dict | BaseGeometry | None = None,
    ground_elevation_km: float | None = None,
    atmosphere_profile: str = "midlatitude_summer",
    aerosol_profile: str = "maritime",
    aot550: float = 0.2,
    water_vapor: float = 2.5,
    ozone: float = 0.3,
    visibility_km: float | None = None,
    input_scale_factor: float = 1.0,
    output_scale_factor: float | None = None,
    output_dtype: str = "float32",
    clip_reflectance: bool = True,
    custom_nodata_value: float | None = None,
    sixs_executable: str | None = None,
    log_to_console: bool = False,
    scene_basename: str | None = None,
) -> Py6SRunResult:
    """Correct a raster using mapped variables and explicit spectral calibration.

    Args:
        input_raster: Input raster path.
        output_raster: Output corrected raster path.
        variables: Mapped JSON object; supplies geometry and omitted sensor fields.
        solar_zenith: Solar zenith in degrees; otherwise variables.solar_zenith.
        solar_azimuth: Solar azimuth in degrees; otherwise variables.solar_azimuth.
        view_zenith: View zenith in degrees; otherwise variables.view_zenith.
        view_azimuth: View azimuth in degrees; otherwise variables.view_azimuth.
        day: UTC acquisition day of month; otherwise variables.day.
        month: UTC acquisition month; otherwise variables.month.
        band_wavelengths_um: Band centers in micrometers, in raster band order.
        dn_to_radiance_factors: Per-band multiplicative calibration; otherwise mapped variables.
        dn_to_radiance_offsets: Per-band additive calibration; otherwise mapped variables.
        dem_file_path: Optional DEM in meters for sampling ground elevation.
        dem_ground_percentile: DEM percentile in [0, 100].
        footprint_geometry: GeoJSON or Shapely footprint in WGS84; otherwise variables.geometry.
        ground_elevation_km: Explicit elevation; None samples the DEM, or uses 0 without a DEM.
        atmosphere_profile: user, tropical, midlatitude_summer, midlatitude_winter, subarctic_summer, subarctic_winter.
        aerosol_profile: maritime, continental, urban, desert, biomass_burning, stratospheric.
        aot550: Aerosol optical thickness at 550 nm; used when visibility_km is None.
        water_vapor: Water vapor in g/cm^2 for the user atmosphere profile.
        ozone: Ozone in cm-atm for the user atmosphere profile.
        visibility_km: Visibility override in kilometers.
        input_scale_factor: Multiplicative input radiance scale.
        output_scale_factor: Reflectance scale before writing; None leaves reflectance unscaled.
        output_dtype: Raster output dtype.
        clip_reflectance: Clip reflectance to [0, 1] before scaling.
        custom_nodata_value: Output nodata override; None inherits the source.
        sixs_executable: 6S executable; None checks SIXS_EXECUTABLE and PATH.
        log_to_console: Emit processing logs.
        scene_basename: Optional scene name in logs.
    """
    if not output_raster:
        raise ValueError("output_raster must be an explicit output path")
    output_raster = str(output_raster)
    if not 0 <= dem_ground_percentile <= 100:
        raise ValueError("dem_ground_percentile must be in [0, 100].")
    effective = dict(locals())
    for key in (
        "input_raster",
        "output_raster",
        "variables",
        "dem_file_path",
        "dem_ground_percentile",
        "footprint_geometry",
        "log_to_console",
        "scene_basename",
    ):
        effective.pop(key)
    for key in (
        "solar_zenith",
        "solar_azimuth",
        "view_zenith",
        "view_azimuth",
        "day",
        "month",
        "band_wavelengths_um",
        "dn_to_radiance_factors",
        "dn_to_radiance_offsets",
    ):
        if effective[key] is None:
            effective[key] = (variables or {}).get(key)
    if ground_elevation_km is None:
        effective["ground_elevation_km"] = 0.0
        if dem_file_path:
            from vhrharmonize.io.metadata import materialize_geometry

            geometry = (
                footprint_geometry
                if footprint_geometry is not None
                else (variables or {}).get("geometry")
            )
            if isinstance(geometry, dict):
                geometry = materialize_geometry(geometry)
            effective["ground_elevation_km"] = (
                get_image_percentile_value(
                    dem_file_path,
                    percentile=dem_ground_percentile,
                    mask=geometry,
                )
                / 1000.0
            )
    Py6SCorrector().run(input_raster=input_raster, output_raster=output_raster, **effective)
    return Py6SRunResult(
        output_raster=output_raster, effective_params=effective, auto_atmos_estimate=None
    )


class Py6SCorrector:
    """Atmospheric correction adapter using Py6S."""

    def run(self, input_raster: str, output_raster: str, **kwargs: Any) -> str:
        """Run block-wise Py6S atmospheric correction for multispectral rasters."""
        required = (
            "solar_zenith",
            "solar_azimuth",
            "view_zenith",
            "view_azimuth",
            "day",
            "month",
            "band_wavelengths_um",
        )
        missing = [k for k in required if kwargs.get(k) is None]
        if missing:
            raise ValueError(f"Missing required Py6S kwargs: {', '.join(missing)}")
        if os.path.exists(output_raster) and os.path.samefile(input_raster, output_raster):
            raise ValueError("Py6S input and output rasters must be different files.")

        solar_zenith = float(kwargs["solar_zenith"])
        solar_azimuth = float(kwargs["solar_azimuth"])
        view_zenith = float(kwargs["view_zenith"])
        view_azimuth = float(kwargs["view_azimuth"])
        day = int(kwargs["day"])
        month = int(kwargs["month"])

        ground_elevation_km = float(kwargs.get("ground_elevation_km", 0.0))
        atmosphere_profile = str(kwargs.get("atmosphere_profile", "midlatitude_summer"))
        aerosol_profile = str(kwargs.get("aerosol_profile", "maritime"))
        aot550 = float(kwargs.get("aot550", 0.2))
        water_vapor = float(kwargs.get("water_vapor", 2.5))
        ozone = float(kwargs.get("ozone", 0.3))

        band_wavelengths_um = list(kwargs["band_wavelengths_um"])
        input_scale_factor = float(kwargs.get("input_scale_factor", 1.0))
        dn_to_radiance_factors = kwargs.get("dn_to_radiance_factors")
        dn_to_radiance_offsets = kwargs.get("dn_to_radiance_offsets")
        output_scale_factor = kwargs.get("output_scale_factor")
        output_dtype = str(kwargs.get("output_dtype", "float32"))
        clip_reflectance = bool(kwargs.get("clip_reflectance", True))
        nodata_override = kwargs.get("custom_nodata_value")
        visibility_km = kwargs.get("visibility_km")
        sixs_executable = kwargs.get("sixs_executable")

        if sixs_executable is None:
            sixs_executable = (
                os.environ.get("SIXS_EXECUTABLE")
                or shutil.which("sixs")
                or shutil.which("sixsV1.1")
            )
        if sixs_executable is None:
            raise RuntimeError(
                "6S executable not found. Install it with `conda install conda-forge::sixs`, "
                "expose it in PATH, or set sixs_executable in the step parameters."
            )

        from Py6S import AeroProfile, AtmosCorr, AtmosProfile, Geometry, SixS, Wavelength

        s = SixS(path=sixs_executable)
        s.geometry = Geometry.User()
        s.geometry.solar_z = solar_zenith
        s.geometry.solar_a = solar_azimuth
        s.geometry.view_z = view_zenith
        s.geometry.view_a = view_azimuth
        s.geometry.day = day
        s.geometry.month = month
        s.altitudes.set_sensor_satellite_level()
        s.altitudes.set_target_custom_altitude(ground_elevation_km)
        s.atmos_profile = AtmosProfile.UserWaterAndOzone(water_vapor, ozone)
        if atmosphere_profile.strip().lower() != "user":
            s.atmos_profile = _py6s_atmosphere_profile(AtmosProfile, atmosphere_profile)
        s.aero_profile = _py6s_aerosol_profile(AeroProfile, aerosol_profile)
        if visibility_km is not None:
            s.aot550 = None
            s.visibility = float(visibility_km)
        else:
            s.visibility = None
            s.aot550 = aot550
        s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(0.2)

        with rasterio.open(input_raster) as src:
            if src.count > len(band_wavelengths_um):
                raise ValueError(
                    f"Not enough band wavelengths for source bands. bands={src.count}, "
                    f"wavelengths={len(band_wavelengths_um)}"
                )
            if dn_to_radiance_factors is not None and src.count > len(dn_to_radiance_factors):
                raise ValueError(
                    "Not enough dn_to_radiance_factors for source bands. "
                    f"bands={src.count}, factors={len(dn_to_radiance_factors)}"
                )
            if dn_to_radiance_offsets is not None and src.count > len(dn_to_radiance_offsets):
                raise ValueError(
                    "Not enough dn_to_radiance_offsets for source bands. "
                    f"bands={src.count}, offsets={len(dn_to_radiance_offsets)}"
                )

            nodata = src.nodata if nodata_override is None else float(nodata_override)
            profile = src.profile.copy()
            profile.update(
                dtype=output_dtype,
                nodata=nodata,
                compress="lzw",
                BIGTIFF="IF_SAFER",
            )
            input_tags = src.tags()
            input_rpc_tags = src.tags(ns="RPC")
            input_gcps, input_gcps_crs = src.gcps

            coeffs = []
            for band_idx in progress(range(src.count), desc="Calculating Py6S coefficients", unit="bands"):
                s.wavelength = Wavelength(float(band_wavelengths_um[band_idx]))
                s.run()
                coeffs.append(
                    (float(s.outputs.coef_xa), float(s.outputs.coef_xb), float(s.outputs.coef_xc))
                )

            os.makedirs(os.path.dirname(output_raster) or ".", exist_ok=True)
            with rasterio.open(output_raster, "w", **profile) as dst:
                if input_tags:
                    dst.update_tags(**input_tags)
                if input_rpc_tags:
                    dst.update_tags(ns="RPC", **input_rpc_tags)
                if input_gcps:
                    dst.gcps = (input_gcps, input_gcps_crs)

                for band_idx in range(1, src.count + 1):
                    xa, xb, xc = coeffs[band_idx - 1]
                    if dn_to_radiance_factors is not None:
                        band_scale = float(dn_to_radiance_factors[band_idx - 1])
                    else:
                        band_scale = input_scale_factor
                    band_offset = (
                        float(dn_to_radiance_offsets[band_idx - 1])
                        if dn_to_radiance_offsets is not None
                        else 0.0
                    )
                    for _, window in raster_windows(src, band_idx, desc=f"Correcting band {band_idx}"):
                        in_block = src.read(band_idx, window=window).astype(np.float32)
                        mask = src.read_masks(band_idx, window=window) > 0
                        if nodata is not None:
                            mask &= in_block != nodata

                        radiance = in_block * band_scale + band_offset
                        y = xa * radiance - xb
                        corrected = y / (1.0 + xc * y)

                        if clip_reflectance:
                            corrected = np.clip(corrected, 0.0, 1.0)

                        if output_scale_factor is not None:
                            corrected = corrected * float(output_scale_factor)

                        if nodata is not None:
                            corrected = np.where(mask, corrected, nodata)

                        dst.write(corrected.astype(output_dtype), band_idx, window=window)

        return output_raster


# ---------------------------------------------------------------------------
# FLAASH
# ---------------------------------------------------------------------------

FLAASH_ALLOWED_PARAMS = {
    "SENSOR_TYPE",
    "INPUT_SCALE",
    "OUTPUT_SCALE",
    "CALIBRATION_FILE",
    "CALIBRATION_FORMAT",
    "CALIBRATION_UNITS",
    "LAT_LONG",
    "SENSOR_ALTITUDE",
    "DATE_TIME",
    "USE_ADJACENCY",
    "DEFAULT_VISIBILITY",
    "USE_POLISHING",
    "POLISHING_RESOLUTION",
    "SENSOR_AUTOCALIBRATION",
    "SENSOR_CAL_PRECISION",
    "SENSOR_CAL_FEATURE_LIST",
    "GROUND_ELEVATION",
    "SOLAR_AZIMUTH",
    "SOLAR_ZENITH",
    "LOS_AZIMUTH",
    "LOS_ZENITH",
    "IFOV",
    "MODTRAN_ATM",
    "MODTRAN_AER",
    "MODTRAN_RES",
    "MODTRAN_MSCAT",
    "CO2_MIXING",
    "WATER_ABS_CHOICE",
    "WATER_MULT",
    "WATER_VAPOR_PRESET",
    "USE_AEROSOL",
    "AEROSOL_SCALE_HT",
    "AER_BAND_RATIO",
    "AER_BAND_WAVL",
    "AER_REFERENCE_VALUE",
    "AER_REFERENCE_PIXEL",
    "AER_BANDLOW_WAVL",
    "AER_BANDLOW_MAXREFL",
    "AER_BANDHIGH_WAVL",
    "AER_BANDHIGH_MAXREFL",
    "INPUT_RASTER",
    "OUTPUT_RASTER_URI",
    "CLOUD_RASTER_URI",
    "WATER_RASTER_URI",
}


def _init_envi_engine(envi_engine_path: str) -> Any:
    """Initialize the ENVI task engine used by FLAASH."""

    import envipyengine.config
    from envipyengine import Engine

    envipyengine.config.set("engine", envi_engine_path)
    envi_engine = Engine("ENVI")
    envi_engine.tasks()
    return envi_engine


def _validate_flaash_params(flaash_params: Dict[str, Any]) -> Dict[str, Any]:
    """Validate FLAASH params against the supported parameter list."""
    unknown_params = sorted(set(flaash_params) - FLAASH_ALLOWED_PARAMS)
    if unknown_params:
        raise ValueError("Unsupported FLAASH parameter(s): " + ", ".join(unknown_params))
    return flaash_params


def _build_flaash_kwargs_from_variables(
    input_raster: str,
    dem_file_path: str,
    footprint_geometry: BaseGeometry,
    variables: Any,
    output_raster: str,
    *,
    dem_ground_percentile: float,
    modtran_atm: str,
    modtran_aer: str,
    use_aerosol: str,
    default_visibility: Optional[float],
    custom_params: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Build validated shared FLAASH parameters.
    Args:
        input_raster: Input raster path.
        dem_file_path: DEM raster path.
        footprint_geometry: Shapely footprint in EPSG:4326 for DEM sampling.
        variables: Mapped variables JSON object.
        output_raster: Output raster path.
        dem_ground_percentile: DEM percentile used for ground elevation estimation.
        modtran_atm: MODTRAN atmosphere profile name.
        modtran_aer: MODTRAN aerosol profile name.
        use_aerosol: FLAASH aerosol handling mode.
        default_visibility: Optional default visibility override.
        custom_params: Optional custom FLAASH parameter overrides.
    Returns:
        Validated FLAASH parameter dictionary.
    """

    ground_elevation_m = get_image_percentile_value(
        dem_file_path,
        percentile=dem_ground_percentile,
        mask=footprint_geometry,
    )
    flaash_params = {
        "INPUT_RASTER": {"url": input_raster, "factory": "URLRaster"},
        "MODTRAN_ATM": modtran_atm,
        "MODTRAN_AER": modtran_aer,
        "MODTRAN_RES": 5.0,
        "MODTRAN_MSCAT": "DISORT",
        "USE_AEROSOL": use_aerosol,
        "DEFAULT_VISIBILITY": default_visibility,
        "AER_BAND_RATIO": 0.5,
        "AER_BANDLOW_WAVL": 425,
        "AER_BANDHIGH_WAVL": 660,
        "AER_BANDHIGH_MAXREFL": 0.2,
        "GROUND_ELEVATION": ground_elevation_m / 1000.0,
        "SOLAR_AZIMUTH": variables["solar_azimuth"],
        "SOLAR_ZENITH": variables["solar_zenith"],
        "LOS_AZIMUTH": variables["line_of_sight_azimuth"],
        "LOS_ZENITH": variables["line_of_sight_zenith"],
        "OUTPUT_RASTER_URI": output_raster,
    }
    if custom_params:
        flaash_params.update(custom_params)
    flaash_params = {key: value for key, value in flaash_params.items() if value is not None}
    return _validate_flaash_params(flaash_params)


def _wsl_path_to_windows_for_envi(path: str) -> str:
    """Convert `/mnt/<drive>/...` WSL paths into Windows drive paths for ENVI."""
    match = re.match(r"^/mnt/([a-zA-Z])/(.*)$", path)
    if not match:
        raise ValueError(f"Path is not Windows-drive-backed via /mnt/<drive>/: {path}")
    drive = match.group(1).upper()
    rest = match.group(2).replace("/", "\\")
    return f"{drive}:\\{rest}"


def _convert_flaash_params_paths_for_windows(flaash_params: Dict[str, Any]) -> Dict[str, Any]:
    """Convert FLAASH path params to Windows form for ENVI-on-Windows execution."""
    converted = dict(flaash_params)
    input_raster = converted.get("INPUT_RASTER")
    if isinstance(input_raster, dict) and "url" in input_raster:
        updated_input = dict(input_raster)
        updated_input["url"] = _wsl_path_to_windows_for_envi(updated_input["url"])
        converted["INPUT_RASTER"] = updated_input

    for key in ("OUTPUT_RASTER_URI", "CLOUD_RASTER_URI", "WATER_RASTER_URI"):
        if converted.get(key):
            converted[key] = _wsl_path_to_windows_for_envi(converted[key])
    return converted


def _execute_flaash_task(
    flaash_params: Dict[str, Any],
    output_params_path: str,
    envi_engine: Any,
    output_image_path_to_delete: Optional[str] = None,
    *,
    log_to_console: bool = False,
) -> None:
    """Execute ENVI FLAASH with the provided parameter dictionary."""
    _log("Processing", enabled=log_to_console, step="flaash")
    if output_image_path_to_delete:
        remove_output_files([output_image_path_to_delete])
    task = envi_engine.task("FLAASH")
    task.execute(flaash_params)
    with open(output_params_path, "w", encoding="utf-8") as file:
        file.write(str(flaash_params))
    _log("Wrote output", enabled=log_to_console, step="flaash")


@reports_progress
def _run_flaash_wrapper(args: tuple[Dict[str, Any], str, Any]) -> str:
    """Executor wrapper for FLAASH grid runs."""
    test_params, test_output_params_path, envi_engine = args
    _execute_flaash_task(test_params, test_output_params_path, envi_engine)
    return test_params["OUTPUT_RASTER_URI"]


@reports_progress(worker_progress=True)
def parallel_flaash(
    test_flaash_params_array: List[tuple[Dict[str, Any], str]],
    envi_engine: Any,
    max_workers: int = 4,
) -> List[str]:
    """Execute `run_flaash` in parallel on parameterized FLAASH runs."""
    all_output_paths = []
    tasks = [
        (test_params, test_output_params_path, envi_engine)
        for test_params, test_output_params_path in test_flaash_params_array
    ]

    with (
        local_worker_progress(current_callback(), processes=True) as reporter,
        ProcessPoolExecutor(max_workers=max_workers) as executor,
    ):
        callback_kwargs = {"progress_callback": reporter} if reporter is not None else {}
        futures = [executor.submit(_run_flaash_wrapper, task, **callback_kwargs) for task in tasks]
        for future in progress(as_completed(futures), total=len(futures), desc="FLAASH", disable=False):
            output_uri = future.result()
            all_output_paths.append(output_uri)

    return all_output_paths


class FLAASHCorrector:
    """Adapter around existing ENVI FLAASH execution helper."""

    def __init__(self, envi_engine: Any) -> None:
        """Initialize the FLAASH adapter.
        Args:
            envi_engine: Initialized ENVI engine instance.
        Returns:
            None.
        """
        self.envi_engine = envi_engine

    def run(self, input_raster: str, output_raster: str, **kwargs: Any) -> str:
        """Run FLAASH correction and return the output raster path."""
        params: Dict[str, Any] = dict(kwargs)
        params.setdefault("INPUT_RASTER", {"url": input_raster, "factory": "URLRaster"})
        params.setdefault("OUTPUT_RASTER_URI", output_raster)
        params_path = kwargs.get("params_path")
        if not params_path:
            raise ValueError("params_path is required for FLAASHCorrector.run")
        _execute_flaash_task(params, params_path, self.envi_engine, output_raster)
        return output_raster


@_logged_operation("atmospheric_correction", inputs=("input_raster",), outputs=("output_raster",))
def run_flaash(
    input_raster: str,
    output_raster: str,
    *,
    variables: Any = None,
    dem_file_path: Optional[str] = None,
    footprint_geometry: dict | BaseGeometry | None = None,
    envi_engine_path: Optional[str] = None,
    envi_engine: Any = None,
    output_params_path: Optional[str] = None,
    convert_paths_for_windows: bool = False,
    delete_output_before_run: Optional[str] = None,
    params: Optional[Dict[str, Any]] = None,
    dem_ground_percentile: float = 50.0,
    modtran_atm: str = "Mid-Latitude Summer",
    modtran_aer: str = "Rural",
    use_aerosol: str = "AUTO",
    default_visibility: Optional[float] = None,
    custom_params: Optional[Dict[str, Any]] = None,
    log_to_console: bool = False,
) -> FLAASHRunResult:
    """Run FLAASH using mapped variables JSON or explicit params.
    Args:
        input_raster: Input raster path.
        output_raster: Output raster path.
        variables: Optional mapped variables JSON object.
        dem_file_path: Optional DEM raster path.
        footprint_geometry: Optional Shapely footprint in EPSG:4326.
        envi_engine_path: Optional ENVI engine executable path.
        envi_engine: Optional initialized ENVI engine.
        output_params_path: Optional executed-params output path.
        convert_paths_for_windows: Whether to convert paths for Windows ENVI execution.
        delete_output_before_run: Optional output path to delete before execution.
        params: Optional explicit FLAASH parameter dictionary.
        dem_ground_percentile: DEM percentile used for ground elevation estimation.
        modtran_atm: MODTRAN atmosphere profile name.
        modtran_aer: MODTRAN aerosol profile name.
        use_aerosol: FLAASH aerosol handling mode.
        default_visibility: Optional default visibility override.
        custom_params: Optional custom FLAASH parameter overrides.
        log_to_console: Whether to emit console logs.
    Returns:
        FLAASH run summary.
    """
    if not output_raster:
        raise ValueError("output_raster must be an explicit output path")
    output_raster = str(output_raster)
    if footprint_geometry is None and variables is not None:
        footprint_geometry = variables.get("geometry")
    if isinstance(footprint_geometry, dict):
        from vhrharmonize.io.metadata import materialize_geometry

        footprint_geometry = materialize_geometry(footprint_geometry)
    if params is None:
        if variables is None or dem_file_path is None or footprint_geometry is None:
            raise ValueError(
                "variables, dem_file_path, and footprint_geometry are required when params is not provided."
            )
        params = _build_flaash_kwargs_from_variables(
            input_raster=input_raster,
            dem_file_path=dem_file_path,
            footprint_geometry=footprint_geometry,
            variables=variables,
            output_raster=output_raster,
            dem_ground_percentile=dem_ground_percentile,
            modtran_atm=modtran_atm,
            modtran_aer=modtran_aer,
            use_aerosol=use_aerosol,
            default_visibility=default_visibility,
            custom_params=custom_params,
        )
    else:
        params = _validate_flaash_params(dict(params))
        params.setdefault("INPUT_RASTER", {"url": input_raster, "factory": "URLRaster"})
        params.setdefault("OUTPUT_RASTER_URI", output_raster)

    params_output_path = output_params_path or f"{output_raster}.flaash_params.txt"
    params_to_run = (
        _convert_flaash_params_paths_for_windows(params) if convert_paths_for_windows else params
    )

    resolved_engine = envi_engine
    if resolved_engine is None:
        if not envi_engine_path:
            raise ValueError("envi_engine or envi_engine_path is required for FLAASH.")
        resolved_engine = _init_envi_engine(envi_engine_path)

    _execute_flaash_task(
        params_to_run,
        params_output_path,
        resolved_engine,
        delete_output_before_run or output_raster,
        log_to_console=log_to_console,
    )
    return FLAASHRunResult(
        output_raster=output_raster,
        params=params,
        params_output_path=params_output_path,
    )


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


@reports_progress(worker_progress=True)
def atmospheric_correction(
    input_raster: str, output_raster: str, method: str = "flaash", **kwargs: Any
) -> str:
    """Run atmospheric correction using the requested backend."""
    method_norm = method.strip().lower()
    if method_norm == "flaash":
        envi_engine = kwargs.pop("envi_engine", None)
        if envi_engine is None:
            raise ValueError("`envi_engine` is required when method='flaash'.")
        runner: AtmosphericCorrector = FLAASHCorrector(envi_engine=envi_engine)
        return runner.run(input_raster=input_raster, output_raster=output_raster, **kwargs)

    if method_norm == "py6s":
        runner = Py6SCorrector()
        return runner.run(input_raster=input_raster, output_raster=output_raster, **kwargs)

    raise ValueError(f"Unsupported atmospheric correction method: {method}")


__all__ = [
    "AtmosphericCorrection",
    "Py6SRunResult",
    "FLAASHRunResult",
    "AtmosphericCorrector",
    "Py6SCorrector",
    "run_py6s",
    "FLAASH_ALLOWED_PARAMS",
    "run_flaash",
    "parallel_flaash",
    "FLAASHCorrector",
    "atmospheric_correction",
]


class AtmosphericCorrection(FunctionPlugin):
    input_dependency_paths = frozenset({"dem_file_path", "input_raster"})
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset({"output_params_path", "output_raster"})
    output_dependency_paths = output_path_resolution_paths
    output_parent_creation_paths = output_path_resolution_paths
    output_collision_check_paths = output_path_resolution_paths
    output_reuse_paths = output_path_resolution_paths
    output_validation_paths = output_path_resolution_paths
    output_invalid_removal_paths = output_path_resolution_paths
    output_overview_calculation_paths = frozenset({"output_raster"})
    output_temporary_cleanup_paths = output_path_resolution_paths
    output_hpc_staging_paths = output_path_resolution_paths
    output_hpc_download_paths = output_path_resolution_paths

    options = {
        "solar_zenith",
        "solar_azimuth",
        "view_zenith",
        "view_azimuth",
        "day",
        "month",
        "ground_elevation_km",
        "atmosphere_profile",
        "aerosol_profile",
        "aot550",
        "water_vapor",
        "ozone",
        "band_wavelengths_um",
        "input_scale_factor",
        "dn_to_radiance_factors",
        "dn_to_radiance_offsets",
        "output_scale_factor",
        "output_dtype",
        "clip_reflectance",
        "custom_nodata_value",
        "visibility_km",
        "sixs_executable",
    }

    def run(self, *, params, shared):
        params = dict(params)
        method = params.pop("method", "py6s")
        if method not in {"py6s", "flaash"}:
            raise ValueError(f"Unknown atmospheric correction method: {method}")
        self.target = f"vhrharmonize.plugins.atmospheric_correction:run_{method}"
        if method == "flaash":
            self.options = set()
            self.aliases = {}
        return super().run(params=params, shared=shared)
