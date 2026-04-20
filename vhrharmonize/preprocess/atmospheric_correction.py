"""Atmospheric correction step module."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import os
import re
import shutil
from typing import Any, Dict, List, Optional, Protocol, Tuple

import numpy as np
import rasterio
from tqdm import tqdm

from vhrharmonize.preprocess.fetch_external_data import fetch_power_atmosphere_for_bbox
from vhrharmonize.preprocess.helpers import log


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


def _py6s_atmosphere_profile(AtmosProfile, value: str):
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


def _py6s_aerosol_profile(AeroProfile, value: str):
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


def build_py6s_kwargs_from_standardized_metadata(
    metadata,
    *,
    ground_elevation_km: float,
    atmosphere_profile: str,
    aerosol_profile: str,
    aot550: float,
    visibility_km: Optional[float],
    water_vapor: float,
    ozone: float,
    sixs_executable: Optional[str],
    output_scale_factor: Optional[float],
    output_dtype: str,
    use_imd_radiance_calibration: bool,
    use_worldview_gain_offset_adjustment: bool,
) -> Dict[str, Any]:
    """Build shared Py6S kwargs from standardized metadata."""
    scene_dt = metadata.resolve_scene_datetime()
    py6s_kwargs = {
        "solar_zenith": metadata.sun_zenith,
        "solar_azimuth": metadata.sun_azimuth,
        "view_zenith": (
            metadata.view_zenith if metadata.view_zenith is not None else metadata.line_of_sight_zenith
        ),
        "view_azimuth": metadata.line_of_sight_azimuth,
        "day": scene_dt.day,
        "month": scene_dt.month,
        "ground_elevation_km": ground_elevation_km,
        "atmosphere_profile": atmosphere_profile,
        "aerosol_profile": aerosol_profile,
        "aot550": aot550,
        "visibility_km": visibility_km,
        "water_vapor": water_vapor,
        "ozone": ozone,
        "sixs_executable": sixs_executable,
        "output_scale_factor": output_scale_factor,
        "output_dtype": output_dtype,
    }
    dn_to_radiance_factors = list(getattr(metadata, "dn_to_radiance_factors", ()) or ())
    if use_imd_radiance_calibration and dn_to_radiance_factors:
        if use_worldview_gain_offset_adjustment:
            gains = list(getattr(metadata, "dn_to_radiance_gains", ()) or ())
            if gains and len(gains) >= len(dn_to_radiance_factors):
                dn_to_radiance_factors = [
                    float(dn_to_radiance_factors[idx]) * float(gains[idx])
                    for idx in range(len(dn_to_radiance_factors))
                ]
        py6s_kwargs["dn_to_radiance_factors"] = dn_to_radiance_factors
        dn_to_radiance_offsets = list(getattr(metadata, "dn_to_radiance_offsets", ()) or ())
        if use_worldview_gain_offset_adjustment and dn_to_radiance_offsets:
            py6s_kwargs["dn_to_radiance_offsets"] = dn_to_radiance_offsets
    return py6s_kwargs


def run_py6s(
    input_raster: str,
    output_raster: str,
    metadata,
    *,
    ground_elevation_km: float,
    atmosphere_profile: str,
    aerosol_profile: str,
    aot550: float,
    visibility_km: Optional[float],
    water_vapor: float,
    ozone: float,
    sixs_executable: Optional[str],
    output_scale_factor: Optional[float],
    output_dtype: str,
    use_imd_radiance_calibration: bool,
    use_worldview_gain_offset_adjustment: bool,
    auto_atmos_source: Optional[str] = None,
    bbox_wgs84: Optional[Tuple[float, float, float, float]] = None,
    auto_atmos_grid_size: int = 3,
    auto_atmos_search_days: int = 1,
    auto_atmos_timeout_s: float = 30.0,
    auto_atmos_power_endpoint: str = "https://power.larc.nasa.gov/api/temporal/daily/point",
    log_to_console: bool = False,
) -> Py6SRunResult:
    """Run Py6S correction using provider-neutral standardized metadata."""
    effective_aot550 = aot550
    effective_water_vapor = water_vapor
    effective_ozone = ozone
    auto_estimate = None

    if auto_atmos_source == "nasa_power" and atmosphere_profile.strip().lower() == "user":
        if bbox_wgs84 is None:
            raise ValueError("bbox_wgs84 is required for Py6S auto atmosphere fetch.")

        scene_dt = metadata.resolve_scene_datetime()
        estimate = fetch_power_atmosphere_for_bbox(
            day_utc=scene_dt.date(),
            min_lon=bbox_wgs84[0],
            min_lat=bbox_wgs84[1],
            max_lon=bbox_wgs84[2],
            max_lat=bbox_wgs84[3],
            grid_size=auto_atmos_grid_size,
            search_days=auto_atmos_search_days,
            timeout_s=auto_atmos_timeout_s,
            endpoint=auto_atmos_power_endpoint,
            log_to_console=log_to_console,
        )
        if estimate.aot550 is not None:
            effective_aot550 = float(estimate.aot550)
        if estimate.water_vapor is not None:
            effective_water_vapor = float(estimate.water_vapor)
        if estimate.ozone_cm_atm is not None:
            effective_ozone = float(estimate.ozone_cm_atm)
        auto_estimate = {
            "source": estimate.source,
            "date_used": estimate.date_used,
            "sample_count": estimate.sample_count,
            "aot550": estimate.aot550,
            "water_vapor": estimate.water_vapor,
            "ozone_cm_atm": estimate.ozone_cm_atm,
        }

    py6s_kwargs = build_py6s_kwargs_from_standardized_metadata(
        metadata,
        ground_elevation_km=ground_elevation_km,
        atmosphere_profile=atmosphere_profile,
        aerosol_profile=aerosol_profile,
        aot550=effective_aot550,
        visibility_km=visibility_km,
        water_vapor=effective_water_vapor,
        ozone=effective_ozone,
        sixs_executable=sixs_executable,
        output_scale_factor=output_scale_factor,
        output_dtype=output_dtype,
        use_imd_radiance_calibration=use_imd_radiance_calibration,
        use_worldview_gain_offset_adjustment=use_worldview_gain_offset_adjustment,
    )
    log("Running Py6S", enabled=log_to_console, step="py6s")
    Py6SCorrector().run(input_raster=input_raster, output_raster=output_raster, **py6s_kwargs)
    log("Wrote output", enabled=log_to_console, step="py6s")
    return Py6SRunResult(
        output_raster=output_raster,
        effective_params={
            "atmosphere_profile": atmosphere_profile,
            "aerosol_profile": aerosol_profile,
            "aot550": effective_aot550,
            "visibility_km": visibility_km,
            "water_vapor": effective_water_vapor,
            "ozone": effective_ozone,
            "output_scale_factor": output_scale_factor,
            "output_dtype": output_dtype,
            "use_imd_radiance_calibration": use_imd_radiance_calibration,
            "use_worldview_gain_offset_adjustment": use_worldview_gain_offset_adjustment,
            "auto_atmos_source": auto_atmos_source,
        },
        auto_atmos_estimate=auto_estimate,
    )


class Py6SCorrector:
    """Atmospheric correction adapter using Py6S."""

    def run(self, input_raster: str, output_raster: str, **kwargs: Any) -> str:
        """Run block-wise Py6S atmospheric correction for multispectral rasters."""
        try:
            from Py6S import AeroProfile, AtmosCorr, AtmosProfile, Geometry, SixS, Wavelength
        except ImportError as exc:
            raise RuntimeError(
                "Py6S is not installed. Install extras with `pip install -e \".[py6s]\"`."
            ) from exc

        required = ("solar_zenith", "solar_azimuth", "view_zenith", "view_azimuth", "day", "month")
        missing = [k for k in required if kwargs.get(k) is None]
        if missing:
            raise ValueError(f"Missing required Py6S kwargs: {', '.join(missing)}")

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

        band_wavelengths_um: List[float] = list(
            kwargs.get(
                "band_wavelengths_um",
                [0.4273, 0.4779, 0.5462, 0.6086, 0.6588, 0.7237, 0.8313, 0.9080],
            )
        )
        input_scale_factor = float(kwargs.get("input_scale_factor", 1.0))
        dn_to_radiance_factors = kwargs.get("dn_to_radiance_factors")
        dn_to_radiance_offsets = kwargs.get("dn_to_radiance_offsets")
        output_scale_factor = kwargs.get("output_scale_factor")
        output_dtype = str(kwargs.get("output_dtype", "float32"))
        clip_reflectance = bool(kwargs.get("clip_reflectance", True))
        nodata_override = kwargs.get("nodata_value")
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
                "6S executable not found. Install 6S and expose it in PATH, or pass "
                "`sixs_executable` (CLI: --py6s-executable)."
            )

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
            for band_idx in range(src.count):
                s.wavelength = Wavelength(float(band_wavelengths_um[band_idx]))
                s.run()
                coeffs.append((float(s.outputs.coef_xa), float(s.outputs.coef_xb), float(s.outputs.coef_xc)))

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
                    for _, window in src.block_windows(band_idx):
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


def init_envi_engine(envi_engine_path: str):
    """Initialize the ENVI task engine used by FLAASH."""
    import envipyengine.config
    from envipyengine import Engine

    envipyengine.config.set("engine", envi_engine_path)
    envi_engine = Engine("ENVI")
    envi_engine.tasks()
    return envi_engine


def validate_flaash_params(flaash_params: Dict[str, Any]) -> Dict[str, Any]:
    """Validate FLAASH params against the supported parameter list."""
    unknown_params = sorted(set(flaash_params) - FLAASH_ALLOWED_PARAMS)
    if unknown_params:
        raise ValueError(
            "Unsupported FLAASH parameter(s): "
            + ", ".join(unknown_params)
        )
    return flaash_params


def build_flaash_kwargs_from_standardized_metadata(
    input_raster: str,
    dem_file_path: str,
    footprint_vector_path: str,
    metadata,
    output_raster: str,
    *,
    dem_ground_percentile: float,
    modtran_atm: str,
    modtran_aer: str,
    use_aerosol: str,
    default_visibility: Optional[float],
    custom_params: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Build validated shared FLAASH params from standardized metadata."""
    from vhrharmonize.io.geospatial import get_image_percentile_value

    ground_elevation_m = get_image_percentile_value(
        dem_file_path,
        percentile=dem_ground_percentile,
        mask=footprint_vector_path,
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
        "SOLAR_AZIMUTH": metadata.sun_azimuth,
        "SOLAR_ZENITH": metadata.sun_zenith,
        "LOS_AZIMUTH": metadata.line_of_sight_azimuth,
        "LOS_ZENITH": metadata.line_of_sight_zenith,
        "OUTPUT_RASTER_URI": output_raster,
    }
    if custom_params:
        flaash_params.update(custom_params)
    flaash_params = {key: value for key, value in flaash_params.items() if value is not None}
    return validate_flaash_params(flaash_params)


def wsl_path_to_windows_for_envi(path: str) -> str:
    """Convert `/mnt/<drive>/...` WSL paths into Windows drive paths for ENVI."""
    match = re.match(r"^/mnt/([a-zA-Z])/(.*)$", path)
    if not match:
        raise ValueError(f"Path is not Windows-drive-backed via /mnt/<drive>/: {path}")
    drive = match.group(1).upper()
    rest = match.group(2).replace("/", "\\")
    return f"{drive}:\\{rest}"


def convert_flaash_params_paths_for_windows(flaash_params: Dict[str, Any]) -> Dict[str, Any]:
    """Convert FLAASH path params to Windows form for ENVI-on-Windows execution."""
    converted = dict(flaash_params)
    input_raster = converted.get("INPUT_RASTER")
    if isinstance(input_raster, dict) and "url" in input_raster:
        updated_input = dict(input_raster)
        updated_input["url"] = wsl_path_to_windows_for_envi(updated_input["url"])
        converted["INPUT_RASTER"] = updated_input

    for key in ("OUTPUT_RASTER_URI", "CLOUD_RASTER_URI", "WATER_RASTER_URI"):
        if converted.get(key):
            converted[key] = wsl_path_to_windows_for_envi(converted[key])
    return converted


def _execute_flaash_task(
    flaash_params,
    output_params_path,
    envi_engine,
    output_image_path_to_delete=None,
    *,
    log_to_console: bool = False,
):
    """Execute ENVI FLAASH with the provided parameter dictionary."""
    log("Running FLAASH", enabled=log_to_console, step="flaash")
    if output_image_path_to_delete:
        if os.path.exists(output_image_path_to_delete):
            os.remove(output_image_path_to_delete)
    task = envi_engine.task("FLAASH")
    task.execute(flaash_params)
    with open(output_params_path, "w", encoding="utf-8") as file:
        file.write(str(flaash_params))
    log("Wrote output", enabled=log_to_console, step="flaash")


def run_flaash_wrapper(args):
    """Executor wrapper for FLAASH grid runs."""
    test_params, test_output_params_path, envi_engine = args
    _execute_flaash_task(test_params, test_output_params_path, envi_engine)
    return test_params["OUTPUT_RASTER_URI"]


def parallel_flaash(test_flaash_params_array, envi_engine, max_workers=4):
    """Execute `run_flaash` in parallel on parameterized FLAASH runs."""
    all_output_paths = []
    tasks = [
        (test_params, test_output_params_path, envi_engine)
        for test_params, test_output_params_path in test_flaash_params_array
    ]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_flaash_wrapper, task) for task in tasks]
        for future in tqdm(as_completed(futures), total=len(futures), desc="FLAASH"):
            output_uri = future.result()
            all_output_paths.append(output_uri)

    return all_output_paths


class FLAASHCorrector:
    """Adapter around existing ENVI FLAASH execution helper."""

    def __init__(self, envi_engine) -> None:
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


def run_flaash(
    input_raster: str,
    output_raster: str,
    *,
    metadata=None,
    dem_file_path: Optional[str] = None,
    footprint_vector_path: Optional[str] = None,
    envi_engine_path: Optional[str] = None,
    envi_engine=None,
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
    """Run FLAASH using provider-neutral standardized metadata or explicit params."""
    if params is None:
        if metadata is None or dem_file_path is None or footprint_vector_path is None:
            raise ValueError(
                "metadata, dem_file_path, and footprint_vector_path are required when params is not provided."
            )
        params = build_flaash_kwargs_from_standardized_metadata(
            input_raster=input_raster,
            dem_file_path=dem_file_path,
            footprint_vector_path=footprint_vector_path,
            metadata=metadata,
            output_raster=output_raster,
            dem_ground_percentile=dem_ground_percentile,
            modtran_atm=modtran_atm,
            modtran_aer=modtran_aer,
            use_aerosol=use_aerosol,
            default_visibility=default_visibility,
            custom_params=custom_params,
        )
    else:
        params = validate_flaash_params(dict(params))
        params.setdefault("INPUT_RASTER", {"url": input_raster, "factory": "URLRaster"})
        params.setdefault("OUTPUT_RASTER_URI", output_raster)

    params_output_path = output_params_path or f"{output_raster}.flaash_params.txt"
    params_to_run = convert_flaash_params_paths_for_windows(params) if convert_paths_for_windows else params

    resolved_engine = envi_engine
    if resolved_engine is None:
        if not envi_engine_path:
            raise ValueError("envi_engine or envi_engine_path is required for FLAASH.")
        resolved_engine = init_envi_engine(envi_engine_path)

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


def atmospheric_correction(input_raster: str, output_raster: str, method: str = "flaash", **kwargs: Any) -> str:
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
    "Py6SRunResult",
    "FLAASHRunResult",
    "AtmosphericCorrector",
    "Py6SCorrector",
    "run_py6s",
    "build_py6s_kwargs_from_standardized_metadata",
    "init_envi_engine",
    "FLAASH_ALLOWED_PARAMS",
    "validate_flaash_params",
    "build_flaash_kwargs_from_standardized_metadata",
    "wsl_path_to_windows_for_envi",
    "convert_flaash_params_paths_for_windows",
    "run_flaash",
    "run_flaash_wrapper",
    "parallel_flaash",
    "FLAASHCorrector",
    "atmospheric_correction",
]
