"""Py6S atmospheric correction with explicit spectral calibration."""

from __future__ import annotations

import csv
import os
import shutil
from dataclasses import dataclass
from typing import Any, Dict, Optional

import numpy as np
import rasterio
from shapely.geometry.base import BaseGeometry

from vhrharmonize.io.geospatial import get_image_percentile_value
from vhrharmonize.io.logging import _logged_operation
from vhrharmonize.io.progress import progress, raster_windows
from .base import FunctionPlugin


@dataclass(frozen=True)
class Py6SRunResult:
    """Summary of a shared Py6S run."""

    output_raster: str
    effective_params: Dict[str, Any]
    auto_atmos_estimate: Optional[Dict[str, Any]]


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


_PY6S_SPECTRAL_PARAMETERS = (
    "wavelengths_center_um", "wavelengths_bounds_um", "wavelengths_response_csv",
)


def _py6s_spectral_parameter(params):
    """Require one spectral representation, including fields supplied by variables."""
    selected = [name for name in _PY6S_SPECTRAL_PARAMETERS if params.get(name) is not None]
    if len(selected) != 1:
        raise ValueError("Specify exactly one of " + ", ".join(_PY6S_SPECTRAL_PARAMETERS) + ".")
    return selected[0]


def _py6s_wavelength_number(value):
    if isinstance(value, (bool, str)) or not np.isscalar(value):
        raise ValueError("Spectral wavelengths must be finite numbers in micrometers.")
    value = float(value)
    if not np.isfinite(value) or not 0.2 <= value <= 4.0:
        raise ValueError("Spectral wavelengths must be between 0.2 and 4.0 micrometers.")
    return value


def _py6s_response_csv(filename):
    """Read original wavelength samples and interpolate each response for 6S."""
    with open(filename, newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, strict=True)
        header = [name.strip() for name in next(reader, [])]
        if (len(header) < 2 or header[0] != "wavelength_nm"
                or any(not name for name in header) or len(set(header)) != len(header)):
            raise ValueError("Response CSV header must be wavelength_nm followed by unique band names in raster order.")
        rows = []
        for row in reader:
            if not row:
                continue
            try:
                if len(row) != len(header):
                    raise ValueError("wrong number of columns")
                rows.append([float(value) for value in row])
            except ValueError as error:
                raise ValueError(f"Invalid response CSV row {reader.line_num}: {error}") from error
    data = np.asarray(rows, dtype=float)
    if len(rows) < 2 or not np.isfinite(data).all():
        raise ValueError("Response CSV requires at least two rows of finite numeric values.")
    wavelengths = data[:, 0] / 1000.0
    response = data[:, 1:]
    if (wavelengths[0] < 0.2 or wavelengths[-1] > 4.0
            or not (np.diff(wavelengths) > 0).all()):
        raise ValueError("Response CSV wavelengths must be strictly increasing within 200–4000 nm.")
    if (response < 0).any() or not (response > 0).any(axis=0).all():
        raise ValueError("Response CSV bands must have nonnegative responses and at least one positive value each.")
    # Keep the CSV's original sampling. Only the runtime filter uses the 6S
    # 2.5 nm grid; interpolation uses the actual source wavelength coordinates.
    start = np.floor(wavelengths[0] / 0.0025 + 1e-10) * 0.0025
    stop = np.ceil(wavelengths[-1] / 0.0025 - 1e-10) * 0.0025
    grid = np.round(start + np.arange(int(round((stop - start) / 0.0025)) + 1) * 0.0025, 12)
    spectra = []
    for name, values in zip(header[1:], response.T):
        sampled = np.interp(grid, wavelengths, values, left=0.0, right=0.0)
        if not (sampled > 0).any():
            raise ValueError(f"Response CSV band {name} has no positive samples on the 6S 2.5 nm grid.")
        sampled /= sampled.max()
        spectra.append((float(start), float(stop), sampled.tolist()))
    return spectra


def _py6s_spectral_bands(params):
    """Return native Wavelength arguments in the input raster's band order."""
    name = _py6s_spectral_parameter(params)
    value = params[name]
    if name == "wavelengths_response_csv":
        return _py6s_response_csv(value)
    if name == "wavelengths_center_um":
        values = [value] if np.isscalar(value) else value
        if not isinstance(values, (list, tuple)) or not values:
            raise ValueError("wavelengths_center_um must contain one center per raster band.")
        return [(_py6s_wavelength_number(center),) for center in values]
    if not isinstance(value, (list, tuple)) or not value:
        raise ValueError("wavelengths_bounds_um must contain [lower, upper] for each raster band.")
    values = [value] if len(value) == 2 and all(np.isscalar(v) for v in value) else value
    spectra = []
    for bounds in values:
        if not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
            raise ValueError("Each wavelength interval must contain [lower, upper].")
        lower, upper = map(_py6s_wavelength_number, bounds)
        if lower >= upper:
            raise ValueError("Spectral lower wavelength must be less than upper wavelength.")
        spectra.append((lower, upper))
    return spectra


def _py6s_output_nodata(source_nodata, override, dtype, scale):
    """Keep output invalid pixels distinct from valid, possibly zero reflectance."""
    integer = np.issubdtype(dtype, np.integer)
    limits = np.iinfo(dtype) if integer else np.finfo(dtype)
    candidate = override if override is not None else source_nodata
    if candidate is not None:
        candidate = float(candidate)
        representable = (np.isnan(candidate) and not integer) or (
            np.isfinite(candidate) and limits.min <= candidate <= limits.max
            and (not integer or candidate.is_integer()))
        if override is not None and not representable:
            raise ValueError(f"custom_nodata_value is not representable in {dtype.name}.")
        # Inherited DN nodata (often zero) can be a valid reflectance value.
        collision = np.isfinite(candidate) and 0 <= candidate <= scale
        if representable and (override is not None or not collision):
            return candidate
    if not integer:
        return float("nan")
    for candidate in (limits.min, limits.max):
        if not 0 <= candidate <= scale:
            return int(candidate)
    raise ValueError("No integer NoData value remains outside the output range; use a wider dtype.")


@_logged_operation("py6s_atmospheric_correction", inputs=("input_raster",), outputs=("output_raster",), worker_progress=True)
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
    wavelengths_center_um: float | list[float] | None = None,
    wavelengths_bounds_um: list[float] | list[list[float]] | None = None,
    wavelengths_response_csv: str | None = None,
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
        wavelengths_center_um: Band centers in micrometers, in raster order. Choose exactly one spectral parameter.
        wavelengths_bounds_um: Per-band [lower, upper] bounds in micrometers with uniform response.
        wavelengths_response_csv: CSV path with wavelength_nm and response columns in raster order; resampled to the 6S 2.5 nm grid at runtime.
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
        custom_nodata_value: Output nodata override. None inherits a usable source value, otherwise selects NaN for float outputs or an integer sentinel outside the scaled reflectance range.
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
        "wavelengths_center_um",
        "wavelengths_bounds_um",
        "wavelengths_response_csv",
        "dn_to_radiance_factors",
        "dn_to_radiance_offsets",
    ):
        if effective[key] is None:
            effective[key] = (variables or {}).get(key)
    _py6s_spectral_parameter(effective)
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
                    crop=True,
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

        spectra = _py6s_spectral_bands(kwargs)
        input_scale_factor = float(kwargs.get("input_scale_factor", 1.0))
        dn_to_radiance_factors = kwargs.get("dn_to_radiance_factors")
        dn_to_radiance_offsets = kwargs.get("dn_to_radiance_offsets")
        output_scale_factor = kwargs.get("output_scale_factor")
        output_dtype = str(kwargs.get("output_dtype", "float32"))
        clip_reflectance = bool(kwargs.get("clip_reflectance", True))
        nodata_override = kwargs.get("custom_nodata_value")
        visibility_km = kwargs.get("visibility_km")
        sixs_executable = kwargs.get("sixs_executable")

        dtype = np.dtype(output_dtype)
        if dtype.kind not in "ifu" or not rasterio.dtypes.check_dtype(output_dtype):
            raise ValueError("output_dtype must be a real floating-point or integer raster dtype.")
        scale = 1.0 if output_scale_factor is None else float(output_scale_factor)
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError("output_scale_factor must be finite and positive.")
        if dtype.kind in "iu" and (output_scale_factor is None or scale <= 1):
            raise ValueError("Integer reflectance output requires output_scale_factor greater than 1; use float32 for unscaled reflectance.")
        scalars = dict(solar_zenith=solar_zenith, solar_azimuth=solar_azimuth,
                       view_zenith=view_zenith, view_azimuth=view_azimuth,
                       ground_elevation_km=ground_elevation_km, aot550=aot550,
                       water_vapor=water_vapor, ozone=ozone, input_scale_factor=input_scale_factor)
        if not all(np.isfinite(value) for value in scalars.values()):
            raise ValueError("Py6S geometry, atmosphere, elevation and input scaling must be finite.")
        if not 0 <= solar_zenith < 90 or not 0 <= view_zenith < 90:
            raise ValueError("Solar and view zenith must be in [0, 90) degrees.")
        if not 1 <= month <= 12 or not 1 <= day <= 31:
            raise ValueError("Py6S acquisition month/day are outside their valid ranges.")
        if min(aot550, water_vapor, ozone) < 0:
            raise ValueError("AOT, water vapor and ozone must be nonnegative.")
        if visibility_km is not None and (not np.isfinite(float(visibility_km)) or float(visibility_km) <= 0):
            raise ValueError("visibility_km must be finite and positive.")

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
            if src.count != len(spectra):
                raise ValueError(
                    f"Spectral definitions must match the source band count exactly. bands={src.count}, "
                    f"spectral_bands={len(spectra)}"
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

            for name, values in (("dn_to_radiance_factors", dn_to_radiance_factors),
                                 ("dn_to_radiance_offsets", dn_to_radiance_offsets)):
                if values is not None and not np.isfinite(np.asarray(values, dtype=float)).all():
                    raise ValueError(f"{name} must contain finite values.")
            nodata = _py6s_output_nodata(src.nodata, nodata_override, dtype, scale)
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
                s.wavelength = Wavelength(*spectra[band_idx])
                s.run()
                coefficients = (float(s.outputs.coef_xa), float(s.outputs.coef_xb), float(s.outputs.coef_xc))
                if not np.isfinite(coefficients).all():
                    raise ValueError(f"6S produced nonfinite coefficients for band {band_idx + 1}.")
                coeffs.append(coefficients)

            os.makedirs(os.path.dirname(output_raster) or ".", exist_ok=True)
            with rasterio.open(output_raster, "w", **profile) as dst:
                if input_tags:
                    dst.update_tags(**input_tags)
                if input_rpc_tags:
                    dst.update_tags(ns="RPC", **input_rpc_tags)
                if input_gcps:
                    dst.gcps = (input_gcps, input_gcps_crs)
                dst.descriptions = src.descriptions
                dst.scales = (1.0 / scale,) * src.count
                dst.offsets = (0.0,) * src.count
                dst.update_tags(REFLECTANCE_SCALE_FACTOR=scale, ATMOSPHERIC_CORRECTION="Py6S")

                for band_idx in range(1, src.count + 1):
                    xa, xb, xc = coeffs[band_idx - 1]
                    dst.set_band_unit(band_idx, "reflectance")
                    dst.update_tags(band_idx, PY6S_XA=xa, PY6S_XB=xb, PY6S_XC=xc)
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
                        mask = (src.read_masks(band_idx, window=window) > 0) & np.isfinite(in_block)
                        # Output NoData never determines which input DN values are valid.
                        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                            radiance = in_block * band_scale + band_offset
                            y = xa * radiance - xb
                            denominator = 1.0 + xc * y
                            corrected = y / denominator
                        mask &= np.isfinite(corrected) & (denominator > 0)
                        if clip_reflectance:
                            corrected = np.clip(corrected, 0.0, 1.0)
                        corrected = corrected * scale
                        if dtype.kind in "iu":
                            corrected = np.rint(corrected)
                        limits = np.iinfo(dtype) if dtype.kind in "iu" else np.finfo(dtype)
                        valid_values = corrected[mask]
                        if valid_values.size and (np.any(valid_values < limits.min) or np.any(valid_values > limits.max)):
                            raise ValueError(f"Corrected band {band_idx} exceeds {dtype.name}; use a wider dtype or a smaller output scale.")
                        corrected = np.where(mask, corrected, nodata)
                        dst.write(corrected.astype(dtype), band_idx, window=window)

        return output_raster


class Py6SAtmosphericCorrection(FunctionPlugin):
    target = "vhrharmonize.plugins.py6s_atmospheric_correction:run_py6s"
    input_dependency_paths = frozenset({"dem_file_path", "input_raster", "wavelengths_response_csv"})
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset({"output_raster"})
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


__all__ = ["Py6SAtmosphericCorrection", "Py6SRunResult", "Py6SCorrector", "run_py6s"]
