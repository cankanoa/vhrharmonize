"""Atmospheric correction step module."""

from concurrent.futures import ProcessPoolExecutor, as_completed
import itertools
import os
import shutil
from typing import Any, Dict, List, Protocol

import numpy as np
import rasterio
from tqdm import tqdm


class AtmosphericCorrector(Protocol):
    """Atmospheric correction adapter interface."""

    def run(self, input_raster: str, output_raster: str, **kwargs: Any) -> str:
        """Run correction and return output path."""
        ...


def run_flaash(
    flaash_params,
    output_params_path,
    envi_engine,
    output_image_path_to_delete=None
    ):
    """Execute ENVI FLAASH with the provided parameter dictionary."""

    print('Processing photo with params:', flaash_params)
    if output_image_path_to_delete:
        try:
            if os.path.exists(output_image_path_to_delete):
                os.remove(output_image_path_to_delete)
                print(f"Deleted existing output image: {output_image_path_to_delete}")
        except Exception as e:
            print(f"Error deleting output image {output_image_path_to_delete}: {e}")

    try:
        task = envi_engine.task('FLAASH')
        task.execute(flaash_params)
        with open(output_params_path, 'w') as file:
            file.write(str(flaash_params))
        print("Photo complete")

    except Exception as e:
        print(f"Error processing: {flaash_params}: {e}")


def run_flaash_wrapper(
    args
    ):

    """
    Wrapper function (needed because executors only picklable callables).
    Returns the 'OUTPUT_RASTER_URI' so we can collect it later.
    """
    test_params, test_output_params_path, envi_engine = args
    run_flaash(test_params, test_output_params_path, envi_engine)
    return test_params["OUTPUT_RASTER_URI"]


def parallel_flaash(
    test_flaash_params_array,
    envi_engine,
    max_workers=4
    ):

    """
    Executes run_flaash in parallel on the items in test_flaash_params_array.
    """
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


def create_test_flaash_params(
    flaash_params,
    output_params_path
    ):
    """Expand selected FLAASH parameter grids into per-run parameter/output tuples."""

    test_params = {
                 'SENSOR_TYPE': None,
                 'INPUT_SCALE': None,
                 'OUTPUT_SCALE': None,
                 'CALIBRATION_FILE': None,
                 'CALIBRATION_FORMAT': None,
                 'CALIBRATION_UNITS': None,
                 'LAT_LONG': None,
                 'SENSOR_ALTITUDE': None,
                 'DATE_TIME': None,
                 'USE_ADJACENCY': None,
                 'DEFAULT_VISIBILITY': None,
                 'USE_POLISHING': None,
                 'POLISHING_RESOLUTION': None,
                 'SENSOR_AUTOCALIBRATION': None,
                 'SENSOR_CAL_PRECISION': None,
                 'SENSOR_CAL_FEATURE_LIST': None,
                 'GROUND_ELEVATION': {None, 0, 1.6, 2, 3, 4, 5},
                 'SOLAR_AZIMUTH': None,
                 'SOLAR_ZENITH': None,
                 'LOS_AZIMUTH': None,
                 'LOS_ZENITH': None,
                 'IFOV': None,
                 'MODTRAN_ATM': None,
                 'MODTRAN_AER': None,
                 'MODTRAN_RES': {15},
                 'MODTRAN_MSCAT': None,
                 'CO2_MIXING': None,
                 'WATER_ABS_CHOICE': None,
                 'WATER_MULT': None,
                 'WATER_VAPOR_PRESET': None,
                 'USE_AEROSOL': {'Disabled'},
                 'AEROSOL_SCALE_HT': {2, 1, 2.5, 0.5},
                 'AER_BAND_RATIO': {0.5, 0.2, 0.9},
                 'AER_BAND_WAVL': None,
                 'AER_REFERENCE_VALUE': None,
                 'AER_REFERENCE_PIXEL': None,
                 'AER_BANDLOW_WAVL': None,
                 'AER_BANDLOW_MAXREFL': None,
                 'AER_BANDHIGH_WAVL': None,
                 'AER_BANDHIGH_MAXREFL': None,
                 }

    testable_params = {
        k: v
        for k, v in test_params.items()
        if k != "OUTPUT_RASTER_URI" and v is not None and len(v) > 0
    }

    if not testable_params:
        return [(flaash_params, output_params_path)]

    param_items = list(testable_params.items())
    param_names = [item[0] for item in param_items]
    param_values_lists = [item[1] for item in param_items]
    all_combinations = list(itertools.product(*param_values_lists))

    base_output_raster_uri = flaash_params.get("OUTPUT_RASTER_URI", "")
    base_raster_dir = os.path.dirname(base_output_raster_uri)
    base_raster_filename = os.path.basename(base_output_raster_uri)
    raster_name, raster_ext = os.path.splitext(base_raster_filename)

    base_params_dir = os.path.dirname(output_params_path)
    base_params_filename = os.path.basename(output_params_path)
    params_name, params_ext = os.path.splitext(base_params_filename)

    result_list = []

    for combo in all_combinations:
        combo_params = dict(flaash_params)
        short_suffix_parts = []
        for name, val in zip(param_names, combo):
            short_suffix_parts.append(f"{name[:3]}-{str(val)[:3]}")
            combo_params[name] = val

        short_suffix = "_".join(short_suffix_parts)

        new_raster_filename = f"{raster_name}_{short_suffix}{raster_ext}"
        new_raster_path = os.path.join(base_raster_dir, new_raster_filename)
        combo_params["OUTPUT_RASTER_URI"] = new_raster_path

        new_params_filename = f"{params_name}_{short_suffix}{params_ext}"
        new_params_path = os.path.join(base_params_dir, new_params_filename)

        result_list.append((combo_params, new_params_path))

    return result_list


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
        run_flaash(params, params_path, self.envi_engine, output_raster)
        return output_raster


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


class Py6SCorrector:
    """Atmospheric correction adapter using Py6S.

    This is a scaffold: wire sensor response, geometry, aerosol/atmosphere models,
    and per-band correction strategy before production use.
    """

    def run(self, input_raster: str, output_raster: str, **kwargs: Any) -> str:
        """Run block-wise Py6S atmospheric correction for multispectral rasters.

        Required kwargs:
            solar_zenith, solar_azimuth, view_zenith, view_azimuth, day, month

        Optional kwargs:
            ground_elevation_km (float, default 0.0)
            atmosphere_profile (str, default "midlatitude_summer")
            aerosol_profile (str, default "maritime")
            aot550 (float, default 0.2)
            water_vapor (float, default 2.5)
            ozone (float, default 0.3)
            band_wavelengths_um (list[float], default WV3 multispectral centers)
            input_scale_factor (float, default 1.0)
            output_scale_factor (float | None, default None)
            output_dtype (str, default "float32")
            clip_reflectance (bool, default True)
            nodata_value (float | None, default source nodata)
            dn_to_radiance_factors (list[float] | None, default None)
            dn_to_radiance_offsets (list[float] | None, default None)
        """
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
    "AtmosphericCorrector",
    "run_flaash",
    "run_flaash_wrapper",
    "parallel_flaash",
    "create_test_flaash_params",
    "FLAASHCorrector",
    "Py6SCorrector",
    "atmospheric_correction",
]
