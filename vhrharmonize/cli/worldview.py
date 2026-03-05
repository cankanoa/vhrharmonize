#!/usr/bin/env python3
"""CLI wrapper for the WorldView preprocessing workflow."""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from typing import Dict, Iterable, List, Optional


def _parse_filter_basenames(raw_values: Optional[List[str]]) -> List[str]:
    if not raw_values:
        return []
    parsed: List[str] = []
    for raw in raw_values:
        for value in raw.split(","):
            value = value.strip()
            if value:
                parsed.append(value)
    return parsed


def _init_envi_engine(envi_engine_path: str):
    import envipyengine.config
    from envipyengine import Engine

    envipyengine.config.set("engine", envi_engine_path)
    envi_engine = Engine("ENVI")
    envi_engine.tasks()
    return envi_engine


def _build_flaash_params(
    mul_tif_file: str,
    dem_file_path: str,
    mul_gpkg_path: str,
    mul_imd_data: Dict[str, float],
    mul_flaash_image_path: str,
    *,
    dem_ground_percentile: float,
    modtran_atm: str,
    modtran_aer: str,
    use_aerosol: str,
    default_visibility: Optional[float],
) -> Dict:
    from vhrharmonize.io.geospatial import get_image_percentile_value

    ground_elevation_m = get_image_percentile_value(
        dem_file_path,
        percentile=dem_ground_percentile,
        mask=mul_gpkg_path,
    )

    return {
        "INPUT_RASTER": {"url": mul_tif_file, "factory": "URLRaster"},
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
        "GROUND_ELEVATION": ground_elevation_m / 1000,
        "SOLAR_AZIMUTH": mul_imd_data.get("SOLAR_AZIMUTH"),
        "SOLAR_ZENITH": mul_imd_data.get("SOLAR_ZENITH"),
        "LOS_AZIMUTH": mul_imd_data.get("LOS_AZIMUTH"),
        "LOS_ZENITH": mul_imd_data.get("LOS_ZENITH"),
        "OUTPUT_RASTER_URI": mul_flaash_image_path,
    }


def _parse_worldview_datetime_from_basename(photo_basename: str) -> datetime:
    """Parse WorldView acquisition datetime from basename prefix."""
    month_map = {
        "JAN": 1,
        "FEB": 2,
        "MAR": 3,
        "APR": 4,
        "MAY": 5,
        "JUN": 6,
        "JUL": 7,
        "AUG": 8,
        "SEP": 9,
        "OCT": 10,
        "NOV": 11,
        "DEC": 12,
    }
    # Example: 17OCT19211717-M1BS-...
    m = re.match(r"^(\d{2})([A-Z]{3})(\d{2})(\d{2})(\d{2})(\d{2})-", photo_basename)
    if not m:
        raise ValueError(f"Could not parse acquisition datetime from basename: {photo_basename}")
    return datetime(
        year=2000 + int(m.group(3)),
        month=month_map[m.group(2)],
        day=int(m.group(1)),
        hour=int(m.group(4)),
        minute=int(m.group(5)),
        second=int(m.group(6)),
    )


def _resolve_scene_datetime(mul_imd_data: Dict, mul_photo_basename: str) -> datetime:
    raw = mul_imd_data.get("ACQUISITION_DATETIME_UTC")
    if isinstance(raw, str) and raw.strip():
        try:
            return datetime.fromisoformat(raw.replace("Z", "+00:00"))
        except ValueError:
            pass
    return _parse_worldview_datetime_from_basename(mul_photo_basename)


def _build_py6s_kwargs(
    mul_imd_data: Dict[str, float],
    mul_photo_basename: str,
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
) -> Dict:
    dt = _resolve_scene_datetime(mul_imd_data, mul_photo_basename)
    py6s_kwargs = {
        "solar_zenith": mul_imd_data.get("SOLAR_ZENITH"),
        "solar_azimuth": mul_imd_data.get("SOLAR_AZIMUTH"),
        "view_zenith": mul_imd_data.get("VIEW_ZENITH", mul_imd_data.get("LOS_ZENITH")),
        "view_azimuth": mul_imd_data.get("LOS_AZIMUTH"),
        "day": dt.day,
        "month": dt.month,
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
    dn_to_radiance_factors = mul_imd_data.get("worldview_dn_to_radiance_factors")
    if use_imd_radiance_calibration and dn_to_radiance_factors:
        if use_worldview_gain_offset_adjustment:
            gains = mul_imd_data.get("worldview_dn_to_radiance_gains")
            if gains and len(gains) >= len(dn_to_radiance_factors):
                dn_to_radiance_factors = [
                    float(dn_to_radiance_factors[idx]) * float(gains[idx])
                    for idx in range(len(dn_to_radiance_factors))
                ]
        py6s_kwargs["dn_to_radiance_factors"] = dn_to_radiance_factors
        dn_to_radiance_offsets = mul_imd_data.get("worldview_dn_to_radiance_offsets")
        if use_worldview_gain_offset_adjustment and dn_to_radiance_offsets:
            py6s_kwargs["dn_to_radiance_offsets"] = dn_to_radiance_offsets
    return py6s_kwargs


def _scene_bbox_wgs84_from_shp(shp_path: str) -> tuple[float, float, float, float]:
    import geopandas as gpd

    gdf = gpd.read_file(shp_path)
    if gdf.empty:
        raise ValueError(f"Empty scene footprint shapefile: {shp_path}")
    if gdf.crs is None:
        raise ValueError(f"Scene footprint shapefile has no CRS: {shp_path}")
    gdf = gdf.to_crs(epsg=4326)
    minx, miny, maxx, maxy = gdf.total_bounds
    return float(minx), float(miny), float(maxx), float(maxy)


def _resolve_output_resolution_for_epsg(output_epsg: int, product_res: Optional[float]) -> Optional[float]:
    """Return safe output resolution for target CRS.

    WorldView IMD `product_res` is in meters. If output EPSG is geographic
    (for example 4326, degrees), passing meter resolution directly would produce
    invalid/extremely coarse output. In that case, return None so GDAL chooses
    default resolution from source georeferencing.
    """
    if product_res is None:
        return None
    try:
        import pyproj

        crs = pyproj.CRS.from_epsg(int(output_epsg))
        if crs.is_geographic:
            print(
                f"Output EPSG:{output_epsg} is geographic; ignoring meter-based product_res={product_res} "
                "for orthorectification."
            )
            return None
    except Exception:
        # If CRS lookup fails, preserve previous behavior.
        pass
    return product_res


def _wsl_path_to_windows_for_envi(path: str) -> str:
    """
    Convert /mnt/<drive>/... WSL path into <Drive>:\\... for ENVI on Windows.
    """
    match = re.match(r"^/mnt/([a-zA-Z])/(.*)$", path)
    if not match:
        raise ValueError(
            f"Path is not Windows-drive-backed via /mnt/<drive>/: {path}"
        )
    drive = match.group(1).upper()
    rest = match.group(2).replace("/", "\\")
    return f"{drive}:\\{rest}"


def _convert_flaash_params_paths_for_windows(flaash_params: Dict) -> Dict:
    converted = dict(flaash_params)
    input_raster = converted.get("INPUT_RASTER")
    if isinstance(input_raster, dict) and "url" in input_raster:
        input_raster = dict(input_raster)
        input_raster["url"] = _wsl_path_to_windows_for_envi(input_raster["url"])
        converted["INPUT_RASTER"] = input_raster

    for key in ("OUTPUT_RASTER_URI", "CLOUD_RASTER_URI", "WATER_RASTER_URI"):
        if converted.get(key):
            converted[key] = _wsl_path_to_windows_for_envi(converted[key])

    return converted


def _required_keys_present(found_default_file: Dict, required_keys: Iterable[str]) -> bool:
    return all(found_default_file.get(key) for key in required_keys)


def _build_cloud_mask_output_path(input_image_path: str, suffix: str) -> str:
    base, ext = os.path.splitext(input_image_path)
    return f"{base}{suffix}{ext}"


def _parse_int_csv(raw_values: str) -> List[int]:
    parsed = []
    for value in raw_values.split(","):
        value = value.strip()
        if value:
            parsed.append(int(value))
    return parsed


def _parse_json_dict(raw_json: Optional[object]) -> Dict:
    if raw_json is None or raw_json == "":
        return {}
    if isinstance(raw_json, dict):
        return raw_json
    if not isinstance(raw_json, str):
        raise ValueError(
            "Expected omnicloud kwargs as a JSON string or dictionary mapping."
        )
    parsed = json.loads(raw_json)
    if not isinstance(parsed, dict):
        raise ValueError("Expected a JSON object for omnicloud kwargs.")
    return parsed


def _load_yaml_config(config_yaml_path: str) -> Dict:
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError(
            "PyYAML is required for --config-yaml. Install it with `pip install pyyaml`."
        ) from exc

    with open(config_yaml_path, "r", encoding="utf-8") as f:
        loaded = yaml.safe_load(f) or {}
    if not isinstance(loaded, dict):
        raise ValueError("Config YAML root must be a mapping/dictionary.")

    # Support either arg-dest style (input_dir) or CLI style (input-dir).
    normalized = {}
    for key, value in loaded.items():
        normalized[key.replace("-", "_")] = value
    return normalized


def _normalize_config_defaults(config_defaults: Dict) -> Dict:
    normalized = dict(config_defaults)
    for list_key in ("input_dir", "filter_basename"):
        if list_key in normalized and isinstance(normalized[list_key], str):
            normalized[list_key] = [normalized[list_key]]
    return normalized


def _run_cloud_mask_command(
    command_template: str,
    input_image_path: str,
    output_image_path: str,
    scene_root_path: str,
    image_basename: str,
) -> None:
    command = command_template.format(
        input=input_image_path,
        output=output_image_path,
        scene_root=scene_root_path,
        image_basename=image_basename,
    )
    print(f"Running cloud mask command: {command}")
    subprocess.run(command, shell=True, check=True)


def _write_scene_metadata_report(report_path: str, payload: Dict) -> None:
    os.makedirs(os.path.dirname(report_path) or ".", exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


def _count_mask_pixels(mask_path: str, mask_value: int = 1) -> int:
    import rasterio

    count = 0
    with rasterio.open(mask_path) as src:
        for _, window in src.block_windows(1):
            block = src.read(1, window=window)
            count += int((block == mask_value).sum())
    return count


def run_workflow(args: argparse.Namespace) -> int:
    """Run full-scene WorldView preprocessing for one or more input directories."""
    from vhrharmonize.providers.worldview.files import (
        find_files,
        find_roots,
        get_metadata_from_files,
    )
    from vhrharmonize.preprocess.atmospheric_correction import atmospheric_correction, run_flaash
    from vhrharmonize.preprocess.cloudmasking import (
        apply_binary_cloud_mask_to_image,
        create_cloud_mask_with_omnicloudmask,
    )
    from vhrharmonize.preprocess.orthorectification import gcp_refined_rpc_orthorectification
    from vhrharmonize.preprocess.pansharpening import pansharpen_image
    from vhrharmonize.io.geospatial import shp_to_gpkg

    filter_basenames = _parse_filter_basenames(args.filter_basename)
    use_existing_mul_ortho = bool(args.existing_mul_ortho_input)
    use_existing_pan_ortho = bool(args.existing_pan_ortho_input)
    use_existing_ortho_inputs = bool(use_existing_mul_ortho and use_existing_pan_ortho)
    cloud_classes = _parse_int_csv(args.cloud_mask_classes)
    omnicloud_kwargs = _parse_json_dict(args.cloud_mask_omnicloud_kwargs_json)
    os.makedirs(args.output_dir, exist_ok=True)

    envi_engine = None
    if (
        args.run_atmospheric_correction
        and args.atmospheric_method == "flaash"
        and not args.skip_flaash
        and not use_existing_ortho_inputs
    ):
        envi_engine = _init_envi_engine(args.envi_engine_path)

    for input_folder in args.input_dir:
        print(f"Scanning input folder: {input_folder}")
        root_paths = find_roots(input_folder)

        for root_folder_path, root_file_path in root_paths:
            found_default_files = find_files(root_folder_path, root_file_path, filter_basenames)

            for _, found_default_file in found_default_files.items():
                mul_photo_basename = found_default_file.get("mul_photo_basename")
                pan_photo_basename = found_default_file.get("pan_photo_basename")
                print(f"Processing: {mul_photo_basename or pan_photo_basename}")
                scene_started_utc = datetime.utcnow().isoformat() + "Z"

                required = (
                    "mul_imd_file",
                    "mul_tif_file",
                    "pan_imd_file",
                    "pan_tif_file",
                    "mul_shp_file",
                )
                if not _required_keys_present(found_default_file, required):
                    print("Skipping scene: required files missing for workflow.")
                    continue

                mul_imd_file = found_default_file["mul_imd_file"]
                pan_imd_file = found_default_file["pan_imd_file"]
                mul_tif_file = found_default_file["mul_tif_file"]
                pan_tif_file = found_default_file["pan_tif_file"]
                mul_shp_path = found_default_file["mul_shp_file"]

                params_overrides_scene, mul_params_overrides_photo, mul_imd_data = (
                    get_metadata_from_files(root_file_path, mul_imd_file, mul_photo_basename)
                )
                _, _, pan_imd_data = get_metadata_from_files(
                    root_file_path, pan_imd_file, pan_photo_basename
                )

                with tempfile.TemporaryDirectory(dir=args.scratch_dir, prefix="vhr_scene_") as work_dir:
                    py6s_effective_params = None
                    py6s_auto_atmos_estimate = None
                    mask_output_path = None
                    cloud_mask_pixel_count = None
                    masked_image_output_path = None
                    ortho_folder_name = "OrthoFromDefaultRPC"
                    if args.run_atmospheric_correction:
                        if args.atmospheric_method == "py6s":
                            atmos_label = "PY6S"
                        elif args.atmospheric_method == "flaash":
                            atmos_label = "FLAASH"
                        else:
                            atmos_label = "RAW"
                    else:
                        atmos_label = "RAW"
                    if use_existing_mul_ortho:
                        print("Using existing orthorectified multispectral input")
                        mul_flaash_ortho_path = args.existing_mul_ortho_input
                        if args.run_pansharpen:
                            print("Using existing orthorectified panchromatic input")
                            pan_ortho_path = args.existing_pan_ortho_input
                    else:
                        print("Creating footprint GeoPackage for atmospheric terrain statistics")
                        mul_gpkg_path = os.path.join(work_dir, f"{mul_photo_basename}.gpkg")
                        shp_to_gpkg(mul_shp_path, mul_gpkg_path, args.footprint_epsg)

                        from vhrharmonize.io.geospatial import get_image_percentile_value

                        ground_elevation_m = get_image_percentile_value(
                            args.dem_file_path,
                            percentile=args.flaash_dem_ground_percentile,
                            mask=mul_gpkg_path,
                        )
                        ground_elevation_km = ground_elevation_m / 1000.0

                        mul_flaash_input_path = mul_tif_file
                        if args.run_atmospheric_correction:
                            if args.atmospheric_method == "none":
                                mul_flaash_input_path = mul_tif_file
                            elif args.atmospheric_method == "flaash" and args.skip_flaash:
                                mul_flaash_input_path = args.existing_flaash_input
                            elif args.atmospheric_method == "flaash":
                                mul_flaash_image_path = os.path.join(work_dir, f"{mul_photo_basename}_FLAASH.dat")
                                mul_flaash_params_path = os.path.join(work_dir, f"{mul_photo_basename}_FLAASH_Params.txt")
                                print("Running FLAASH")
                                flaash_params = _build_flaash_params(
                                    mul_tif_file,
                                    args.dem_file_path,
                                    mul_gpkg_path,
                                    mul_imd_data,
                                    mul_flaash_image_path,
                                    dem_ground_percentile=args.flaash_dem_ground_percentile,
                                    modtran_atm=args.flaash_modtran_atm,
                                    modtran_aer=args.flaash_modtran_aer,
                                    use_aerosol=args.flaash_use_aerosol,
                                    default_visibility=args.flaash_default_visibility,
                                )
                                if params_overrides_scene is not None:
                                    flaash_params.update(params_overrides_scene)
                                if mul_params_overrides_photo is not None:
                                    flaash_params.update(mul_params_overrides_photo)
                                flaash_params = {k: v for k, v in flaash_params.items() if v is not None}
                                flaash_params_for_envi = _convert_flaash_params_paths_for_windows(flaash_params)
                                run_flaash(
                                    flaash_params_for_envi,
                                    mul_flaash_params_path,
                                    envi_engine,
                                    mul_flaash_image_path,
                                )
                                if not os.path.exists(mul_flaash_image_path):
                                    print(
                                        "Skipping scene: FLAASH did not produce output. "
                                        "Check ENVI path handling and FLAASH parameters."
                                    )
                                    continue
                                mul_flaash_input_path = mul_flaash_image_path
                            else:
                                print("Running Py6S atmospheric correction")
                                mul_py6s_image_path = os.path.join(work_dir, f"{mul_photo_basename}_PY6S.tif")
                                auto_aot550 = args.py6s_aot550
                                auto_water_vapor = args.py6s_water_vapor
                                auto_ozone = args.py6s_ozone
                                if (
                                    args.py6s_auto_atmos_source == "nasa_power"
                                    and args.py6s_atmosphere_profile.strip().lower() == "user"
                                ):
                                    try:
                                        from vhrharmonize.preprocess.atmosphere_nasa import (
                                            fetch_power_atmosphere_for_bbox,
                                        )

                                        scene_dt = _resolve_scene_datetime(mul_imd_data, mul_photo_basename)
                                        min_lon, min_lat, max_lon, max_lat = _scene_bbox_wgs84_from_shp(mul_shp_path)
                                        estimate = fetch_power_atmosphere_for_bbox(
                                            day_utc=scene_dt.date(),
                                            min_lon=min_lon,
                                            min_lat=min_lat,
                                            max_lon=max_lon,
                                            max_lat=max_lat,
                                            grid_size=args.py6s_auto_atmos_grid_size,
                                            search_days=args.py6s_auto_atmos_search_days,
                                            timeout_s=args.py6s_auto_atmos_timeout_s,
                                            endpoint=args.py6s_auto_atmos_power_endpoint,
                                        )
                                        if estimate.aot550 is not None:
                                            auto_aot550 = float(estimate.aot550)
                                        if estimate.water_vapor is not None:
                                            auto_water_vapor = float(estimate.water_vapor)
                                        if estimate.ozone_cm_atm is not None:
                                            auto_ozone = float(estimate.ozone_cm_atm)
                                        py6s_auto_atmos_estimate = {
                                            "source": estimate.source,
                                            "date_used": estimate.date_used,
                                            "sample_count": estimate.sample_count,
                                            "aot550": estimate.aot550,
                                            "water_vapor": estimate.water_vapor,
                                            "ozone_cm_atm": estimate.ozone_cm_atm,
                                        }
                                        print(
                                            "Auto-updated Py6S atmosphere from NASA POWER "
                                            f"(date={estimate.date_used}, samples={estimate.sample_count}): "
                                            f"aot550={auto_aot550:.4f}, water_vapor={auto_water_vapor:.4f}, ozone={auto_ozone:.4f}"
                                        )
                                    except Exception as exc:
                                        print(
                                            "Warning: failed to auto-fetch atmosphere from NASA POWER; "
                                            f"using configured values. ({exc})"
                                        )
                                        py6s_auto_atmos_estimate = {"error": str(exc)}
                                py6s_kwargs = _build_py6s_kwargs(
                                    mul_imd_data,
                                    mul_photo_basename,
                                    ground_elevation_km=ground_elevation_km,
                                    atmosphere_profile=args.py6s_atmosphere_profile,
                                    aerosol_profile=args.py6s_aerosol_profile,
                                    aot550=auto_aot550,
                                    visibility_km=args.py6s_visibility,
                                    water_vapor=auto_water_vapor,
                                    ozone=auto_ozone,
                                    sixs_executable=args.py6s_executable,
                                    output_scale_factor=args.py6s_output_scale_factor,
                                    output_dtype=args.py6s_output_dtype,
                                    use_imd_radiance_calibration=args.py6s_use_imd_radiance_calibration,
                                    use_worldview_gain_offset_adjustment=args.py6s_use_worldview_gain_offset_adjustment,
                                )
                                py6s_effective_params = {
                                    "atmosphere_profile": args.py6s_atmosphere_profile,
                                    "aerosol_profile": args.py6s_aerosol_profile,
                                    "aot550": auto_aot550,
                                    "visibility_km": args.py6s_visibility,
                                    "water_vapor": auto_water_vapor,
                                    "ozone": auto_ozone,
                                    "output_scale_factor": args.py6s_output_scale_factor,
                                    "output_dtype": args.py6s_output_dtype,
                                    "use_imd_radiance_calibration": args.py6s_use_imd_radiance_calibration,
                                    "use_worldview_gain_offset_adjustment": args.py6s_use_worldview_gain_offset_adjustment,
                                    "auto_atmos_source": args.py6s_auto_atmos_source,
                                }
                                if "dn_to_radiance_factors" in py6s_kwargs:
                                    sensor_id = mul_imd_data.get("sensor_id")
                                    version = mul_imd_data.get("worldview_radiometric_adjustment_version")
                                    if version and args.py6s_use_worldview_gain_offset_adjustment:
                                        print(
                                            "Using WorldView IMD radiometric calibration "
                                            f"(absCal/EBW + gain/offset {version}, sensor={sensor_id})"
                                        )
                                    else:
                                        print(
                                            "Using WorldView IMD radiometric calibration "
                                            "(absCal/EBW only; no gain/offset table available)"
                                        )
                                atmospheric_correction(
                                    mul_tif_file,
                                    mul_py6s_image_path,
                                    method="py6s",
                                    **py6s_kwargs,
                                )
                                mul_flaash_input_path = mul_py6s_image_path
                        else:
                            print("Skipping atmospheric correction")

                        print("Using default RPC model for orthorectification (no GCP refinement)")
                        print("Orthorectifying multispectral image")
                        mul_flaash_ortho_path = os.path.join(
                            work_dir,
                            f"{mul_photo_basename}_{atmos_label}_{ortho_folder_name}.tif",
                        )
                        gcp_refined_rpc_orthorectification(
                            mul_flaash_input_path,
                            mul_flaash_ortho_path,
                            args.dem_file_path,
                            args.epsg,
                            output_nodata_value=args.nodata_value,
                            dtype=args.dtype,
                            output_resolution=_resolve_output_resolution_for_epsg(
                                args.epsg,
                                mul_imd_data.get("product_res"),
                            ),
                        )

                        if args.run_pansharpen:
                            print("Orthorectifying panchromatic image")
                            pan_ortho_path = os.path.join(
                                work_dir,
                                f"{pan_photo_basename}_{ortho_folder_name}.tif",
                            )
                            gcp_refined_rpc_orthorectification(
                                pan_tif_file,
                                pan_ortho_path,
                                args.dem_file_path,
                                args.epsg,
                                output_nodata_value=args.nodata_value,
                                dtype=args.dtype,
                                output_resolution=_resolve_output_resolution_for_epsg(
                                    args.epsg,
                                    pan_imd_data.get("product_res"),
                                ),
                            )

                    if not args.run_pansharpen:
                        print("Skipping pansharpen step")
                        final_scene_path = mul_flaash_ortho_path
                    else:
                        print("Pansharpening multispectral image")
                        mul_pansharp_path = os.path.join(
                            work_dir,
                            f"{mul_photo_basename}_{atmos_label}_{ortho_folder_name}_Pansharps.tif",
                        )
                        pansharpen_image(
                            mul_flaash_ortho_path,
                            pan_ortho_path,
                            mul_pansharp_path,
                            change_nodata_value=args.nodata_value,
                        )
                        final_scene_path = mul_pansharp_path

                    print("Writing final scene output")

                    if args.run_cloud_mask and args.cloud_mask_command:
                        cloud_mask_output_path = _build_cloud_mask_output_path(
                            final_scene_path,
                            args.cloud_mask_output_suffix,
                        )
                        _run_cloud_mask_command(
                            args.cloud_mask_command,
                            final_scene_path,
                            cloud_mask_output_path,
                            root_folder_path,
                            mul_photo_basename,
                        )
                        final_scene_path = cloud_mask_output_path

                    if args.run_cloud_mask and args.cloud_mask_method == "omnicloudmask":
                        cloud_mask_input_path = mul_flaash_ortho_path
                        mask_output_path = _build_cloud_mask_output_path(
                            final_scene_path,
                            args.cloud_mask_mask_suffix,
                        )
                        masked_image_output_path = _build_cloud_mask_output_path(
                            final_scene_path,
                            args.cloud_mask_output_suffix,
                        )
                        print("Running built-in OmniCloudMask cloud masking on orthorectified multispectral image")
                        create_cloud_mask_with_omnicloudmask(
                            cloud_mask_input_path,
                            mask_output_path,
                            red_band_index=args.cloud_mask_red_band_index,
                            green_band_index=args.cloud_mask_green_band_index,
                            nir_band_index=args.cloud_mask_nir_band_index,
                            cloud_classes=cloud_classes,
                            buffer_pixels=args.cloud_buffer_pixels,
                            omnicloud_kwargs=omnicloud_kwargs,
                            inference_resolution_m=args.cloud_mask_inference_resolution_m,
                        )
                        cloud_mask_pixel_count = _count_mask_pixels(mask_output_path, mask_value=1)
                        print(f"Cloud mask pixels (value=1): {cloud_mask_pixel_count}")
                        apply_binary_cloud_mask_to_image(
                            final_scene_path,
                            mask_output_path,
                            masked_image_output_path,
                            output_nodata_value=args.nodata_value,
                            allow_mask_reprojection=True,
                        )
                        final_scene_path = masked_image_output_path
                    elif not args.run_cloud_mask:
                        print("Skipping cloud masking")

                    scene_output_filename = f"{mul_photo_basename}{args.output_suffix}.tif"
                    scene_output_path = os.path.join(args.output_dir, scene_output_filename)
                    shutil.copy2(final_scene_path, scene_output_path)
                    print(f"Wrote: {scene_output_path}")

                    scene_metadata_path = os.path.join(
                        args.output_dir,
                        f"{mul_photo_basename}{args.output_suffix}_metadata.json",
                    )
                    scene_metadata = {
                        "scene": {
                            "input_dir": input_folder,
                            "scene_root": root_folder_path,
                            "mul_photo_basename": mul_photo_basename,
                            "pan_photo_basename": pan_photo_basename,
                            "started_utc": scene_started_utc,
                            "completed_utc": datetime.utcnow().isoformat() + "Z",
                        },
                        "inputs": {
                            "mul_imd_file": mul_imd_file,
                            "mul_tif_file": mul_tif_file,
                            "pan_imd_file": pan_imd_file,
                            "pan_tif_file": pan_tif_file,
                            "mul_shp_file": mul_shp_path,
                            "dem_file_path": args.dem_file_path,
                            "root_file_path": root_file_path,
                        },
                        "workflow": {
                            "atmospheric_method": args.atmospheric_method,
                            "run_atmospheric_correction": args.run_atmospheric_correction,
                            "run_pansharpen": args.run_pansharpen,
                            "run_cloud_mask": args.run_cloud_mask,
                            "epsg": args.epsg,
                            "nodata_value": args.nodata_value,
                            "dtype": args.dtype,
                            "output_suffix": args.output_suffix,
                            "flaash_dem_ground_percentile": args.flaash_dem_ground_percentile,
                        },
                        "py6s": {
                            "configured": {
                                "atmosphere_profile": args.py6s_atmosphere_profile,
                                "aerosol_profile": args.py6s_aerosol_profile,
                                "aot550": args.py6s_aot550,
                                "visibility_km": args.py6s_visibility,
                                "water_vapor": args.py6s_water_vapor,
                                "ozone": args.py6s_ozone,
                                "output_scale_factor": args.py6s_output_scale_factor,
                                "output_dtype": args.py6s_output_dtype,
                                "use_imd_radiance_calibration": args.py6s_use_imd_radiance_calibration,
                                "use_worldview_gain_offset_adjustment": args.py6s_use_worldview_gain_offset_adjustment,
                                "auto_atmos_source": args.py6s_auto_atmos_source,
                            },
                            "effective": py6s_effective_params,
                            "auto_atmos_estimate": py6s_auto_atmos_estimate,
                        },
                        "cloud_mask": {
                            "enabled": args.run_cloud_mask,
                            "method": args.cloud_mask_method,
                            "command": args.cloud_mask_command,
                            "inference_resolution_m": args.cloud_mask_inference_resolution_m,
                            "classes": cloud_classes,
                            "buffer_pixels": args.cloud_buffer_pixels,
                            "kwargs": omnicloud_kwargs,
                            "mask_output_path": mask_output_path,
                            "mask_pixel_count": (
                                cloud_mask_pixel_count
                                if args.run_cloud_mask and args.cloud_mask_method == "omnicloudmask"
                                else None
                            ),
                            "masked_image_output_path": masked_image_output_path,
                        },
                        "outputs": {
                            "final_scene_path": scene_output_path,
                            "scene_metadata_path": scene_metadata_path,
                        },
                    }
                    _write_scene_metadata_report(scene_metadata_path, scene_metadata)
                    print(f"Wrote metadata: {scene_metadata_path}")
                print("Scene complete")

    print("All processing complete")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the WorldView preprocessing CLI."""
    parser = argparse.ArgumentParser(
        description="Run WorldView preprocessing (atmospheric correction + orthorectify + pansharpen) on full scenes.",
    )
    parser.add_argument(
        "--config-yaml",
        help="Optional YAML config file. Keys should match argument names (use '_' or '-').",
    )
    parser.add_argument(
        "--input-dir",
        action="append",
        help=(
            "Input directory to scan for scene folders. "
            "If Root*.txt files exist they are used for overrides; otherwise scene roots are auto-discovered."
        ),
    )
    parser.add_argument("--dem-file-path", help="DEM GeoTIFF path in WGS84 ellipsoidal height.")
    parser.add_argument("--envi-engine-path", help="Path to ENVI taskengine executable.")
    parser.add_argument(
        "--atmospheric-method",
        choices=["flaash", "py6s", "none"],
        default="py6s",
        help="Atmospheric correction backend for multispectral input before orthorectification.",
    )
    parser.add_argument("--epsg", type=int, default=4326, help="Output EPSG code for orthorectification.")
    parser.add_argument("--nodata-value", type=float, default=-9999, help="Output NoData value.")
    parser.add_argument("--dtype", default="int16", help="GDAL dtype used during orthorectification.")
    parser.add_argument(
        "--flaash-dem-ground-percentile",
        type=float,
        default=50.0,
        help="DEM percentile used for FLAASH GROUND_ELEVATION (50=median).",
    )
    parser.add_argument(
        "--flaash-modtran-atm",
        default="Mid-Latitude Summer",
        help="Base MODTRAN_ATM value for FLAASH (can be overridden per scene/photo in Root_WV.txt).",
    )
    parser.add_argument(
        "--flaash-modtran-aer",
        default="Maritime",
        help="Base MODTRAN_AER value for FLAASH (can be overridden per scene/photo in Root_WV.txt).",
    )
    parser.add_argument(
        "--flaash-use-aerosol",
        default="Disabled",
        help="Base USE_AEROSOL for FLAASH (e.g., Disabled or Automatic Selection).",
    )
    parser.add_argument(
        "--flaash-default-visibility",
        type=float,
        help="Optional base DEFAULT_VISIBILITY for FLAASH in km.",
    )
    parser.add_argument(
        "--py6s-atmosphere-profile",
        default="midlatitude_summer",
        help="Py6S atmosphere profile (e.g., midlatitude_summer, tropical, user).",
    )
    parser.add_argument(
        "--py6s-aerosol-profile",
        default="maritime",
        help="Py6S aerosol profile (e.g., maritime, continental, urban, desert).",
    )
    parser.add_argument(
        "--py6s-aot550",
        type=float,
        default=0.2,
        help="Py6S aerosol optical thickness at 550nm.",
    )
    parser.add_argument(
        "--py6s-visibility",
        type=float,
        help="Optional Py6S visibility in km (used instead of --py6s-aot550 when provided).",
    )
    parser.add_argument(
        "--py6s-water-vapor",
        type=float,
        default=2.5,
        help="Py6S water vapor (g/cm^2).",
    )
    parser.add_argument(
        "--py6s-ozone",
        type=float,
        default=0.3,
        help="Py6S ozone (cm-atm).",
    )
    parser.add_argument(
        "--py6s-output-scale-factor",
        type=float,
        default=10000.0,
        help="Optional scale factor applied to Py6S reflectance output.",
    )
    parser.add_argument(
        "--py6s-output-dtype",
        default="int16",
        help="Output dtype written by Py6S atmospheric correction.",
    )
    parser.add_argument(
        "--py6s-executable",
        help="Optional full path to the 6S executable (e.g. sixs or sixsV1.1).",
    )
    parser.add_argument(
        "--py6s-use-imd-radiance-calibration",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Use WorldView IMD absCalFactor/effectiveBandwidth for per-band DN->radiance "
            "before Py6S. Enabled by default; disable with --no-py6s-use-imd-radiance-calibration."
        ),
    )
    parser.add_argument(
        "--py6s-use-worldview-gain-offset-adjustment",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Apply built-in WV02/WV03 gain/offset adjustment table on top of IMD absCalFactor/effectiveBandwidth. "
            "Enabled by default."
        ),
    )
    parser.add_argument(
        "--py6s-auto-atmos-source",
        choices=["none", "nasa_power"],
        default="none",
        help=(
            "Optional auto-source for Py6S aot550/water_vapor/ozone. "
            "Applied when --py6s-atmosphere-profile=user."
        ),
    )
    parser.add_argument(
        "--py6s-auto-atmos-grid-size",
        type=int,
        default=3,
        help="NxN sample grid size over scene bbox for auto atmosphere fetch (default: 3).",
    )
    parser.add_argument(
        "--py6s-auto-atmos-search-days",
        type=int,
        default=1,
        help="Search +/- N days around scene date for auto atmosphere fetch (default: 1).",
    )
    parser.add_argument(
        "--py6s-auto-atmos-timeout-s",
        type=float,
        default=30.0,
        help="HTTP timeout in seconds for auto atmosphere API calls (default: 30).",
    )
    parser.add_argument(
        "--py6s-auto-atmos-power-endpoint",
        default="https://power.larc.nasa.gov/api/temporal/daily/point",
        help="NASA POWER endpoint for auto atmosphere fetch.",
    )
    parser.add_argument(
        "--footprint-epsg",
        type=int,
        default=4326,
        help="EPSG assigned to Maxar shapefile footprints before conversion to GeoPackage.",
    )
    parser.add_argument(
        "--filter-basename",
        action="append",
        help=(
            "Optional image basenames to process. Repeat or pass comma-separated values. "
            "Accepts mul or pan photo basename."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default="Scene_Output",
        help="Directory for final full-scene outputs.",
    )
    parser.add_argument(
        "--output-suffix",
        default="_final",
        help="Suffix appended to final scene output filenames before .tif extension.",
    )
    parser.add_argument(
        "--scratch-dir",
        default="/tmp",
        help="Scratch directory for temporary intermediate files.",
    )
    parser.add_argument(
        "--run-atmospheric-correction",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable/disable atmospheric correction stage (default: enabled).",
    )
    parser.add_argument(
        "--run-pansharpen",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable/disable pansharpen stage (default: enabled).",
    )
    parser.add_argument(
        "--skip-flaash",
        action="store_true",
        help="Skip FLAASH and use --existing-flaash-input as multispectral input (flaash method only).",
    )
    parser.add_argument(
        "--existing-flaash-input",
        help="Path to existing multispectral .dat file used when --skip-flaash is set.",
    )
    parser.add_argument(
        "--existing-mul-ortho-input",
        help="Optional path to existing orthorectified multispectral raster. If set with --existing-pan-ortho-input, workflow resumes at pansharpen.",
    )
    parser.add_argument(
        "--existing-pan-ortho-input",
        help="Optional path to existing orthorectified panchromatic raster. Requires --existing-mul-ortho-input.",
    )
    parser.add_argument(
        "--cloud-mask-command",
        help=(
            "Optional shell command template to run on the full-scene image. "
            "Template variables: {input}, {output}, {scene_root}, {image_basename}."
        ),
    )
    parser.add_argument(
        "--run-cloud-mask",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Enable/disable cloud masking stage. Enabled by default; use "
            "--no-run-cloud-mask to skip without changing config."
        ),
    )
    parser.add_argument(
        "--cloud-mask-method",
        choices=["omnicloudmask"],
        default="omnicloudmask",
        help=(
            "Optional built-in cloud masking method. "
            "Mask is inferred from orthorectified multispectral and applied to current workflow output."
        ),
    )
    parser.add_argument(
        "--cloud-mask-red-band-index",
        type=int,
        default=5,
        help="1-based red band index for built-in omnicloudmask mode (WorldView default: 5).",
    )
    parser.add_argument(
        "--cloud-mask-green-band-index",
        type=int,
        default=3,
        help="1-based green band index for built-in omnicloudmask mode (WorldView default: 3).",
    )
    parser.add_argument(
        "--cloud-mask-nir-band-index",
        type=int,
        default=7,
        help="1-based NIR band index for built-in omnicloudmask mode (WorldView default: 7).",
    )
    parser.add_argument(
        "--cloud-mask-classes",
        default="1,2,3",
        help="Comma-separated class ids treated as cloudy/cloud-shadow in omnicloudmask output.",
    )
    parser.add_argument(
        "--cloud-buffer-pixels",
        type=int,
        default=10,
        help="Optional binary dilation size in pixels applied around cloud mask regions.",
    )
    parser.add_argument(
        "--cloud-mask-inference-resolution-m",
        type=float,
        default=10.0,
        help=(
            "Target map resolution for OmniCloudMask inference input (default: 10.0). "
            "Mask is reprojected back to final output grid when applied."
        ),
    )
    parser.add_argument(
        "--cloud-mask-omnicloud-kwargs-json",
        help='Optional JSON object passed directly to omnicloudmask predict_from_array, e.g. \'{"batch_size": 8}\'.',
    )
    parser.add_argument(
        "--cloud-mask-output-suffix",
        default="_cloudmasked",
        help="Suffix added before extension for cloud-masked image output.",
    )
    parser.add_argument(
        "--cloud-mask-mask-suffix",
        default="_cloudmask",
        help="Suffix added before extension for binary cloud mask output.",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """Parse CLI/config arguments, validate them, and execute the workflow."""
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument("--config-yaml")
    config_args, _ = config_parser.parse_known_args(argv)

    config_defaults: Dict = {}
    if config_args.config_yaml:
        try:
            config_defaults = _normalize_config_defaults(_load_yaml_config(config_args.config_yaml))
        except Exception as exc:
            raise SystemExit(f"Failed to load --config-yaml: {exc}")

    parser = build_parser()
    if config_defaults:
        parser.set_defaults(**config_defaults)
    args = parser.parse_args(argv)

    if not args.input_dir:
        parser.error("--input-dir is required (via CLI or --config-yaml).")
    if not args.dem_file_path and not args.existing_mul_ortho_input:
        parser.error(
            "--dem-file-path is required unless --existing-mul-ortho-input is provided."
        )
    if (
        args.run_atmospheric_correction
        and args.atmospheric_method == "flaash"
        and not args.skip_flaash
        and not args.envi_engine_path
        and not (args.existing_mul_ortho_input and (args.existing_pan_ortho_input or not args.run_pansharpen))
    ):
        parser.error("--envi-engine-path is required when running FLAASH atmospheric correction.")
    if args.flaash_dem_ground_percentile < 0 or args.flaash_dem_ground_percentile > 100:
        parser.error("--flaash-dem-ground-percentile must be between 0 and 100.")
    if args.run_atmospheric_correction and args.atmospheric_method == "flaash" and not args.skip_flaash:
        if not re.match(r"^/mnt/[a-zA-Z]/", args.scratch_dir):
            parser.error(
                "--scratch-dir must be under /mnt/<drive>/... when running FLAASH with Windows ENVI."
            )
    if (
        args.run_atmospheric_correction
        and args.atmospheric_method == "flaash"
        and args.skip_flaash
        and not args.existing_flaash_input
    ):
        if not (args.existing_mul_ortho_input and (args.existing_pan_ortho_input or not args.run_pansharpen)):
            parser.error("--existing-flaash-input is required when --skip-flaash is set.")
    if args.run_pansharpen and args.existing_mul_ortho_input and not args.existing_pan_ortho_input:
        parser.error("--existing-pan-ortho-input is required when --run-pansharpen is enabled with --existing-mul-ortho-input.")
    for path_arg, path_value in (
        ("--existing-mul-ortho-input", args.existing_mul_ortho_input),
        ("--existing-pan-ortho-input", args.existing_pan_ortho_input if args.run_pansharpen else None),
    ):
        if path_value and not os.path.isfile(path_value):
            parser.error(f"{path_arg} does not exist: {path_value}")
    if not os.path.isdir(args.scratch_dir):
        parser.error(f"--scratch-dir does not exist or is not a directory: {args.scratch_dir}")
    if args.run_cloud_mask and args.cloud_mask_method and args.cloud_mask_command:
        parser.error("Use either --cloud-mask-method or --cloud-mask-command, not both.")
    if args.cloud_mask_inference_resolution_m <= 0:
        parser.error("--cloud-mask-inference-resolution-m must be > 0.")
    if args.py6s_auto_atmos_grid_size < 1:
        parser.error("--py6s-auto-atmos-grid-size must be >= 1.")
    if args.py6s_auto_atmos_search_days < 0:
        parser.error("--py6s-auto-atmos-search-days must be >= 0.")
    if args.py6s_auto_atmos_timeout_s <= 0:
        parser.error("--py6s-auto-atmos-timeout-s must be > 0.")
    try:
        _parse_json_dict(args.cloud_mask_omnicloud_kwargs_json)
    except Exception as exc:
        parser.error(f"Invalid --cloud-mask-omnicloud-kwargs-json: {exc}")

    return run_workflow(args)


if __name__ == "__main__":
    sys.exit(main())
