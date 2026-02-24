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
from typing import Dict, Iterable, List, Optional

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)


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
    from vhrharmonize import get_image_percentile_value

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


def _parse_json_dict(raw_json: Optional[str]) -> Dict:
    if not raw_json:
        return {}
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


def _sanitize_for_filename(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in str(value))


def _load_tiles(tile_gpkg_path: str, tile_layer: Optional[str], tile_id_field: str):
    import geopandas as gpd

    tile_gdf = gpd.read_file(tile_gpkg_path, layer=tile_layer)
    if tile_gdf.empty:
        raise ValueError(f"Tile layer is empty: {tile_gpkg_path}")
    if tile_id_field not in tile_gdf.columns:
        raise ValueError(f"Tile id field '{tile_id_field}' not found in {tile_gpkg_path}")
    if tile_gdf.crs is None:
        raise ValueError("Tile GPKG has no CRS. Please define CRS on tile geometries.")
    return tile_gdf


def _buffer_tiles(tile_gdf, tile_buffer_m: float):
    if tile_buffer_m == 0:
        return tile_gdf
    buffered = tile_gdf.copy()
    buffered["geometry"] = buffered.geometry.buffer(tile_buffer_m)
    return buffered


def _get_scene_footprint_geometry(mul_shp_path: str, footprint_epsg: int, target_epsg: int):
    import geopandas as gpd

    scene_gdf = gpd.read_file(mul_shp_path)
    if scene_gdf.empty:
        return None
    if scene_gdf.crs is None:
        print(
            f"CRS update: scene footprint has no CRS, assigning EPSG:{footprint_epsg} "
            f"for {os.path.basename(mul_shp_path)}"
        )
        scene_gdf = scene_gdf.set_crs(epsg=footprint_epsg)
    scene_epsg = scene_gdf.crs.to_epsg()
    if scene_epsg != target_epsg:
        print(
            f"CRS update: reprojecting scene footprint from {scene_gdf.crs} to EPSG:{target_epsg} "
            f"for {os.path.basename(mul_shp_path)}"
        )
        scene_gdf = scene_gdf.to_crs(epsg=target_epsg)
    # GeoPandas deprecates unary_union in favor of union_all.
    if hasattr(scene_gdf.geometry, "union_all"):
        return scene_gdf.geometry.union_all()
    return scene_gdf.geometry.unary_union


def _find_overlapping_tiles(scene_geom, tile_gdf):
    if scene_geom is None:
        return tile_gdf.iloc[0:0]

    sindex = tile_gdf.sindex
    if sindex is not None:
        candidate_idx = list(sindex.intersection(scene_geom.bounds))
        if not candidate_idx:
            return tile_gdf.iloc[0:0]
        candidates = tile_gdf.iloc[candidate_idx]
    else:
        candidates = tile_gdf
    return candidates[candidates.geometry.intersects(scene_geom)]


def _write_clipped_tile(
    input_raster_path: str,
    tile_geometry,
    output_tile_path: str,
):
    import rasterio
    import rasterio.mask

    os.makedirs(os.path.dirname(output_tile_path), exist_ok=True)
    with rasterio.open(input_raster_path) as src:
        clipped_data, clipped_transform = rasterio.mask.mask(
            src,
            [tile_geometry.__geo_interface__],
            crop=True,
            nodata=src.nodata,
        )
        profile = src.profile.copy()
        profile.update(
            height=clipped_data.shape[1],
            width=clipped_data.shape[2],
            transform=clipped_transform,
        )
        with rasterio.open(output_tile_path, "w", **profile) as dst:
            dst.write(clipped_data)


def _reproject_geometry(geometry, source_crs: str, target_crs: str):
    import geopandas as gpd

    geom_gs = gpd.GeoSeries([geometry], crs=source_crs)
    return geom_gs.to_crs(target_crs).iloc[0]


def _resolve_dimension_name(points_dtype_names, requested_dimension: str) -> str:
    if requested_dimension in points_dtype_names:
        return requested_dimension
    lower_map = {name.lower(): name for name in points_dtype_names}
    resolved = lower_map.get(requested_dimension.lower())
    if resolved:
        return resolved
    raise ValueError(
        f"Requested LiDAR dimension '{requested_dimension}' not found. "
        f"Available dimensions: {sorted(points_dtype_names)}"
    )


def _write_lidar_mean_raster_from_points(
    points,
    dimension_name: str,
    bounds_xy,
    resolution: float,
    output_raster_path: str,
    output_crs: str,
    nodata_value: float,
) -> None:
    import numpy as np
    import rasterio
    from rasterio.transform import from_origin

    min_x, min_y, max_x, max_y = bounds_xy
    if max_x <= min_x or max_y <= min_y:
        raise ValueError(f"Invalid bounds for LiDAR rasterization: {bounds_xy}")

    width = max(1, int(np.ceil((max_x - min_x) / resolution)))
    height = max(1, int(np.ceil((max_y - min_y) / resolution)))

    mean_grid = np.full((height, width), nodata_value, dtype=np.float32)

    if points is not None and points.size > 0:
        x = points["X"].astype(np.float64, copy=False)
        y = points["Y"].astype(np.float64, copy=False)
        values = points[dimension_name].astype(np.float64, copy=False)

        valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(values)
        x = x[valid]
        y = y[valid]
        values = values[valid]

        if x.size > 0:
            col = np.floor((x - min_x) / resolution).astype(np.int64)
            row = np.floor((max_y - y) / resolution).astype(np.int64)
            in_bounds = (col >= 0) & (col < width) & (row >= 0) & (row < height)
            row = row[in_bounds]
            col = col[in_bounds]
            values = values[in_bounds]

            if values.size > 0:
                sums = np.zeros((height, width), dtype=np.float64)
                counts = np.zeros((height, width), dtype=np.uint32)
                np.add.at(sums, (row, col), values)
                np.add.at(counts, (row, col), 1)
                nonzero = counts > 0
                mean_grid[nonzero] = (sums[nonzero] / counts[nonzero]).astype(np.float32)

    transform = from_origin(min_x, max_y, resolution, resolution)
    os.makedirs(os.path.dirname(output_raster_path), exist_ok=True)
    with rasterio.open(
        output_raster_path,
        "w",
        driver="GTiff",
        width=width,
        height=height,
        count=1,
        dtype="float32",
        crs=output_crs,
        transform=transform,
        nodata=nodata_value,
    ) as dst:
        dst.write(mean_grid, 1)


def _run_lidar_tile_means(args: argparse.Namespace, tile_gdf_buffered, tile_crs) -> None:
    from pyforestscan.handlers import read_lidar
    from pyforestscan.utils import get_srs_from_ept

    lidar_srs = args.lidar_srs or get_srs_from_ept(args.lidar_ept_path)
    if not lidar_srs:
        raise ValueError(
            "Unable to determine LiDAR SRS from EPT metadata. Set --lidar-srs explicitly."
        )

    lidar_out_dir = args.lidar_output_dir or os.path.join(args.tile_output_dir, "lidar_mean")
    os.makedirs(lidar_out_dir, exist_ok=True)
    print(f"Running LiDAR mean raster step from EPT: {args.lidar_ept_path}")
    print(f"LiDAR SRS: {lidar_srs}; requested dimension: {args.lidar_dimension}")

    for _, tile_row in tile_gdf_buffered.iterrows():
        tile_id = _sanitize_for_filename(tile_row[args.tile_id_field])
        tile_geom = tile_row.geometry
        tile_geom_lidar = _reproject_geometry(tile_geom, tile_crs, lidar_srs)
        min_x, min_y, max_x, max_y = tile_geom_lidar.bounds
        tile_bounds = ([min_x, max_x], [min_y, max_y])
        print(f"LiDAR tile {tile_id}: reading EPT points within bounds {tile_bounds}")

        arrays = read_lidar(
            input_file=args.lidar_ept_path,
            srs=lidar_srs,
            bounds=tile_bounds,
            thin_radius=None,
            hag=False,
            hag_dtm=False,
            dtm=None,
            crop_poly=False,
            poly=None,
        )
        points = arrays[0] if arrays and len(arrays) > 0 else None

        dimension_name = None
        if points is not None and points.size > 0:
            dimension_name = _resolve_dimension_name(points.dtype.names, args.lidar_dimension)
        else:
            # Keep requested name for output naming if tile has no points.
            dimension_name = args.lidar_dimension
            print(f"LiDAR tile {tile_id}: no points found in bounds; writing nodata raster.")

        output_name = f"{tile_id}_{dimension_name}{args.lidar_output_suffix}.tif"
        output_path = os.path.join(lidar_out_dir, output_name)
        if os.path.exists(output_path) and not args.lidar_overwrite:
            print(f"LiDAR tile {tile_id}: output exists, skipping ({output_path})")
            continue
        _write_lidar_mean_raster_from_points(
            points=points,
            dimension_name=dimension_name,
            bounds_xy=(min_x, min_y, max_x, max_y),
            resolution=args.lidar_resolution,
            output_raster_path=output_path,
            output_crs=lidar_srs,
            nodata_value=args.lidar_nodata,
        )
        print(f"LiDAR tile {tile_id}: wrote {output_path}")


def run_workflow(args: argparse.Namespace) -> int:
    from vhrharmonize import (
        apply_binary_cloud_mask_to_image,
        create_cloud_mask_with_omnicloudmask,
        find_files,
        find_roots,
        gcp_refined_rpc_orthorectification,
        get_metadata_from_files,
        pansharpen_image,
        run_flaash,
        shp_to_gpkg,
    )

    filter_basenames = _parse_filter_basenames(args.filter_basename)
    use_existing_ortho_inputs = bool(
        args.existing_mul_ortho_input and args.existing_pan_ortho_input
    )
    cloud_classes = _parse_int_csv(args.cloud_mask_classes)
    omnicloud_kwargs = _parse_json_dict(args.cloud_mask_omnicloud_kwargs_json)
    tile_gdf = _load_tiles(args.tile_gpkg, args.tile_layer, args.tile_id_field)
    print(f"Tile CRS: {tile_gdf.crs}; requested output CRS: EPSG:{args.epsg}")
    tile_epsg = tile_gdf.crs.to_epsg()
    if tile_epsg != args.epsg:
        print(f"CRS update: reprojecting tile polygons from {tile_gdf.crs} to EPSG:{args.epsg}")
        tile_gdf = tile_gdf.to_crs(epsg=args.epsg)
    else:
        print("CRS check: tile CRS already matches output CRS.")
    tile_crs = tile_gdf.crs
    tile_gdf_buffered = tile_gdf
    if args.tile_buffer_m != 0:
        print(f"Applying tile buffer: {args.tile_buffer_m} map units")
        tile_gdf_buffered = _buffer_tiles(tile_gdf, args.tile_buffer_m)
    os.makedirs(args.tile_output_dir, exist_ok=True)

    if args.lidar_ept_path:
        _run_lidar_tile_means(args, tile_gdf_buffered, str(tile_crs))

    envi_engine = None
    if not args.skip_flaash and not use_existing_ortho_inputs:
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

                required = (
                    "mul_imd_file",
                    "mul_tif_file",
                    "pan_imd_file",
                    "pan_tif_file",
                    "mul_shp_file",
                    "pan_shp_file",
                )
                if not _required_keys_present(found_default_file, required):
                    print("Skipping scene: required files missing for workflow.")
                    continue

                mul_imd_file = found_default_file["mul_imd_file"]
                pan_imd_file = found_default_file["pan_imd_file"]
                mul_tif_file = found_default_file["mul_tif_file"]
                pan_tif_file = found_default_file["pan_tif_file"]
                mul_shp_path = found_default_file["mul_shp_file"]
                pan_shp_path = found_default_file["pan_shp_file"]

                scene_geom = _get_scene_footprint_geometry(
                    mul_shp_path,
                    footprint_epsg=args.footprint_epsg,
                    target_epsg=args.epsg,
                )

                params_overrides_scene, mul_params_overrides_photo, mul_imd_data = (
                    get_metadata_from_files(root_file_path, mul_imd_file, mul_photo_basename)
                )
                _, _, pan_imd_data = get_metadata_from_files(
                    root_file_path, pan_imd_file, pan_photo_basename
                )

                overlapping_tiles = _find_overlapping_tiles(scene_geom, tile_gdf_buffered)
                if overlapping_tiles.empty:
                    print("Skipping scene: no overlap with buffered tile polygons.")
                    continue

                with tempfile.TemporaryDirectory(dir=args.scratch_dir, prefix="vhr_scene_") as work_dir:
                    ortho_folder_name = "OrthoFromDefaultRPC"
                    if use_existing_ortho_inputs:
                        print("Using existing orthorectified inputs; skipping FLAASH and orthorectification")
                        mul_flaash_ortho_path = args.existing_mul_ortho_input
                        pan_ortho_path = args.existing_pan_ortho_input
                    else:
                        print("Creating footprint GeoPackage for FLAASH terrain statistics")
                        mul_gpkg_path = os.path.join(work_dir, f"{mul_photo_basename}.gpkg")
                        shp_to_gpkg(mul_shp_path, mul_gpkg_path, args.footprint_epsg)

                        mul_flaash_image_path = os.path.join(work_dir, f"{mul_photo_basename}_FLAASH.dat")
                        mul_flaash_params_path = os.path.join(work_dir, f"{mul_photo_basename}_FLAASH_Params.txt")

                        if args.skip_flaash:
                            mul_flaash_input_path = args.existing_flaash_input
                        else:
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

                        print("Using default RPC model for orthorectification (no GCP refinement)")

                        print("Orthorectifying multispectral image")
                        mul_flaash_ortho_path = os.path.join(
                            work_dir,
                            f"{mul_photo_basename}_FLAASH_{ortho_folder_name}.tif",
                        )
                        gcp_refined_rpc_orthorectification(
                            mul_flaash_input_path,
                            mul_flaash_ortho_path,
                            args.dem_file_path,
                            args.epsg,
                            output_nodata_value=args.nodata_value,
                            dtype=args.dtype,
                            output_resolution=mul_imd_data.get("product_res"),
                        )

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
                            output_resolution=pan_imd_data.get("product_res"),
                        )

                    print("Pansharpening multispectral image")
                    mul_pansharp_path = os.path.join(
                        work_dir,
                        f"{mul_photo_basename}_FLAASH_{ortho_folder_name}_Pansharps.tif",
                    )
                    pansharpen_image(
                        mul_flaash_ortho_path,
                        pan_ortho_path,
                        mul_pansharp_path,
                        change_nodata_value=args.nodata_value,
                    )

                    print(f"Writing final clipped outputs for {len(overlapping_tiles)} overlapping tile(s)")
                    for _, tile_row in overlapping_tiles.iterrows():
                        tile_id = _sanitize_for_filename(tile_row[args.tile_id_field])
                        tile_filename = (
                            f"{tile_id}_{mul_photo_basename}{args.tile_output_suffix}.tif"
                        )
                        tile_output_path = os.path.join(args.tile_output_dir, tile_filename)
                        tile_work_input_path = os.path.join(
                            work_dir,
                            f"{tile_id}_{mul_photo_basename}_tile_input.tif",
                        )
                        try:
                            _write_clipped_tile(
                                mul_pansharp_path,
                                tile_row.geometry,
                                tile_work_input_path,
                            )
                        except ValueError:
                            print(
                                f"Skipping tile {tile_id}: no overlap with pansharpened raster extent."
                            )
                            continue

                        final_tile_path = tile_work_input_path

                        if args.cloud_mask_command:
                            cloud_mask_output_path = _build_cloud_mask_output_path(
                                tile_work_input_path,
                                args.cloud_mask_output_suffix,
                            )
                            _run_cloud_mask_command(
                                args.cloud_mask_command,
                                tile_work_input_path,
                                cloud_mask_output_path,
                                root_folder_path,
                                mul_photo_basename,
                            )
                            final_tile_path = cloud_mask_output_path

                        if args.cloud_mask_method == "omnicloudmask":
                            mask_output_path = _build_cloud_mask_output_path(
                                tile_work_input_path,
                                args.cloud_mask_mask_suffix,
                            )
                            masked_image_output_path = _build_cloud_mask_output_path(
                                tile_work_input_path,
                                args.cloud_mask_output_suffix,
                            )
                            print(f"Running built-in OmniCloudMask cloud masking for tile {tile_id}")
                            create_cloud_mask_with_omnicloudmask(
                                tile_work_input_path,
                                mask_output_path,
                                red_band_index=args.cloud_mask_red_band_index,
                                green_band_index=args.cloud_mask_green_band_index,
                                nir_band_index=args.cloud_mask_nir_band_index,
                                cloud_classes=cloud_classes,
                                buffer_pixels=args.cloud_buffer_pixels,
                                omnicloud_kwargs=omnicloud_kwargs,
                            )
                            apply_binary_cloud_mask_to_image(
                                tile_work_input_path,
                                mask_output_path,
                                masked_image_output_path,
                                output_nodata_value=args.nodata_value,
                            )
                            final_tile_path = masked_image_output_path

                        shutil.copy2(final_tile_path, tile_output_path)
                print("Scene complete")

    print("All processing complete")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run WorldView preprocessing (FLAASH + orthorectify + pansharpen) from CLI.",
    )
    parser.add_argument(
        "--config-yaml",
        help="Optional YAML config file. Keys should match argument names (use '_' or '-').",
    )
    parser.add_argument(
        "--input-dir",
        action="append",
        help="Input directory to scan for scene folders and Root*.txt files. Repeat for multiple directories.",
    )
    parser.add_argument("--dem-file-path", help="DEM GeoTIFF path in WGS84 ellipsoidal height.")
    parser.add_argument("--envi-engine-path", help="Path to ENVI taskengine executable.")
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
    parser.add_argument("--tile-gpkg", help="Tile GeoPackage path for tile-driven clipping mode.")
    parser.add_argument("--tile-layer", help="Optional tile layer name inside --tile-gpkg.")
    parser.add_argument(
        "--tile-id-field",
        default="tile_id",
        help="Tile attribute field used to name output tiles.",
    )
    parser.add_argument(
        "--tile-output-dir",
        default="Tile_Output",
        help="Directory for final clipped tile outputs when --tile-gpkg is used.",
    )
    parser.add_argument(
        "--tile-output-suffix",
        default="_final",
        help="Suffix appended to final output filenames before .tif extension.",
    )
    parser.add_argument(
        "--tile-buffer-m",
        type=float,
        default=100.0,
        help=(
            "Buffer distance applied to tile geometries before overlap checks and clipping "
            "(in tile CRS map units; for EPSG:6635 this is meters)."
        ),
    )
    parser.add_argument(
        "--lidar-ept-path",
        help="Optional EPT metadata path (e.g., /path/to/ept.json) for LiDAR mean-raster tile outputs.",
    )
    parser.add_argument(
        "--lidar-srs",
        help="Optional LiDAR SRS override for EPT read bounds (e.g., EPSG:6635).",
    )
    parser.add_argument(
        "--lidar-dimension",
        default="Intensity",
        help="LiDAR dimension to aggregate as mean (e.g., Intensity, Z, HeightAboveGround).",
    )
    parser.add_argument(
        "--lidar-resolution",
        type=float,
        help="Output resolution for LiDAR mean rasters in LiDAR CRS map units.",
    )
    parser.add_argument(
        "--lidar-output-dir",
        help="Output directory for LiDAR mean tile rasters (default: <tile-output-dir>/lidar_mean).",
    )
    parser.add_argument(
        "--lidar-output-suffix",
        default="_mean",
        help="Suffix appended to LiDAR mean tile filenames before extension.",
    )
    parser.add_argument(
        "--lidar-nodata",
        type=float,
        default=-9999.0,
        help="NoData value for LiDAR mean raster outputs.",
    )
    parser.add_argument(
        "--lidar-overwrite",
        action="store_true",
        help="Overwrite existing LiDAR mean tile rasters. By default existing outputs are skipped.",
    )
    parser.add_argument(
        "--scratch-dir",
        default="/tmp",
        help="Scratch directory for temporary intermediate files.",
    )
    parser.add_argument(
        "--skip-flaash",
        action="store_true",
        help="Skip FLAASH and use --existing-flaash-input as the multispectral input for orthorectification.",
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
            "Optional shell command template to run per tile after clipping. "
            "Template variables: {input}, {output}, {scene_root}, {image_basename}."
        ),
    )
    parser.add_argument(
        "--cloud-mask-method",
        choices=["omnicloudmask"],
        help="Optional built-in cloud masking method to run after pansharpening.",
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
        default=0,
        help="Optional binary dilation size in pixels applied around cloud mask regions.",
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
    if not args.dem_file_path:
        parser.error("--dem-file-path is required (via CLI or --config-yaml).")
    if not args.tile_gpkg:
        parser.error("--tile-gpkg is required (tile-driven mode is mandatory).")
    if not args.skip_flaash and not args.envi_engine_path:
        if not (args.existing_mul_ortho_input and args.existing_pan_ortho_input):
            parser.error("--envi-engine-path is required unless --skip-flaash is set.")
    if args.tile_buffer_m < 0:
        parser.error("--tile-buffer-m must be >= 0.")
    if args.lidar_ept_path and (args.lidar_resolution is None or args.lidar_resolution <= 0):
        parser.error("--lidar-resolution must be > 0 when --lidar-ept-path is set.")
    if args.flaash_dem_ground_percentile < 0 or args.flaash_dem_ground_percentile > 100:
        parser.error("--flaash-dem-ground-percentile must be between 0 and 100.")
    if not args.skip_flaash:
        if not re.match(r"^/mnt/[a-zA-Z]/", args.scratch_dir):
            parser.error(
                "--scratch-dir must be under /mnt/<drive>/... when running FLAASH with Windows ENVI."
            )
    if args.skip_flaash and not args.existing_flaash_input:
        if not (args.existing_mul_ortho_input and args.existing_pan_ortho_input):
            parser.error("--existing-flaash-input is required when --skip-flaash is set.")
    if bool(args.existing_mul_ortho_input) ^ bool(args.existing_pan_ortho_input):
        parser.error(
            "--existing-mul-ortho-input and --existing-pan-ortho-input must be provided together."
        )
    for path_arg, path_value in (
        ("--existing-mul-ortho-input", args.existing_mul_ortho_input),
        ("--existing-pan-ortho-input", args.existing_pan_ortho_input),
    ):
        if path_value and not os.path.isfile(path_value):
            parser.error(f"{path_arg} does not exist: {path_value}")
    if not os.path.isdir(args.scratch_dir):
        parser.error(f"--scratch-dir does not exist or is not a directory: {args.scratch_dir}")
    if args.cloud_mask_method and args.cloud_mask_command:
        parser.error("Use either --cloud-mask-method or --cloud-mask-command, not both.")
    try:
        _parse_json_dict(args.cloud_mask_omnicloud_kwargs_json)
    except Exception as exc:
        parser.error(f"Invalid --cloud-mask-omnicloud-kwargs-json: {exc}")

    return run_workflow(args)


if __name__ == "__main__":
    sys.exit(main())
