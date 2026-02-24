#!/usr/bin/env python3
"""Tile a source image, generate tile-wise LiDAR mean rasters from EPT, and run elastix alignment."""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional

import geopandas as gpd
import numpy as np
import rasterio
import rasterio.mask
from rasterio.transform import from_origin


def _load_yaml_config(config_yaml_path: str) -> Dict:
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError("PyYAML is required for --config-yaml.") from exc

    with open(config_yaml_path, "r", encoding="utf-8") as f:
        loaded = yaml.safe_load(f) or {}
    if not isinstance(loaded, dict):
        raise ValueError("Config YAML root must be a mapping/dictionary.")
    return {k.replace("-", "_"): v for k, v in loaded.items()}


def _sanitize_for_filename(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in str(value))


def _load_tiles(tile_gpkg_path: str, tile_layer: Optional[str], tile_id_field: str):
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


def _write_clipped_tile(input_raster_path: str, tile_geometry, output_tile_path: str):
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


def _count_valid_pixels(raster_path: str, band_index: int = 1) -> int:
    with rasterio.open(raster_path) as src:
        data = src.read(band_index)
        nodata = src.nodata
    valid = np.isfinite(data)
    if nodata is not None:
        valid &= data != nodata
    return int(valid.sum())


def _reproject_geometry(geometry, source_crs: str, target_crs: str):
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
    min_x, min_y, max_x, max_y = bounds_xy
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
        if values.size > 0:
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


def _extract_single_band(input_path: str, output_path: str, band_index: int) -> None:
    cmd = ["gdal_translate", "-q", "-b", str(band_index), input_path, output_path]
    subprocess.run(cmd, check=True)


def _wsl_to_windows_path(path: str) -> str:
    if path.startswith("/mnt/") and len(path) > 6 and path[5].isalpha() and path[6] == "/":
        drive = path[5].upper()
        rest = path[7:].replace("/", "\\")
        return f"{drive}:\\{rest}"
    return path


def _get_fixed_grid_params(fixed_path: str):
    with rasterio.open(fixed_path) as src:
        bounds = src.bounds
        transform = src.transform
        crs = src.crs
    if crs is None:
        raise ValueError(f"Fixed raster has no CRS: {fixed_path}")
    epsg = crs.to_epsg()
    if epsg is None:
        raise ValueError(f"Fixed raster CRS has no EPSG code: {fixed_path} ({crs})")
    return epsg, bounds.left, bounds.top, bounds.right, bounds.bottom, abs(transform.a)


def _run_elastix_pair(
    *,
    elastix_script: str,
    fixed_path: str,
    moving_path: str,
    outdir: str,
    moving_band_index: int,
    resolution_override: Optional[float] = None,
) -> None:
    epsg, ulx, uly, lrx, lry, res = _get_fixed_grid_params(fixed_path)
    if resolution_override is not None:
        res = float(resolution_override)
    with tempfile.TemporaryDirectory(prefix="elx_tile_") as temp_dir:
        moving_band_path = os.path.join(temp_dir, "moving_band.tif")
        _extract_single_band(moving_path, moving_band_path, moving_band_index)
        cmd = [
            elastix_script,
            "--fixed",
            fixed_path,
            "--moving",
            moving_band_path,
            "--outdir",
            outdir,
            "--epsg",
            str(epsg),
            "--ullr",
            str(ulx),
            str(uly),
            str(lrx),
            str(lry),
            "--res",
            str(res),
        ]
        subprocess.run(cmd, check=True)


def _apply_transform_to_all_bands(
    *,
    moving_image_path: str,
    work_outdir: str,
    epsg: int,
    ulx: float,
    uly: float,
    lrx: float,
    lry: float,
    res: float,
) -> str:
    transformix_exe = os.environ["TRANSFORMIX_EXE"]
    tp_rigid = os.path.join(work_outdir, "elastix_rigid", "TransformParameters.0.txt")
    if not os.path.exists(tp_rigid):
        raise FileNotFoundError(f"Missing rigid transform file: {tp_rigid}")

    with rasterio.open(moving_image_path) as src:
        band_count = src.count

    band_results = []
    for band_idx in range(1, band_count + 1):
        band_dir = os.path.join(work_outdir, "warp_all_bands", f"band_{band_idx:02d}")
        tmp_dir = os.path.join(band_dir, "tmp")
        os.makedirs(tmp_dir, exist_ok=True)
        moving_band = os.path.join(tmp_dir, f"moving_band_{band_idx:02d}.tif")
        moving_band_u16 = os.path.join(tmp_dir, f"moving_band_{band_idx:02d}_u16.tif")
        moving_roi = os.path.join(tmp_dir, f"moving_band_{band_idx:02d}_roi.tif")
        result_georef = os.path.join(band_dir, "result_georef.tif")

        subprocess.run(
            ["gdal_translate", "-q", "-b", str(band_idx), moving_image_path, moving_band],
            check=True,
        )
        subprocess.run(
            ["gdal_translate", "-q", "-of", "GTiff", "-ot", "UInt16", "-co", "TILED=YES", "-co", "COMPRESS=NONE", moving_band, moving_band_u16],
            check=True,
        )
        subprocess.run(
            [
                "gdalwarp",
                "-q",
                "-te",
                str(ulx),
                str(lry),
                str(lrx),
                str(uly),
                "-tap",
                "-tr",
                str(res),
                str(res),
                "-r",
                "near",
                moving_band_u16,
                moving_roi,
            ],
            check=True,
        )

        subprocess.run(
            [
                transformix_exe,
                "-in",
                _wsl_to_windows_path(moving_roi),
                "-tp",
                _wsl_to_windows_path(tp_rigid),
                "-out",
                _wsl_to_windows_path(band_dir),
            ],
            check=True,
        )

        raw_result = os.path.join(band_dir, "result.tif")
        subprocess.run(
            [
                "gdal_translate",
                "-q",
                "-a_srs",
                f"EPSG:{epsg}",
                "-a_ullr",
                str(ulx),
                str(uly),
                str(lrx),
                str(lry),
                raw_result,
                result_georef,
            ],
            check=True,
        )
        band_results.append(result_georef)

    if not band_results:
        raise RuntimeError("No transformed band outputs generated.")

    if len(band_results) == 1:
        return band_results[0]

    vrt_path = os.path.join(work_outdir, "warp_all_bands", "result_aligned.vrt")
    final_tif = os.path.join(work_outdir, "warp_all_bands", "result_aligned.tif")
    subprocess.run(["gdalbuildvrt", "-q", "-separate", vrt_path, *band_results], check=True)
    subprocess.run(["gdal_translate", "-q", vrt_path, final_tif], check=True)
    return final_tif


def run_workflow(args: argparse.Namespace) -> int:
    from pyforestscan.handlers import read_lidar
    from pyforestscan.utils import get_srs_from_ept

    if "ELASTIX_EXE" not in os.environ or "TRANSFORMIX_EXE" not in os.environ:
        raise EnvironmentError("ELASTIX_EXE and TRANSFORMIX_EXE must be set in environment.")

    tile_gdf = _load_tiles(args.tile_gpkg, args.tile_layer, args.tile_id_field)
    image_name = os.path.splitext(os.path.basename(args.image_path))[0]

    with rasterio.open(args.image_path) as src:
        image_crs = src.crs
    if image_crs is None:
        raise ValueError(f"Input image has no CRS: {args.image_path}")
    if tile_gdf.crs != image_crs:
        print(f"Reprojecting tiles from {tile_gdf.crs} to image CRS {image_crs}")
        tile_gdf = tile_gdf.to_crs(image_crs)

    tile_gdf = _buffer_tiles(tile_gdf, args.tile_buffer_m)
    os.makedirs(args.image_tile_output_dir, exist_ok=True)
    os.makedirs(args.lidar_tile_output_dir, exist_ok=True)
    os.makedirs(args.elastix_work_dir, exist_ok=True)

    lidar_srs = args.lidar_srs or get_srs_from_ept(args.ept_path)
    if not lidar_srs:
        raise ValueError("Unable to determine LiDAR SRS from EPT. Set --lidar-srs.")
    print(f"LiDAR SRS: {lidar_srs}")

    total_tiles = 0
    aligned_tiles = 0
    skipped_no_data = 0
    failed_tiles = 0

    for _, tile_row in tile_gdf.iterrows():
        total_tiles += 1
        tile_id = _sanitize_for_filename(tile_row[args.tile_id_field])
        image_tile_path = os.path.join(
            args.image_tile_output_dir,
            f"{tile_id}_{image_name}{args.image_tile_suffix}.tif",
        )
        intensity_tile_path = os.path.join(
            args.lidar_tile_output_dir,
            f"{tile_id}_{args.lidar_dimension}{args.lidar_tile_suffix}.tif",
        )
        aligned_path = os.path.join(
            args.image_tile_output_dir,
            f"{tile_id}_{image_name}{args.image_tile_suffix}{args.aligned_suffix}.tif",
        )

        if not os.path.exists(image_tile_path) or args.overwrite:
            print(f"Tile {tile_id}: writing image tile")
            _write_clipped_tile(args.image_path, tile_row.geometry, image_tile_path)
        else:
            print(f"Tile {tile_id}: image tile exists, skipping")
        moving_valid = _count_valid_pixels(image_tile_path, band_index=args.moving_band_index)
        if moving_valid < args.min_valid_pixels:
            print(
                f"Tile {tile_id}: moving tile has too few valid pixels "
                f"({moving_valid} < {args.min_valid_pixels}), skipping."
            )
            skipped_no_data += 1
            continue

        if not os.path.exists(intensity_tile_path) or args.overwrite:
            tile_geom_lidar = _reproject_geometry(tile_row.geometry, str(tile_gdf.crs), lidar_srs)
            min_x, min_y, max_x, max_y = tile_geom_lidar.bounds
            bounds = ([min_x, max_x], [min_y, max_y])
            print(f"Tile {tile_id}: reading EPT within bounds {bounds}")
            arrays = read_lidar(
                input_file=args.ept_path,
                srs=lidar_srs,
                bounds=bounds,
                thin_radius=None,
                hag=False,
                hag_dtm=False,
                dtm=None,
                crop_poly=False,
                poly=None,
            )
            points = arrays[0] if arrays and len(arrays) > 0 else None
            dimension_name = args.lidar_dimension
            if points is not None and points.size > 0:
                dimension_name = _resolve_dimension_name(points.dtype.names, args.lidar_dimension)
            _write_lidar_mean_raster_from_points(
                points=points,
                dimension_name=dimension_name,
                bounds_xy=(min_x, min_y, max_x, max_y),
                resolution=args.lidar_resolution,
                output_raster_path=intensity_tile_path,
                output_crs=lidar_srs,
                nodata_value=args.lidar_nodata,
            )
            print(f"Tile {tile_id}: wrote intensity raster")
        else:
            print(f"Tile {tile_id}: intensity tile exists, skipping")
        fixed_valid = _count_valid_pixels(intensity_tile_path, band_index=1)
        if fixed_valid < args.min_valid_pixels:
            print(
                f"Tile {tile_id}: intensity tile has too few valid pixels "
                f"({fixed_valid} < {args.min_valid_pixels}), skipping."
            )
            skipped_no_data += 1
            continue

        if os.path.exists(aligned_path) and not args.overwrite:
            print(f"Tile {tile_id}: aligned output exists, skipping alignment")
            aligned_tiles += 1
            continue

        work_outdir = os.path.join(
            args.elastix_work_dir,
            _sanitize_for_filename(tile_id),
            _sanitize_for_filename(os.path.splitext(os.path.basename(image_tile_path))[0]),
        )
        os.makedirs(work_outdir, exist_ok=True)
        print(f"Tile {tile_id}: running elastix alignment")
        try:
            with rasterio.open(image_tile_path) as src_moving:
                moving_res = float(abs(src_moving.transform.a))
            _, _, _, _, _, fixed_res = _get_fixed_grid_params(intensity_tile_path)
            if args.elastix_resolution == "moving":
                reg_res = moving_res
            elif args.elastix_resolution == "fixed":
                reg_res = fixed_res
            else:
                reg_res = float(args.elastix_resolution)

            _run_elastix_pair(
                elastix_script=args.elastix_script,
                fixed_path=intensity_tile_path,
                moving_path=image_tile_path,
                outdir=work_outdir,
                moving_band_index=args.moving_band_index,
                resolution_override=reg_res,
            )
            epsg, ulx, uly, lrx, lry, res = _get_fixed_grid_params(intensity_tile_path)
            aligned_multiband = _apply_transform_to_all_bands(
                moving_image_path=image_tile_path,
                work_outdir=work_outdir,
                epsg=epsg,
                ulx=ulx,
                uly=uly,
                lrx=lrx,
                lry=lry,
                res=reg_res,
            )
            shutil.copy2(aligned_multiband, aligned_path)
            print(f"Tile {tile_id}: wrote {aligned_path}")
            aligned_tiles += 1
        except Exception as exc:
            print(f"Tile {tile_id}: elastix failed, skipping tile. Error: {exc}")
            failed_tiles += 1

    print(
        "All processing complete "
        f"(tiles={total_tiles}, aligned={aligned_tiles}, "
        f"skipped_no_data={skipped_no_data}, failed={failed_tiles})."
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Tile image + LiDAR EPT mean rasters + elastix alignment per tile.",
    )
    parser.add_argument("--config-yaml", help="Optional YAML config file.")
    parser.add_argument("--image-path", help="Path to source image raster (multi-band supported).")
    parser.add_argument("--ept-path", help="Path to EPT metadata file (ept.json).")
    parser.add_argument("--tile-gpkg", help="Tile GeoPackage path.")
    parser.add_argument("--tile-layer", help="Optional tile layer name in GeoPackage.")
    parser.add_argument("--tile-id-field", default="id", help="Tile ID field.")
    parser.add_argument("--tile-buffer-m", type=float, default=100.0, help="Tile buffer in map units.")
    parser.add_argument("--image-tile-output-dir", help="Output directory for clipped image tiles.")
    parser.add_argument("--image-tile-suffix", default="_final", help="Suffix for clipped image tile names.")
    parser.add_argument("--lidar-tile-output-dir", help="Output directory for intensity tiles.")
    parser.add_argument("--lidar-dimension", default="Intensity", help="LiDAR dimension for mean raster.")
    parser.add_argument("--lidar-resolution", type=float, help="LiDAR output raster resolution.")
    parser.add_argument("--lidar-srs", help="Optional LiDAR SRS override.")
    parser.add_argument("--lidar-tile-suffix", default="_mean", help="Suffix for LiDAR tile names.")
    parser.add_argument("--lidar-nodata", type=float, default=-9999.0, help="NoData for LiDAR mean rasters.")
    parser.add_argument(
        "--elastix-script",
        default="/home/milo/repos/elastix-wrapper/run_elastix_wv3_cbl_to_intensity.sh",
        help="Path to elastix wrapper shell script.",
    )
    parser.add_argument("--elastix-work-dir", help="Directory for elastix intermediate outputs.")
    parser.add_argument("--moving-band-index", type=int, default=1, help="Band used from moving image for alignment.")
    parser.add_argument(
        "--elastix-resolution",
        default="moving",
        help=(
            "Registration/output grid resolution. Use 'moving' (default), 'fixed', "
            "or a numeric value in map units."
        ),
    )
    parser.add_argument(
        "--min-valid-pixels",
        type=int,
        default=100,
        help="Minimum valid pixels required in moving/fixed tiles to attempt elastix.",
    )
    parser.add_argument("--aligned-suffix", default="_aligned", help="Suffix for final aligned image tiles.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing tile/intensity/aligned outputs.")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument("--config-yaml")
    config_args, _ = config_parser.parse_known_args(argv)

    config_defaults: Dict = {}
    if config_args.config_yaml:
        config_defaults = _load_yaml_config(config_args.config_yaml)

    parser = build_parser()
    if config_defaults:
        parser.set_defaults(**config_defaults)
    args = parser.parse_args(argv)

    required = (
        "image_path",
        "ept_path",
        "tile_gpkg",
        "image_tile_output_dir",
        "lidar_tile_output_dir",
        "elastix_work_dir",
    )
    for key in required:
        if not getattr(args, key):
            parser.error(f"--{key.replace('_', '-')} is required (or set in --config-yaml).")
    if args.tile_buffer_m < 0:
        parser.error("--tile-buffer-m must be >= 0.")
    if args.lidar_resolution is None or args.lidar_resolution <= 0:
        parser.error("--lidar-resolution must be > 0.")
    if args.moving_band_index < 1:
        parser.error("--moving-band-index must be >= 1.")
    if args.elastix_resolution not in ("moving", "fixed"):
        try:
            if float(args.elastix_resolution) <= 0:
                raise ValueError
        except ValueError:
            parser.error("--elastix-resolution must be 'moving', 'fixed', or a positive number.")
    if args.min_valid_pixels < 0:
        parser.error("--min-valid-pixels must be >= 0.")

    return run_workflow(args)


if __name__ == "__main__":
    sys.exit(main())
