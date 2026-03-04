#!/usr/bin/env python3
"""Create a 0.5 m intensity raster from an EPT source for a reference raster bbox.

This script is intended to run in an environment where `pyforestscan` is installed.
"""

from __future__ import annotations

import argparse
from typing import Optional, Tuple

import numpy as np
import rasterio
from pyproj import Transformer
from scipy.interpolate import griddata

from pyforestscan.calculate import calculate_voxel_stat
from pyforestscan.handlers import create_geotiff, read_lidar


def _parse_z_index_range(start: Optional[int], stop: Optional[int]) -> Optional[Tuple[int, Optional[int]]]:
    if start is None and stop is None:
        return None
    if start is None:
        start = 0
    if start < 0:
        raise ValueError("--z-index-start must be >= 0")
    if stop is not None and stop <= start:
        raise ValueError("--z-index-stop must be > --z-index-start")
    return (start, stop)


def _transform_bounds(left: float, bottom: float, right: float, top: float, src_crs: str, dst_crs: str):
    if src_crs == dst_crs:
        return left, bottom, right, top

    transformer = Transformer.from_crs(src_crs, dst_crs, always_xy=True)
    xs, ys = transformer.transform(
        [left, right, right, left],
        [bottom, bottom, top, top],
    )
    return min(xs), min(ys), max(xs), max(ys)


def _interpolate_nans(arr: np.ndarray, method: str) -> np.ndarray:
    if method == "none":
        return arr

    valid = np.isfinite(arr)
    if not np.any(valid):
        return arr
    if np.all(valid):
        return arr

    rows, cols = np.indices(arr.shape)
    valid_points = np.column_stack([cols[valid], rows[valid]])
    valid_values = arr[valid]
    missing_points = np.column_stack([cols[~valid], rows[~valid]])

    interp_method = method
    filled = griddata(valid_points, valid_values, missing_points, method=interp_method)

    out = arr.copy()
    out[~valid] = filled

    # If linear/cubic leaves edge NaNs, finish with nearest.
    still_nan = ~np.isfinite(out)
    if np.any(still_nan):
        fallback = griddata(valid_points, valid_values, np.column_stack([cols[still_nan], rows[still_nan]]), method="nearest")
        out[still_nan] = fallback

    return out


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Create intensity raster from EPT using PyForestScan voxel stats.")
    parser.add_argument(
        "--ept-file",
        default="/mnt/x/PROJECTS_2/Big_Island/ChangeHI_Trees/Dry_Forest/Data/Lidar/ept-full/ept.json",
        help="Path to input EPT (ept.json).",
    )
    parser.add_argument("--reference-raster", required=True, help="Raster whose bbox is used for clipping.")
    parser.add_argument("--output-raster", required=True, help="Output GeoTIFF path.")
    parser.add_argument(
        "--lidar-srs",
        default="EPSG:6635",
        help="CRS of the EPT source, used by pyforestscan read_lidar (default: EPSG:6635).",
    )
    parser.add_argument("--resolution", type=float, default=0.5, help="XY resolution in map units (default: 0.5).")
    parser.add_argument("--voxel-height", type=float, default=1.0, help="Z voxel size for binning (default: 1.0).")
    parser.add_argument(
        "--dimension",
        default="Intensity",
        help="Point attribute/dimension to summarize (default: Intensity).",
    )
    parser.add_argument(
        "--stat",
        default="mean",
        choices=["mean", "sum", "count", "min", "max", "median", "std"],
        help="Statistic for voxel columns (default: mean).",
    )
    parser.add_argument("--z-index-start", type=int, help="Optional start z-index (inclusive).")
    parser.add_argument("--z-index-stop", type=int, help="Optional stop z-index (exclusive).")
    parser.add_argument(
        "--interpolation",
        default="none",
        choices=["none", "nearest", "linear", "cubic"],
        help="Interpolation method for NaN gaps after voxel stat (default: none).",
    )
    parser.add_argument("--nodata", type=float, default=-9999, help="NoData value for output raster.")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.resolution <= 0:
        parser.error("--resolution must be > 0")
    if args.voxel_height <= 0:
        parser.error("--voxel-height must be > 0")

    z_index_range = _parse_z_index_range(args.z_index_start, args.z_index_stop)

    with rasterio.open(args.reference_raster) as src:
        if src.crs is None:
            raise ValueError(f"Reference raster has no CRS: {args.reference_raster}")
        ref_crs = src.crs.to_string()
        left, bottom, right, top = src.bounds

    clip_left, clip_bottom, clip_right, clip_top = _transform_bounds(
        left,
        bottom,
        right,
        top,
        src_crs=ref_crs,
        dst_crs=args.lidar_srs,
    )

    bounds = ([clip_left, clip_right], [clip_bottom, clip_top])
    arrays = read_lidar(args.ept_file, args.lidar_srs, bounds=bounds)
    if not arrays or arrays[0].size == 0:
        raise RuntimeError("No lidar points returned for requested bounds.")

    pts = arrays[0]
    layer, extent = calculate_voxel_stat(
        pts,
        voxel_resolution=(args.resolution, args.resolution, args.voxel_height),
        dimension=args.dimension,
        stat=args.stat,
        z_index_range=z_index_range,
    )

    layer = _interpolate_nans(layer.astype(np.float32), method=args.interpolation)
    layer = np.where(np.isfinite(layer), layer, float(args.nodata)).astype(np.float32)

    create_geotiff(layer, args.output_raster, args.lidar_srs, extent, nodata=float(args.nodata))

    print("Created raster:", args.output_raster)
    print("Reference raster:", args.reference_raster)
    print("Reference CRS:", ref_crs)
    print("Lidar CRS:", args.lidar_srs)
    print("Clip bounds in lidar CRS:", (clip_left, clip_bottom, clip_right, clip_top))
    print("Output extent:", extent)
    print("Metric:", f"{args.stat}({args.dimension})")
    print("Resolution:", args.resolution)
    print("Interpolation:", args.interpolation)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
