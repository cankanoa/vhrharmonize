#!/usr/bin/env python3
"""Clip a raster to tiles and run OmniCloudMask per tile."""

import argparse
import importlib.util
import json
import os
import sys
from typing import Dict, List, Optional

import geopandas as gpd
import rasterio
import rasterio.mask


def _load_cloudmask_functions():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    module_path = os.path.join(os.path.dirname(script_dir), "vhrharmonize", "cloudmask.py")
    spec = importlib.util.spec_from_file_location("vhr_cloudmask", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load cloudmask module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.create_cloud_mask_with_omnicloudmask, module.apply_binary_cloud_mask_to_image


def _parse_json_dict(raw_json: Optional[str]) -> Dict:
    if not raw_json:
        return {}
    parsed = json.loads(raw_json)
    if not isinstance(parsed, dict):
        raise ValueError("Expected a JSON object for omnicloud kwargs.")
    return parsed


def _sanitize_for_filename(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in str(value))


def _parse_ids(raw_values: Optional[List[str]]) -> List[str]:
    if not raw_values:
        return []
    out: List[str] = []
    for raw in raw_values:
        for v in raw.split(","):
            v = v.strip()
            if v:
                out.append(v)
    return out


def _clip_tile(input_raster: str, geom, output_path: str) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with rasterio.open(input_raster) as src:
        data, transform = rasterio.mask.mask(
            src,
            [geom.__geo_interface__],
            crop=True,
            nodata=src.nodata,
        )
        profile = src.profile.copy()
        profile.update(height=data.shape[1], width=data.shape[2], transform=transform)
        with rasterio.open(output_path, "w", **profile) as dst:
            dst.write(data)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Run cloud masking on raster tiles only.")
    parser.add_argument("--input-raster", required=True, help="Large input raster (e.g. pansharpened tif).")
    parser.add_argument("--tile-gpkg", required=True, help="Tile GeoPackage.")
    parser.add_argument("--tile-layer", help="Optional tile layer name.")
    parser.add_argument("--tile-id-field", default="id", help="Tile ID field in tile layer.")
    parser.add_argument("--tile-buffer-m", type=float, default=0.0, help="Optional tile buffer in map units.")
    parser.add_argument("--output-dir", required=True, help="Output directory for tile products.")
    parser.add_argument("--tile-suffix", default="_final", help="Suffix for clipped tile rasters.")
    parser.add_argument("--masked-suffix", default="_cloudmasked", help="Suffix for masked tile rasters.")
    parser.add_argument("--mask-suffix", default="_cloudmask", help="Suffix for binary mask rasters.")
    parser.add_argument("--nodata-value", type=float, default=-9999, help="NoData for masked output.")
    parser.add_argument("--red-band-index", type=int, default=5)
    parser.add_argument("--green-band-index", type=int, default=3)
    parser.add_argument("--nir-band-index", type=int, default=7)
    parser.add_argument("--cloud-classes", default="1,2,3")
    parser.add_argument("--buffer-pixels", type=int, default=0)
    parser.add_argument("--omnicloud-kwargs-json")
    parser.add_argument("--filter-tile-id", action="append", help="Optional tile ids (repeat or comma-separated).")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs.")
    args = parser.parse_args(argv)

    if args.tile_buffer_m < 0:
        parser.error("--tile-buffer-m must be >= 0.")
    if args.buffer_pixels < 0:
        parser.error("--buffer-pixels must be >= 0.")
    if not os.path.isfile(args.input_raster):
        parser.error(f"--input-raster does not exist: {args.input_raster}")

    tile_ids_filter = set(_parse_ids(args.filter_tile_id))
    cloud_classes = [int(x.strip()) for x in args.cloud_classes.split(",") if x.strip()]
    omnicloud_kwargs = _parse_json_dict(args.omnicloud_kwargs_json)
    create_cloud_mask_with_omnicloudmask, apply_binary_cloud_mask_to_image = _load_cloudmask_functions()

    tiles = gpd.read_file(args.tile_gpkg, layer=args.tile_layer)
    if tiles.empty:
        raise ValueError(f"Tile layer empty: {args.tile_gpkg}")
    if args.tile_id_field not in tiles.columns:
        raise ValueError(f"Missing tile id field '{args.tile_id_field}' in tile layer.")
    if tiles.crs is None:
        raise ValueError("Tile layer CRS is undefined.")

    with rasterio.open(args.input_raster) as src:
        raster_crs = src.crs
    if raster_crs is None:
        raise ValueError("Input raster CRS is undefined.")
    if tiles.crs != raster_crs:
        tiles = tiles.to_crs(raster_crs)
    if args.tile_buffer_m != 0:
        tiles = tiles.copy()
        tiles["geometry"] = tiles.geometry.buffer(args.tile_buffer_m)

    os.makedirs(args.output_dir, exist_ok=True)
    total = 0
    done = 0
    skipped = 0
    for _, row in tiles.iterrows():
        tile_id = _sanitize_for_filename(row[args.tile_id_field])
        if tile_ids_filter and tile_id not in tile_ids_filter:
            continue
        total += 1

        tile_raster = os.path.join(args.output_dir, f"{tile_id}{args.tile_suffix}.tif")
        mask_raster = os.path.join(args.output_dir, f"{tile_id}{args.mask_suffix}.tif")
        masked_raster = os.path.join(args.output_dir, f"{tile_id}{args.masked_suffix}.tif")
        if (not args.overwrite) and os.path.exists(masked_raster):
            print(f"Tile {tile_id}: masked output exists, skipping.")
            skipped += 1
            continue

        try:
            print(f"Tile {tile_id}: clipping")
            _clip_tile(args.input_raster, row.geometry, tile_raster)
            print(f"Tile {tile_id}: cloudmask")
            create_cloud_mask_with_omnicloudmask(
                tile_raster,
                mask_raster,
                red_band_index=args.red_band_index,
                green_band_index=args.green_band_index,
                nir_band_index=args.nir_band_index,
                cloud_classes=cloud_classes,
                buffer_pixels=args.buffer_pixels,
                omnicloud_kwargs=omnicloud_kwargs,
            )
            apply_binary_cloud_mask_to_image(
                tile_raster,
                mask_raster,
                masked_raster,
                output_nodata_value=args.nodata_value,
            )
            done += 1
        except Exception as exc:
            print(f"Tile {tile_id}: failed: {exc}")

    print(f"Done. tiles={total}, masked={done}, skipped={skipped}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
