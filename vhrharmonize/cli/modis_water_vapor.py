"""Canonical MODIS water vapor CLI entrypoint."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from typing import Dict, List, Optional

from vhrharmonize.preprocess.fetch_external_data import (
    fetch_modis_water_vapor_for_bbox,
    init_ee_client,
)
from vhrharmonize.providers.worldview import find_files, parse_worldview_basename

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(MODULE_DIR))
DEFAULT_ENV_FILE = os.path.join(PROJECT_ROOT, "configs", ".env")


def _parse_filter_basenames(raw_values: Optional[List[str]]) -> List[str]:
    """Normalize optional basename filters.
    Args:
        raw_values: Repeated or comma-delimited basename filters.
    Returns:
        Flattened basename filter list.
    """
    if not raw_values:
        return []
    parsed: List[str] = []
    for raw in raw_values:
        for value in raw.split(","):
            value = value.strip()
            if value:
                parsed.append(value)
    return parsed


def _parse_scene_datetime_utc(photo_basename: str) -> object:
    """Parse a scene timestamp from a basename.
    Args:
        photo_basename: WorldView photo basename.
    Returns:
        Parsed acquisition datetime object.
    """
    parts = parse_worldview_basename(photo_basename)
    if parts is None:
        raise ValueError(f"Could not parse acquisition datetime from: {photo_basename}")
    return parts.acquisition_datetime_utc


def _load_scene_bbox_wgs84(shp_path: str) -> tuple[float, float, float, float]:
    """Load a scene footprint bounding box in WGS84.
    Args:
        shp_path: Path to the footprint shapefile.
    Returns:
        Bounding box as min lon, min lat, max lon, max lat.
    """
    import geopandas as gpd

    gdf = gpd.read_file(shp_path)
    if gdf.empty:
        raise ValueError(f"Empty scene footprint: {shp_path}")
    if gdf.crs is None:
        raise ValueError(f"Scene footprint has no CRS: {shp_path}")
    gdf = gdf.to_crs(epsg=4326)
    minx, miny, maxx, maxy = gdf.total_bounds
    return float(minx), float(miny), float(maxx), float(maxy)


def _build_parser() -> argparse.ArgumentParser:
    """Build the MODIS water vapor CLI parser.
    Args:
        None.
    Returns:
        Configured argument parser.
    """
    parser = argparse.ArgumentParser(
        description="Fetch MODIS water vapor per scene from Google Earth Engine and write reports.",
    )
    parser.add_argument("--input-dir", action="append", required=True, help="Input directory to scan for scenes.")
    parser.add_argument(
        "--filter-basename",
        action="append",
        help="Optional image basenames to process. Repeat or pass comma-separated values.",
    )
    parser.add_argument(
        "--env-file",
        default=DEFAULT_ENV_FILE,
        help="Optional .env file for Earth Engine settings (default: configs/.env).",
    )
    parser.add_argument("--ee-project", help="Earth Engine project id (recommended).")
    parser.add_argument("--authenticate", action="store_true", help="Run ee.Authenticate() before ee.Initialize().")
    parser.add_argument("--hours-window", type=int, default=24, help="Hours around scene time to search MODIS.")
    parser.add_argument("--terra-collection", default="MODIS/061/MOD08_D3")
    parser.add_argument("--aqua-collection", default="MODIS/061/MYD08_D3")
    parser.add_argument("--band-name", default="Atmospheric_Water_Vapor_Mean")
    parser.add_argument(
        "--aod-band-candidates",
        default=(
            "Aerosol_Optical_Depth_Land_Ocean_Mean,"
            "Aerosol_Optical_Depth_Land_Ocean_Mean_Mean,"
            "Aerosol_Optical_Depth_Land_Ocean"
        ),
    )
    parser.add_argument("--aod-scale-factor", type=float, default=0.001)
    parser.add_argument("--modis-scale-factor", type=float, default=0.001)
    parser.add_argument("--modtran-baseline-water-vapor", type=float, default=2.92)
    parser.add_argument("--reduce-scale-m", type=int, default=1000)
    parser.add_argument("--output-json", default="outputs/modis_water_vapor_results.json")
    parser.add_argument("--output-csv", default="outputs/modis_water_vapor_results.csv")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """Run the MODIS water vapor CLI.
    Args:
        argv: Optional command line arguments.
    Returns:
        Process exit code.
    """
    args = _build_parser().parse_args(argv)
    ee = init_ee_client(
        ee_project=args.ee_project,
        authenticate=args.authenticate,
        env_file=args.env_file,
    )
    filter_basenames = _parse_filter_basenames(args.filter_basename)
    aod_band_candidates = [value.strip() for value in args.aod_band_candidates.split(",") if value.strip()]

    rows: List[Dict[str, object]] = []
    for input_folder in args.input_dir:
        found = find_files(input_folder, filter_basenames)
        for _, scene in found.items():
            mul_photo_basename = scene.get("mul_photo_basename")
            mul_shp = scene.get("mul_shp_file")
            if not mul_photo_basename or not mul_shp:
                continue

            print(f"Fetching MODIS water vapor for {mul_photo_basename}")
            try:
                scene_dt = _parse_scene_datetime_utc(mul_photo_basename)
                min_lon, min_lat, max_lon, max_lat = _load_scene_bbox_wgs84(mul_shp)
                estimate = fetch_modis_water_vapor_for_bbox(
                    scene_datetime_utc=scene_dt,
                    min_lon=min_lon,
                    min_lat=min_lat,
                    max_lon=max_lon,
                    max_lat=max_lat,
                    ee=ee,
                    hours_window=args.hours_window,
                    terra_collection=args.terra_collection,
                    aqua_collection=args.aqua_collection,
                    band_name=args.band_name,
                    aod_band_candidates=aod_band_candidates,
                    aod_scale_factor=args.aod_scale_factor,
                    modis_scale_factor=args.modis_scale_factor,
                    modtran_baseline_water_vapor=args.modtran_baseline_water_vapor,
                    reduce_scale_m=args.reduce_scale_m,
                )
                row = {
                    "input_folder": input_folder,
                    "photo_basename": mul_photo_basename,
                    "scene_time_utc": scene_dt.isoformat(),
                    **estimate.to_dict(),
                }
                rows.append(row)
            except Exception as exc:
                rows.append(
                    {
                        "input_folder": input_folder,
                        "photo_basename": mul_photo_basename,
                        "scene_time_utc": None,
                        "status": "error",
                        "selected_collection": None,
                        "modis_raw_value": None,
                        "modis_scaled_g_cm2": None,
                        "water_vapor_preset": None,
                        "aod_raw": None,
                        "aod_scaled": None,
                        "default_visibility": None,
                        "modtran_atm": None,
                        "image_time_utc": None,
                        "abs_time_diff_hours": None,
                        "error": str(exc),
                    }
                )

    os.makedirs(os.path.dirname(args.output_json) or ".", exist_ok=True)
    with open(args.output_json, "w", encoding="utf-8") as handle:
        json.dump(rows, handle, indent=2)

    os.makedirs(os.path.dirname(args.output_csv) or ".", exist_ok=True)
    fieldnames = [
        "input_folder",
        "photo_basename",
        "scene_time_utc",
        "status",
        "selected_collection",
        "modis_raw_value",
        "modis_scaled_g_cm2",
        "water_vapor_preset",
        "aod_raw",
        "aod_scaled",
        "default_visibility",
        "modtran_atm",
        "image_time_utc",
        "abs_time_diff_hours",
        "error",
    ]
    with open(args.output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in fieldnames})

    ok_count = sum(1 for row in rows if row.get("status") == "ok")
    print(f"Done. Scenes processed: {len(rows)}; successful MODIS matches: {ok_count}")
    print(f"JSON report: {args.output_json}")
    print(f"CSV report: {args.output_csv}")
    return 0


__all__ = ["main"]


if __name__ == "__main__":
    sys.exit(main())
