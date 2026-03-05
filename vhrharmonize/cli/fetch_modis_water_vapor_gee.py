#!/usr/bin/env python3
"""Fetch per-scene MODIS water vapor from Google Earth Engine."""

import argparse
import csv
import json
import os
import re
from datetime import datetime, timedelta, timezone
from typing import Dict, List, Optional

from vhrharmonize.providers.worldview.files import find_files, find_roots

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(MODULE_DIR))
DEFAULT_ENV_FILE = os.path.join(PROJECT_ROOT, "configs", ".env")


MONTH_MAP = {
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


def _load_env_file(env_file_path: str) -> Dict[str, str]:
    values: Dict[str, str] = {}
    if not env_file_path or not os.path.exists(env_file_path):
        return values
    with open(env_file_path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip().strip('"').strip("'")
            if key:
                values[key] = value
    return values


def _resolve_ee_project(cli_value: Optional[str], env_values: Dict[str, str]) -> Optional[str]:
    if cli_value:
        return cli_value
    for key in ("GEE_PROJECT", "EE_PROJECT", "GOOGLE_EARTH_ENGINE_PROJECT"):
        if env_values.get(key):
            return env_values[key]
    return None


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


def _parse_scene_datetime_utc(photo_basename: str) -> datetime:
    # Example: 17OCT19211717-M1BS-...
    match = re.match(r"^(\d{2})([A-Z]{3})(\d{2})(\d{2})(\d{2})(\d{2})-", photo_basename)
    if not match:
        raise ValueError(f"Could not parse acquisition datetime from: {photo_basename}")
    day = int(match.group(1))
    month = MONTH_MAP[match.group(2)]
    year = 2000 + int(match.group(3))
    hour = int(match.group(4))
    minute = int(match.group(5))
    second = int(match.group(6))
    return datetime(year, month, day, hour, minute, second, tzinfo=timezone.utc)


def _load_scene_geometry_wgs84(shp_path: str):
    import geopandas as gpd

    gdf = gpd.read_file(shp_path)
    if gdf.empty:
        raise ValueError(f"Empty scene footprint: {shp_path}")
    if gdf.crs is None:
        raise ValueError(f"Scene footprint has no CRS: {shp_path}")
    gdf = gdf.to_crs(epsg=4326)
    return gdf.geometry.union_all()


def _init_ee(project: Optional[str], authenticate: bool):
    import ee

    if authenticate:
        ee.Authenticate()
    if project:
        ee.Initialize(project=project)
    else:
        ee.Initialize()
    return ee


def _fetch_collection_value(
    ee,
    collection_id: str,
    band_name: str,
    geom_wgs84,
    scene_dt_utc: datetime,
    hours_window: int,
    reduce_scale_m: int,
):
    start = (scene_dt_utc - timedelta(hours=hours_window)).strftime("%Y-%m-%dT%H:%M:%S")
    end = (scene_dt_utc + timedelta(hours=hours_window)).strftime("%Y-%m-%dT%H:%M:%S")

    geom_json = geom_wgs84.__geo_interface__
    ee_geom = ee.Geometry(geom_json)
    scene_ms = int(scene_dt_utc.timestamp() * 1000)

    col = ee.ImageCollection(collection_id).filterDate(start, end).filterBounds(ee_geom)
    size = col.size().getInfo()
    if size == 0:
        return None

    def _add_time_diff(img):
        diff = ee.Number(img.get("system:time_start")).subtract(scene_ms).abs()
        return img.set("abs_time_diff", diff)

    best = col.map(_add_time_diff).sort("abs_time_diff").first()
    best_info = best.getInfo()
    props = best_info.get("properties", {})
    bands = [b.get("id") for b in best_info.get("bands", [])]
    if band_name not in bands:
        return {
            "collection": collection_id,
            "band_found": False,
            "bands": bands,
            "raw_value": None,
            "image_time_utc": None,
            "abs_time_diff_hours": None,
        }

    stats = best.select(band_name).reduceRegion(
        reducer=ee.Reducer.median(),
        geometry=ee_geom,
        scale=reduce_scale_m,
        bestEffort=True,
        maxPixels=1_000_000_000,
    ).getInfo()

    raw_value = stats.get(band_name) if isinstance(stats, dict) else None
    time_start_ms = props.get("system:time_start")
    abs_diff_ms = abs(int(time_start_ms) - scene_ms) if time_start_ms is not None else None
    image_time_utc = (
        datetime.fromtimestamp(time_start_ms / 1000, tz=timezone.utc).isoformat()
        if time_start_ms is not None
        else None
    )

    return {
        "collection": collection_id,
        "band_found": True,
        "bands": bands,
        "raw_value": raw_value,
        "image_time_utc": image_time_utc,
        "abs_time_diff_hours": (abs_diff_ms / 3_600_000) if abs_diff_ms is not None else None,
    }


def _fetch_first_available_band_value(
    ee,
    collection_id: str,
    band_candidates: List[str],
    geom_wgs84,
    scene_dt_utc: datetime,
    hours_window: int,
    reduce_scale_m: int,
):
    start = (scene_dt_utc - timedelta(hours=hours_window)).strftime("%Y-%m-%dT%H:%M:%S")
    end = (scene_dt_utc + timedelta(hours=hours_window)).strftime("%Y-%m-%dT%H:%M:%S")

    geom_json = geom_wgs84.__geo_interface__
    ee_geom = ee.Geometry(geom_json)
    scene_ms = int(scene_dt_utc.timestamp() * 1000)

    col = ee.ImageCollection(collection_id).filterDate(start, end).filterBounds(ee_geom)
    size = col.size().getInfo()
    if size == 0:
        return None

    def _add_time_diff(img):
        diff = ee.Number(img.get("system:time_start")).subtract(scene_ms).abs()
        return img.set("abs_time_diff", diff)

    best = col.map(_add_time_diff).sort("abs_time_diff").first()
    best_info = best.getInfo()
    props = best_info.get("properties", {})
    bands = [b.get("id") for b in best_info.get("bands", [])]

    selected_band = None
    for band in band_candidates:
        if band in bands:
            selected_band = band
            break
    if selected_band is None:
        return {
            "collection": collection_id,
            "band_found": False,
            "selected_band": None,
            "bands": bands,
            "raw_value": None,
            "image_time_utc": None,
            "abs_time_diff_hours": None,
        }

    stats = best.select(selected_band).reduceRegion(
        reducer=ee.Reducer.median(),
        geometry=ee_geom,
        scale=reduce_scale_m,
        bestEffort=True,
        maxPixels=1_000_000_000,
    ).getInfo()

    raw_value = stats.get(selected_band) if isinstance(stats, dict) else None
    time_start_ms = props.get("system:time_start")
    abs_diff_ms = abs(int(time_start_ms) - scene_ms) if time_start_ms is not None else None
    image_time_utc = (
        datetime.fromtimestamp(time_start_ms / 1000, tz=timezone.utc).isoformat()
        if time_start_ms is not None
        else None
    )

    return {
        "collection": collection_id,
        "band_found": True,
        "selected_band": selected_band,
        "bands": bands,
        "raw_value": raw_value,
        "image_time_utc": image_time_utc,
        "abs_time_diff_hours": (abs_diff_ms / 3_600_000) if abs_diff_ms is not None else None,
    }


def _choose_best_candidate(candidates: List[Dict]) -> Optional[Dict]:
    valid = [
        c
        for c in candidates
        if c
        and c.get("band_found")
        and c.get("raw_value") is not None
        and c.get("abs_time_diff_hours") is not None
    ]
    if not valid:
        return None
    return sorted(valid, key=lambda c: c["abs_time_diff_hours"])[0]


def _choose_modtran_atm_by_lat_month(lat_deg: float, month: int) -> str:
    if abs(lat_deg) <= 23.5:
        return "Tropical Atmosphere"
    if month in (12, 1, 2):
        return "Mid-Latitude Winter"
    return "Mid-Latitude Summer"


def _visibility_from_aod(aod: float) -> float:
    # Conservative lookup to avoid extreme swings in FLAASH behavior.
    if aod < 0.10:
        return 30.0
    if aod < 0.20:
        return 20.0
    if aod < 0.35:
        return 12.0
    if aod < 0.60:
        return 8.0
    return 5.0


def _update_root_file_per_photo(root_file_path: str, per_photo_updates: Dict[str, Dict]):
    if os.path.exists(root_file_path):
        with open(root_file_path, "r", encoding="utf-8") as f:
            content = f.read().strip()
        data = json.loads(content) if content else {}
    else:
        data = {}

    if not isinstance(data, dict):
        data = {}

    params_per_photo = data.get("ParamsOverridesPerPhoto")
    if not isinstance(params_per_photo, dict):
        params_per_photo = {}
        data["ParamsOverridesPerPhoto"] = params_per_photo

    for photo_basename, overrides in per_photo_updates.items():
        existing = params_per_photo.get(photo_basename)
        if not isinstance(existing, dict):
            existing = {}
        existing.update(overrides)
        params_per_photo[photo_basename] = existing

    with open(root_file_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser for MODIS water vapor and atmosphere override extraction."""
    parser = argparse.ArgumentParser(
        description="Fetch MODIS water vapor per scene from Google Earth Engine and compute FLAASH preset.",
    )
    parser.add_argument(
        "--input-dir",
        action="append",
        required=True,
        help="Input directory to scan for Root files and scenes. Repeat for multiple directories.",
    )
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
    parser.add_argument(
        "--authenticate",
        action="store_true",
        help="Run interactive ee.Authenticate() before ee.Initialize().",
    )
    parser.add_argument(
        "--hours-window",
        type=int,
        default=24,
        help="Hours around scene timestamp to search MODIS granules.",
    )
    parser.add_argument(
        "--terra-collection",
        default="MODIS/061/MOD08_D3",
        help="Earth Engine collection id for Terra water vapor.",
    )
    parser.add_argument(
        "--aqua-collection",
        default="MODIS/061/MYD08_D3",
        help="Earth Engine collection id for Aqua water vapor.",
    )
    parser.add_argument(
        "--band-name",
        default="Atmospheric_Water_Vapor_Mean",
        help="Band name to reduce from MODIS product.",
    )
    parser.add_argument(
        "--aod-band-candidates",
        default=(
            "Aerosol_Optical_Depth_Land_Ocean_Mean,"
            "Aerosol_Optical_Depth_Land_Ocean_Mean_Mean,"
            "Aerosol_Optical_Depth_Land_Ocean"
        ),
        help="Comma-separated candidate AOD band names; first available is used.",
    )
    parser.add_argument(
        "--aod-scale-factor",
        type=float,
        default=0.001,
        help="Scale factor applied to AOD band values before visibility heuristic.",
    )
    parser.add_argument(
        "--modis-scale-factor",
        type=float,
        default=0.001,
        help="Scale factor applied to reduced MODIS band value to obtain g/cm^2 (or cm PW).",
    )
    parser.add_argument(
        "--modtran-baseline-water-vapor",
        type=float,
        default=2.92,
        help="Baseline water vapor (g/cm^2) used to compute WATER_VAPOR_PRESET ratio.",
    )
    parser.add_argument(
        "--reduce-scale-m",
        type=int,
        default=1000,
        help="Scale (meters) for reduceRegion median.",
    )
    parser.add_argument(
        "--output-json",
        default="outputs/modis_water_vapor_results.json",
        help="Output JSON report path.",
    )
    parser.add_argument(
        "--output-csv",
        default="outputs/modis_water_vapor_results.csv",
        help="Output CSV report path.",
    )
    parser.add_argument(
        "--write-root-overrides",
        action="store_true",
        help="Write computed WATER_VAPOR_PRESET into Root_WV.txt under ParamsOverridesPerPhoto.",
    )
    parser.add_argument(
        "--write-atmosphere-overrides",
        action="store_true",
        help=(
            "Also write per-photo atmospheric overrides using local scene context and MODIS AOD: "
            "MODTRAN_ATM, USE_AEROSOL, DEFAULT_VISIBILITY."
        ),
    )
    return parser


def main() -> int:
    """Run MODIS/AOD retrieval workflow and write JSON/CSV outputs (and optional root updates)."""
    args = build_parser().parse_args()
    env_values = _load_env_file(args.env_file)
    ee_project = _resolve_ee_project(args.ee_project, env_values)
    if not ee_project:
        print(
            "Warning: no Earth Engine project id provided. "
            "Set --ee-project or add GEE_PROJECT in configs/.env."
        )
    ee = _init_ee(ee_project, args.authenticate)
    filter_basenames = _parse_filter_basenames(args.filter_basename)
    aod_band_candidates = [x.strip() for x in args.aod_band_candidates.split(",") if x.strip()]

    rows = []
    updates_by_root: Dict[str, Dict[str, Dict]] = {}

    for input_folder in args.input_dir:
        root_paths = find_roots(input_folder)
        for root_folder_path, root_file_path in root_paths:
            found = find_files(root_folder_path, root_file_path, filter_basenames)
            for _, scene in found.items():
                mul_photo_basename = scene.get("mul_photo_basename")
                mul_shp = scene.get("mul_shp_file")
                if not mul_photo_basename or not mul_shp:
                    continue

                print(f"Fetching MODIS water vapor for {mul_photo_basename}")
                try:
                    scene_dt = _parse_scene_datetime_utc(mul_photo_basename)
                    geom = _load_scene_geometry_wgs84(mul_shp)
                    centroid_lat = float(geom.centroid.y)

                    terra = _fetch_collection_value(
                        ee,
                        args.terra_collection,
                        args.band_name,
                        geom,
                        scene_dt,
                        args.hours_window,
                        args.reduce_scale_m,
                    )
                    aqua = _fetch_collection_value(
                        ee,
                        args.aqua_collection,
                        args.band_name,
                        geom,
                        scene_dt,
                        args.hours_window,
                        args.reduce_scale_m,
                    )

                    best = _choose_best_candidate([terra, aqua])
                    if best is None:
                        rows.append(
                            {
                                "input_folder": input_folder,
                                "root_file_path": root_file_path,
                                "photo_basename": mul_photo_basename,
                                "scene_time_utc": scene_dt.isoformat(),
                                "status": "no_modis_match",
                                "selected_collection": None,
                                "modis_raw_value": None,
                                "modis_scaled_g_cm2": None,
                                "water_vapor_preset": None,
                                "image_time_utc": None,
                                "abs_time_diff_hours": None,
                            }
                        )
                        continue

                    raw_value = float(best["raw_value"])
                    scaled = raw_value * args.modis_scale_factor
                    preset = scaled / args.modtran_baseline_water_vapor

                    terra_aod = _fetch_first_available_band_value(
                        ee,
                        args.terra_collection,
                        aod_band_candidates,
                        geom,
                        scene_dt,
                        args.hours_window,
                        args.reduce_scale_m,
                    )
                    aqua_aod = _fetch_first_available_band_value(
                        ee,
                        args.aqua_collection,
                        aod_band_candidates,
                        geom,
                        scene_dt,
                        args.hours_window,
                        args.reduce_scale_m,
                    )
                    best_aod = _choose_best_candidate([terra_aod, aqua_aod])
                    aod_raw = float(best_aod["raw_value"]) if best_aod and best_aod.get("raw_value") is not None else None
                    aod_scaled = aod_raw * args.aod_scale_factor if aod_raw is not None else None
                    default_visibility = _visibility_from_aod(aod_scaled) if aod_scaled is not None else None
                    modtran_atm = _choose_modtran_atm_by_lat_month(centroid_lat, scene_dt.month)

                    row = {
                        "input_folder": input_folder,
                        "root_file_path": root_file_path,
                        "photo_basename": mul_photo_basename,
                        "scene_time_utc": scene_dt.isoformat(),
                        "status": "ok",
                        "selected_collection": best["collection"],
                        "modis_raw_value": raw_value,
                        "modis_scaled_g_cm2": scaled,
                        "water_vapor_preset": preset,
                        "aod_raw": aod_raw,
                        "aod_scaled": aod_scaled,
                        "default_visibility": default_visibility,
                        "modtran_atm": modtran_atm,
                        "image_time_utc": best["image_time_utc"],
                        "abs_time_diff_hours": best["abs_time_diff_hours"],
                    }
                    rows.append(row)

                    if args.write_root_overrides or args.write_atmosphere_overrides:
                        updates_by_root.setdefault(root_file_path, {})
                        update_values = {}
                        if args.write_root_overrides:
                            update_values["WATER_VAPOR_PRESET"] = round(preset, 4)
                        if args.write_atmosphere_overrides:
                            update_values.update(
                                {
                                    "MODTRAN_ATM": modtran_atm,
                                    "USE_AEROSOL": "Automatic Selection",
                                }
                            )
                            if default_visibility is not None:
                                update_values["DEFAULT_VISIBILITY"] = round(default_visibility, 2)
                        updates_by_root[root_file_path][mul_photo_basename] = update_values

                except Exception as exc:
                    rows.append(
                        {
                            "input_folder": input_folder,
                            "root_file_path": root_file_path,
                            "photo_basename": mul_photo_basename,
                            "scene_time_utc": None,
                            "status": "error",
                            "error": str(exc),
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
                        }
                    )

    os.makedirs(os.path.dirname(args.output_json) or ".", exist_ok=True)
    with open(args.output_json, "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2)

    os.makedirs(os.path.dirname(args.output_csv) or ".", exist_ok=True)
    fieldnames = [
        "input_folder",
        "root_file_path",
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
    with open(args.output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k) for k in fieldnames})

    if args.write_root_overrides or args.write_atmosphere_overrides:
        for root_file_path, per_photo_updates in updates_by_root.items():
            print(f"Updating {root_file_path} with {len(per_photo_updates)} per-photo override entries")
            _update_root_file_per_photo(root_file_path, per_photo_updates)

    ok_count = sum(1 for r in rows if r.get("status") == "ok")
    print(f"Done. Scenes processed: {len(rows)}; successful MODIS matches: {ok_count}")
    print(f"JSON report: {args.output_json}")
    print(f"CSV report: {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
