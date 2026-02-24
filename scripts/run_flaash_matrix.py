#!/usr/bin/env python3
"""Run a minimal FLAASH aerosol/atmosphere matrix on one scene."""

import argparse
import csv
import os
import re
from typing import Dict, List, Optional

import numpy as np
import rasterio

import envipyengine.config
from envipyengine import Engine

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
if PROJECT_ROOT not in os.sys.path:
    os.sys.path.insert(0, PROJECT_ROOT)

from vhrharmonize import find_files, find_roots, get_metadata_from_files, get_image_percentile_value, run_flaash, shp_to_gpkg


MATRIX = [
    # Initial workflow default pattern
    {"name": "baseline_initial", "MODTRAN_ATM": "Mid-Latitude Summer", "MODTRAN_AER": "Maritime", "USE_AEROSOL": "Disabled", "DEFAULT_VISIBILITY": None},
    # Keep maritime, enable retrieval
    {"name": "maritime_auto_30km", "MODTRAN_ATM": "Mid-Latitude Summer", "MODTRAN_AER": "Maritime", "USE_AEROSOL": "Automatic Selection", "DEFAULT_VISIBILITY": 30.0},
    # Rural families (disabled retrieval)
    {"name": "highvis_rural_disabled", "MODTRAN_ATM": "Mid-Latitude Summer", "MODTRAN_AER": "High-Visibility Rural", "USE_AEROSOL": "Disabled", "DEFAULT_VISIBILITY": None},
    {"name": "lowvis_rural_disabled", "MODTRAN_ATM": "Mid-Latitude Summer", "MODTRAN_AER": "Low-Visibility Rural", "USE_AEROSOL": "Disabled", "DEFAULT_VISIBILITY": None},
    # Broader atmosphere family checks
    {"name": "tropical_maritime_disabled", "MODTRAN_ATM": "Tropical Atmosphere", "MODTRAN_AER": "Maritime", "USE_AEROSOL": "Disabled", "DEFAULT_VISIBILITY": None},
    {"name": "tropical_tropo_auto_20km", "MODTRAN_ATM": "Tropical Atmosphere", "MODTRAN_AER": "Tropospheric", "USE_AEROSOL": "Automatic Selection", "DEFAULT_VISIBILITY": 20.0},
]


def _wsl_to_windows(path: str) -> str:
    m = re.match(r"^/mnt/([a-zA-Z])/(.*)$", path)
    if not m:
        raise ValueError(f"Expected /mnt/<drive>/ path for Windows ENVI: {path}")
    drive = m.group(1).upper()
    rest = m.group(2).replace("/", "\\")
    return f"{drive}:\\{rest}"


def _sample_negative_percent_by_band(path: str) -> List[float]:
    with rasterio.open(path) as src:
        nodata = src.nodata
        h, w = src.height, src.width
        rows = [int(h * x) for x in (0.1, 0.5, 0.9)]
        cols = [int(w * x) for x in (0.1, 0.5, 0.9)]
        wh = min(1024, h // 4 if h > 4096 else h)
        ww = min(1024, w // 4 if w > 4096 else w)
        results: List[float] = []
        for b in range(1, src.count + 1):
            vals = []
            for r in rows:
                for c in cols:
                    r0 = max(0, min(h - wh, r - wh // 2))
                    c0 = max(0, min(w - ww, c - ww // 2))
                    arr = src.read(b, window=((r0, r0 + wh), (c0, c0 + ww)))
                    arr = arr[arr != nodata] if nodata is not None else arr.ravel()
                    if arr.size:
                        vals.append(arr)
            v = np.concatenate(vals)
            results.append(float((v < 0).sum() / v.size * 100.0))
        return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run minimal FLAASH aerosol model matrix for one scene.")
    parser.add_argument("--input-dir", required=True, help="Folder scanned for Root*.txt and scene files.")
    parser.add_argument("--photo-basename", required=True, help="Target M1BS photo basename.")
    parser.add_argument("--dem-file-path", required=True, help="DEM path (ellipsoidal heights).")
    parser.add_argument("--envi-engine-path", required=True, help="Windows ENVI taskengine.exe path.")
    parser.add_argument("--scratch-dir", required=True, help="Scratch output dir under /mnt/<drive>/...")
    parser.add_argument("--footprint-epsg", type=int, default=4326, help="Footprint EPSG assignment.")
    parser.add_argument("--dem-ground-percentile", type=float, default=100.0, help="DEM percentile for GROUND_ELEVATION.")
    parser.add_argument("--output-csv", default="scripts/flaash_matrix_results.csv", help="Results CSV path.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.dem_ground_percentile < 0 or args.dem_ground_percentile > 100:
        raise SystemExit("--dem-ground-percentile must be between 0 and 100.")
    if not re.match(r"^/mnt/[a-zA-Z]/", args.scratch_dir):
        raise SystemExit("--scratch-dir must be under /mnt/<drive>/... for Windows ENVI.")

    os.makedirs(args.scratch_dir, exist_ok=True)

    root_paths = find_roots(args.input_dir)
    if not root_paths:
        raise SystemExit("No Root* files found in input-dir.")

    target = None
    for root_folder_path, root_file_path in root_paths:
        found = find_files(root_folder_path, root_file_path, [args.photo_basename])
        for _, scene in found.items():
            if scene.get("mul_photo_basename") == args.photo_basename:
                target = (scene, root_file_path)
                break
        if target:
            break
    if not target:
        raise SystemExit(f"Target photo basename not found: {args.photo_basename}")

    scene, root_file_path = target
    mul_tif_file = scene["mul_tif_file"]
    mul_shp_path = scene["mul_shp_file"]
    mul_imd_file = scene["mul_imd_file"]
    mul_photo_basename = scene["mul_photo_basename"]

    params_overrides_scene, mul_params_overrides_photo, mul_imd_data = get_metadata_from_files(
        root_file_path, mul_imd_file, mul_photo_basename
    )
    water_vapor_preset = None
    for src in (params_overrides_scene or {}, mul_params_overrides_photo or {}):
        if "WATER_VAPOR_PRESET" in src:
            water_vapor_preset = src["WATER_VAPOR_PRESET"]

    envipyengine.config.set("engine", args.envi_engine_path)
    envi_engine = Engine("ENVI")
    envi_engine.tasks()

    rows = []
    for item in MATRIX:
        run_dir = os.path.join(args.scratch_dir, item["name"])
        os.makedirs(run_dir, exist_ok=True)

        mul_gpkg_path = os.path.join(run_dir, f"{mul_photo_basename}.gpkg")
        shp_to_gpkg(mul_shp_path, mul_gpkg_path, args.footprint_epsg)
        ground_elevation_m = get_image_percentile_value(
            args.dem_file_path, percentile=args.dem_ground_percentile, mask=mul_gpkg_path
        )

        flaash_dat = os.path.join(run_dir, f"{mul_photo_basename}_FLAASH.dat")
        flaash_params_txt = os.path.join(run_dir, f"{mul_photo_basename}_FLAASH_Params.txt")

        params: Dict = {
            "INPUT_RASTER": {"url": _wsl_to_windows(mul_tif_file), "factory": "URLRaster"},
            "MODTRAN_ATM": item["MODTRAN_ATM"],
            "MODTRAN_AER": item["MODTRAN_AER"],
            "MODTRAN_RES": 5.0,
            "MODTRAN_MSCAT": "DISORT",
            "USE_AEROSOL": item["USE_AEROSOL"],
            "DEFAULT_VISIBILITY": item["DEFAULT_VISIBILITY"],
            "AER_BAND_RATIO": 0.5,
            "AER_BANDLOW_WAVL": 425,
            "AER_BANDHIGH_WAVL": 660,
            "AER_BANDHIGH_MAXREFL": 0.2,
            "GROUND_ELEVATION": float(ground_elevation_m) / 1000,
            "SOLAR_AZIMUTH": mul_imd_data.get("SOLAR_AZIMUTH"),
            "SOLAR_ZENITH": mul_imd_data.get("SOLAR_ZENITH"),
            "LOS_AZIMUTH": mul_imd_data.get("LOS_AZIMUTH"),
            "LOS_ZENITH": mul_imd_data.get("LOS_ZENITH"),
            "OUTPUT_RASTER_URI": _wsl_to_windows(flaash_dat),
        }
        if water_vapor_preset is not None:
            params["WATER_VAPOR_PRESET"] = water_vapor_preset

        params = {k: v for k, v in params.items() if v is not None}
        print(f"Running {item['name']}")
        run_flaash(params, flaash_params_txt, envi_engine, flaash_dat)

        if not os.path.exists(flaash_dat):
            rows.append(
                {
                    "scenario": item["name"],
                    "status": "flaash_failed",
                    "water_vapor_preset": water_vapor_preset,
                    "modtran_atm": item["MODTRAN_ATM"],
                    "modtran_aer": item["MODTRAN_AER"],
                    "use_aerosol": item["USE_AEROSOL"],
                    "default_visibility": item["DEFAULT_VISIBILITY"],
                    "dem_ground_percentile": args.dem_ground_percentile,
                    "ground_elevation_km": float(ground_elevation_m) / 1000,
                    "neg_band_1": None,
                    "neg_band_2": None,
                    "neg_band_3": None,
                    "neg_band_4": None,
                    "neg_band_5": None,
                    "neg_band_6": None,
                    "neg_band_7": None,
                    "neg_band_8": None,
                    "neg_avg": None,
                    "flaash_output": flaash_dat,
                }
            )
            continue

        neg = _sample_negative_percent_by_band(flaash_dat)
        rows.append(
            {
                "scenario": item["name"],
                "status": "ok",
                "water_vapor_preset": water_vapor_preset,
                "modtran_atm": item["MODTRAN_ATM"],
                "modtran_aer": item["MODTRAN_AER"],
                "use_aerosol": item["USE_AEROSOL"],
                "default_visibility": item["DEFAULT_VISIBILITY"],
                "dem_ground_percentile": args.dem_ground_percentile,
                "ground_elevation_km": float(ground_elevation_m) / 1000,
                "neg_band_1": neg[0],
                "neg_band_2": neg[1],
                "neg_band_3": neg[2],
                "neg_band_4": neg[3],
                "neg_band_5": neg[4],
                "neg_band_6": neg[5],
                "neg_band_7": neg[6],
                "neg_band_8": neg[7],
                "neg_avg": float(np.mean(neg)),
                "flaash_output": flaash_dat,
            }
        )

    os.makedirs(os.path.dirname(args.output_csv) or ".", exist_ok=True)
    fieldnames = [
        "scenario",
        "status",
        "water_vapor_preset",
        "modtran_atm",
        "modtran_aer",
        "use_aerosol",
        "default_visibility",
        "dem_ground_percentile",
        "ground_elevation_km",
        "neg_band_1",
        "neg_band_2",
        "neg_band_3",
        "neg_band_4",
        "neg_band_5",
        "neg_band_6",
        "neg_band_7",
        "neg_band_8",
        "neg_avg",
        "flaash_output",
    ]
    with open(args.output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Done. Results CSV: {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
