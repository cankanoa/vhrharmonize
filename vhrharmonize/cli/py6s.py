#!/usr/bin/env python3
"""Run Py6S atmospheric correction only (no ortho/pansharpen)."""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime
from typing import Dict, List, Optional


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


def _write_json(path: str, payload: Dict) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


def run_py6s_only(args: argparse.Namespace) -> int:
    from vhrharmonize.preprocess.atmospheric_correction import (
        run_py6s,
    )
    from vhrharmonize.preprocess.orthorectification import (
        gcp_refined_rpc_orthorectification,
        resolve_output_resolution_for_crs,
    )
    from vhrharmonize.providers.worldview.files import find_files
    from vhrharmonize.providers.worldview.metadata import load_worldview_metadata
    from vhrharmonize.providers.standardized import StandardizedMetadata

    filter_basenames = _parse_filter_basenames(args.filter_basename)
    os.makedirs(args.output_dir, exist_ok=True)

    for input_folder in args.input_dir:
        print(f"Scanning input folder: {input_folder}")
        found_default_files = find_files(input_folder, filter_basenames)
        for _, found_default_file in found_default_files.items():
                root_folder_path = found_default_file.get("root_folder_path")
                mul_photo_basename = found_default_file.get("mul_photo_basename")
                print(f"Processing: {mul_photo_basename}")
                required = ("mul_imd_file", "mul_tif_file", "mul_shp_file")
                if not all(found_default_file.get(k) for k in required):
                    print("Skipping scene: required files missing for Py6S-only workflow.")
                    continue

                mul_imd_file = found_default_file["mul_imd_file"]
                mul_tif_file = found_default_file["mul_tif_file"]
                mul_shp_path = found_default_file["mul_shp_file"]
                mul_worldview_metadata = load_worldview_metadata(
                    mul_imd_file,
                    photo_basename=mul_photo_basename,
                )
                mul_metadata = StandardizedMetadata.from_worldview_metadata(mul_worldview_metadata)

                py6s_output_path = os.path.join(args.output_dir, f"{mul_photo_basename}{args.output_suffix}.tif")
                scene_bbox = _scene_bbox_wgs84_from_shp(mul_shp_path)
                py6s_result = run_py6s(
                    input_raster=mul_tif_file,
                    output_raster=py6s_output_path,
                    metadata=mul_metadata,
                    ground_elevation_km=0.0,
                    atmosphere_profile=args.py6s_atmosphere_profile,
                    aerosol_profile=args.py6s_aerosol_profile,
                    aot550=args.py6s_aot550,
                    visibility_km=args.py6s_visibility,
                    water_vapor=args.py6s_water_vapor,
                    ozone=args.py6s_ozone,
                    sixs_executable=args.py6s_executable,
                    output_scale_factor=args.py6s_output_scale_factor,
                    output_dtype=args.py6s_output_dtype,
                    use_imd_radiance_calibration=args.py6s_use_imd_radiance_calibration,
                    use_worldview_gain_offset_adjustment=args.py6s_use_worldview_gain_offset_adjustment,
                    auto_atmos_source=args.py6s_auto_atmos_source,
                    bbox_wgs84=scene_bbox,
                    auto_atmos_grid_size=args.py6s_auto_atmos_grid_size,
                    auto_atmos_search_days=args.py6s_auto_atmos_search_days,
                    auto_atmos_timeout_s=args.py6s_auto_atmos_timeout_s,
                    auto_atmos_power_endpoint=args.py6s_auto_atmos_power_endpoint,
                )
                final_output_path = py6s_output_path
                if args.epsg is not None:
                    if not args.dem_file_path:
                        raise ValueError(
                            "--dem-file-path is required when --epsg is set for vhr-py6s output orthorectification."
                        )
                    ortho_output_path = os.path.join(
                        args.output_dir,
                        f"{mul_photo_basename}{args.output_suffix}_ortho.tif",
                    )
                    print(f"Orthorectifying Py6S output to EPSG:{args.epsg}")
                    gcp_refined_rpc_orthorectification(
                        py6s_output_path,
                        ortho_output_path,
                        args.dem_file_path,
                        args.epsg,
                        output_nodata_value=None,
                        dtype=args.py6s_output_dtype,
                        output_resolution=resolve_output_resolution_for_crs(
                            args.epsg,
                            mul_metadata.product_resolution,
                        ),
                    )
                    if args.keep_intermediate_py6s:
                        final_output_path = ortho_output_path
                    else:
                        os.replace(ortho_output_path, py6s_output_path)
                        final_output_path = py6s_output_path
                print(f"Wrote: {final_output_path}")

                report_path = os.path.join(
                    args.output_dir,
                    f"{mul_photo_basename}{args.output_suffix}_metadata.json",
                )
                report = {
                    "scene": {
                        "input_dir": input_folder,
                        "scene_root": root_folder_path,
                        "mul_photo_basename": mul_photo_basename,
                        "completed_utc": datetime.utcnow().isoformat() + "Z",
                    },
                    "inputs": {
                        "mul_imd_file": mul_imd_file,
                        "mul_tif_file": mul_tif_file,
                        "mul_shp_file": mul_shp_path,
                    },
                    "standardized_metadata": mul_metadata.to_dict(),
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
                        "effective": {
                            **py6s_result.effective_params,
                        },
                        "auto_atmos_estimate": py6s_result.auto_atmos_estimate,
                    },
                    "orthorectification": {
                        "enabled": args.epsg is not None,
                        "output_epsg": args.epsg,
                        "dem_file_path": args.dem_file_path,
                        "keep_intermediate_py6s": args.keep_intermediate_py6s,
                    },
                    "outputs": {
                        "py6s_output_path": py6s_output_path,
                        "final_output_path": final_output_path,
                        "metadata_report_path": report_path,
                    },
                }
                _write_json(report_path, report)
                print(f"Wrote metadata: {report_path}")

    print("All processing complete")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run Py6S-only atmospheric correction on WorldView scenes.")
    parser.add_argument("--config-yaml", help="Optional YAML config file.")
    parser.add_argument("--input-dir", action="append", help="Input directory to scan for scenes.")
    parser.add_argument(
        "--filter-basename",
        action="append",
        help="Optional image basenames to process. Repeat or pass comma-separated values.",
    )
    parser.add_argument("--output-dir", required=False, default="Scene_Output", help="Output directory.")
    parser.add_argument("--output-suffix", default="_py6s", help="Suffix for Py6S output filenames.")
    parser.add_argument(
        "--epsg",
        type=int,
        help=(
            "Optional output EPSG for post-Py6S orthorectification. "
            "If set, --dem-file-path is required."
        ),
    )
    parser.add_argument(
        "--dem-file-path",
        help="DEM path used for RPC orthorectification when --epsg is set.",
    )
    parser.add_argument(
        "--keep-intermediate-py6s",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "When orthorectifying in vhr-py6s, keep both outputs: raw Py6S and *_ortho.tif. "
            "Default false replaces raw file with orthorectified output."
        ),
    )
    parser.add_argument("--py6s-atmosphere-profile", default="user")
    parser.add_argument("--py6s-aerosol-profile", default="maritime")
    parser.add_argument("--py6s-aot550", type=float, default=0.2)
    parser.add_argument("--py6s-visibility", type=float)
    parser.add_argument("--py6s-water-vapor", type=float, default=2.5)
    parser.add_argument("--py6s-ozone", type=float, default=0.3)
    parser.add_argument("--py6s-output-scale-factor", type=float, default=10000.0)
    parser.add_argument("--py6s-output-dtype", default="int16")
    parser.add_argument("--py6s-executable")
    parser.add_argument(
        "--py6s-use-imd-radiance-calibration",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--py6s-use-worldview-gain-offset-adjustment",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--py6s-auto-atmos-source",
        choices=["none", "nasa_power"],
        default="none",
    )
    parser.add_argument("--py6s-auto-atmos-grid-size", type=int, default=3)
    parser.add_argument("--py6s-auto-atmos-search-days", type=int, default=1)
    parser.add_argument("--py6s-auto-atmos-timeout-s", type=float, default=30.0)
    parser.add_argument(
        "--py6s-auto-atmos-power-endpoint",
        default="https://power.larc.nasa.gov/api/temporal/daily/point",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument("--config-yaml")
    config_args, _ = config_parser.parse_known_args(argv)

    config_defaults: Dict = {}
    if config_args.config_yaml:
        config_defaults = _normalize_config_defaults(_load_yaml_config(config_args.config_yaml))

    parser = build_parser()
    if config_defaults:
        parser.set_defaults(**config_defaults)
    args = parser.parse_args(argv)

    if not args.input_dir:
        parser.error("--input-dir is required (via CLI or --config-yaml).")
    if args.py6s_auto_atmos_grid_size < 1:
        parser.error("--py6s-auto-atmos-grid-size must be >= 1.")
    if args.py6s_auto_atmos_search_days < 0:
        parser.error("--py6s-auto-atmos-search-days must be >= 0.")
    if args.py6s_auto_atmos_timeout_s <= 0:
        parser.error("--py6s-auto-atmos-timeout-s must be > 0.")
    if args.epsg is not None and not args.dem_file_path:
        parser.error("--dem-file-path is required when --epsg is set.")
    return run_py6s_only(args)


if __name__ == "__main__":
    sys.exit(main())
