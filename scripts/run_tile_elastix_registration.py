#!/usr/bin/env python3
"""Register output imagery tiles to LiDAR intensity tiles using elastix-wrapper."""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional

import rasterio


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
    return {k.replace("-", "_"): v for k, v in loaded.items()}


def _normalize_config_defaults(config_defaults: Dict) -> Dict:
    normalized = dict(config_defaults)
    for list_key in ("filter_tile_id",):
        if list_key in normalized and isinstance(normalized[list_key], str):
            normalized[list_key] = [normalized[list_key]]
    return normalized


def _sanitize_for_filename(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in str(value))


def _parse_filter_tile_ids(raw_values: Optional[List[str]]) -> List[str]:
    if not raw_values:
        return []
    parsed: List[str] = []
    for raw in raw_values:
        for value in raw.split(","):
            value = value.strip()
            if value:
                parsed.append(value)
    return parsed


def _discover_tiles(directory: str, recursive: bool) -> List[str]:
    pattern = os.path.join(directory, "**", "*.tif") if recursive else os.path.join(directory, "*.tif")
    return sorted(glob.glob(pattern, recursive=recursive))


def _strip_suffix(name: str, suffix: str) -> Optional[str]:
    if name.endswith(suffix):
        return name[: -len(suffix)]
    return None


def _build_intensity_tile_map(
    intensity_files: List[str],
    *,
    intensity_dimension: str,
    intensity_suffix: str,
) -> Dict[str, str]:
    result = {}
    expected_tail = f"_{intensity_dimension}{intensity_suffix}.tif"
    for path in intensity_files:
        basename = os.path.basename(path)
        tile_id = _strip_suffix(basename, expected_tail)
        if tile_id is None:
            continue
        result[tile_id] = path
    return result


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
    ulx = bounds.left
    uly = bounds.top
    lrx = bounds.right
    lry = bounds.bottom
    xres = abs(transform.a)
    if xres <= 0:
        raise ValueError(f"Invalid pixel size for fixed raster: {fixed_path}")
    return epsg, ulx, uly, lrx, lry, xres


def _extract_single_band(input_path: str, output_path: str, band_index: int) -> None:
    cmd = [
        "gdal_translate",
        "-q",
        "-b",
        str(band_index),
        input_path,
        output_path,
    ]
    subprocess.run(cmd, check=True)


def _run_pair(
    *,
    elastix_script: str,
    fixed_path: str,
    moving_path: str,
    outdir: str,
    moving_band_index: int,
) -> None:
    epsg, ulx, uly, lrx, lry, res = _get_fixed_grid_params(fixed_path)

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


def run(args: argparse.Namespace) -> int:
    if not os.path.isfile(args.elastix_script):
        raise FileNotFoundError(f"Elastix wrapper script not found: {args.elastix_script}")
    if not os.access(args.elastix_script, os.X_OK):
        raise PermissionError(f"Elastix wrapper script is not executable: {args.elastix_script}")

    if "ELASTIX_EXE" not in os.environ or "TRANSFORMIX_EXE" not in os.environ:
        raise EnvironmentError("ELASTIX_EXE and TRANSFORMIX_EXE must be set in environment.")

    filter_tile_ids = set(_parse_filter_tile_ids(args.filter_tile_id))
    work_root = args.output_dir or "/tmp/elastix_tile_work"
    os.makedirs(work_root, exist_ok=True)

    intensity_files = _discover_tiles(args.intensity_dir, args.recursive)
    moving_files = _discover_tiles(args.moving_dir, args.recursive)
    intensity_map = _build_intensity_tile_map(
        intensity_files,
        intensity_dimension=args.intensity_dimension,
        intensity_suffix=args.intensity_suffix,
    )

    if not intensity_map:
        print("No intensity tiles matched expected naming pattern.")
        return 0

    total_pairs = 0
    success_pairs = 0
    for tile_id, fixed_path in sorted(intensity_map.items()):
        if filter_tile_ids and tile_id not in filter_tile_ids:
            continue
        moving_paths = [
            p for p in moving_files
            if os.path.basename(p).startswith(f"{tile_id}_")
            and os.path.basename(p).endswith(f"{args.moving_suffix}.tif")
        ]
        if not moving_paths:
            print(f"Tile {tile_id}: no moving imagery tiles found, skipping.")
            continue

        for moving_path in sorted(moving_paths):
            moving_base = os.path.splitext(os.path.basename(moving_path))[0]
            aligned_path = os.path.join(
                os.path.dirname(moving_path),
                f"{moving_base}{args.aligned_suffix}.tif",
            )
            outdir = os.path.join(
                work_root,
                _sanitize_for_filename(tile_id),
                _sanitize_for_filename(moving_base),
            )
            final_georef = os.path.join(outdir, "warp", "result_georef.tif")
            total_pairs += 1

            if not args.overwrite and os.path.exists(aligned_path):
                print(f"Tile {tile_id}: aligned output exists, skipping {moving_base}.")
                success_pairs += 1
                continue

            os.makedirs(outdir, exist_ok=True)
            print(f"Registering tile {tile_id}: {os.path.basename(moving_path)} -> {os.path.basename(fixed_path)}")
            try:
                _run_pair(
                    elastix_script=args.elastix_script,
                    fixed_path=fixed_path,
                    moving_path=moving_path,
                    outdir=outdir,
                    moving_band_index=args.moving_band_index,
                )
                if not os.path.exists(final_georef):
                    raise FileNotFoundError(f"Elastix output missing: {final_georef}")
                shutil.copy2(final_georef, aligned_path)
                print(f"Tile {tile_id}: wrote aligned output {aligned_path}")
                success_pairs += 1
            except Exception as exc:
                print(f"Tile {tile_id}: registration failed for {moving_base}: {exc}")

    print(f"Done. Pair attempts: {total_pairs}; successful/available outputs: {success_pairs}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run elastix-wrapper registration for existing imagery tiles against intensity tiles.",
    )
    parser.add_argument("--config-yaml", help="Optional YAML config file.")
    parser.add_argument("--moving-dir", help="Directory containing moving imagery tiles (e.g., cloudmasked/final tiles).")
    parser.add_argument("--intensity-dir", help="Directory containing fixed intensity mean tiles.")
    parser.add_argument(
        "--output-dir",
        help=(
            "Optional work directory for per-tile elastix intermediates. "
            "Final aligned outputs are written next to each moving tile."
        ),
    )
    parser.add_argument(
        "--elastix-script",
        default="/home/milo/repos/elastix-wrapper/run_elastix_wv3_cbl_to_intensity.sh",
        help="Path to elastix-wrapper shell script.",
    )
    parser.add_argument(
        "--moving-suffix",
        default="_final",
        help="Suffix used by moving imagery tile filenames (before .tif).",
    )
    parser.add_argument(
        "--intensity-dimension",
        default="Intensity",
        help="Dimension token used in intensity tile names.",
    )
    parser.add_argument(
        "--intensity-suffix",
        default="_mean",
        help="Suffix used by intensity tile filenames (before .tif).",
    )
    parser.add_argument(
        "--moving-band-index",
        type=int,
        default=1,
        help="1-based band index extracted from moving tiles for registration (WV3 coastal blue is band 1).",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Recursively search moving and intensity directories for tiles.",
    )
    parser.add_argument(
        "--filter-tile-id",
        action="append",
        help="Optional tile id filter. Repeat or pass comma-separated list.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing registration outputs.",
    )
    parser.add_argument(
        "--aligned-suffix",
        default="_aligned",
        help="Suffix appended to moving tile basename for final aligned output.",
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

    if not args.moving_dir:
        parser.error("--moving-dir is required (or set in --config-yaml).")
    if not args.intensity_dir:
        parser.error("--intensity-dir is required (or set in --config-yaml).")
    if args.moving_band_index < 1:
        parser.error("--moving-band-index must be >= 1.")

    return run(args)


if __name__ == "__main__":
    sys.exit(main())
