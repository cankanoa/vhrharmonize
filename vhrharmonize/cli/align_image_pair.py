#!/usr/bin/env python3
"""Align one image (moving) to another image (fixed) using elastix."""

import argparse
import json
import os
import sys
from typing import Optional

from vhrharmonize.pipelines.alignment import align_image_pair


def build_parser() -> argparse.ArgumentParser:
    """Build CLI parser for pairwise image alignment."""
    parser = argparse.ArgumentParser(
        description=(
            "Align a moving image (A) to a fixed image (B). "
            "Default behavior is structural WV/LiDAR alignment on full extent."
        ),
    )
    parser.add_argument("--moving-image", required=True, help="Path to moving image A (will be warped).")
    parser.add_argument("--fixed-image", required=True, help="Path to fixed/reference image B.")
    parser.add_argument("--output-image", required=True, help="Path to output aligned image.")
    parser.add_argument(
        "--band-index",
        type=int,
        default=0,
        help="0-based band index used for registration metric (default: 0).",
    )
    parser.add_argument(
        "--moving-band-index",
        type=int,
        help="Optional 0-based moving-image band index for registration metric.",
    )
    parser.add_argument(
        "--fixed-band-index",
        type=int,
        help="Optional 0-based fixed-image band index for registration metric.",
    )
    parser.add_argument(
        "--no-tiling",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Disable/enable tiling (default: no tiling, full-extent registration).",
    )
    parser.add_argument(
        "--tile-size",
        type=int,
        default=1000,
        help="Tile size in pixels when tiling is enabled (default: 1000).",
    )
    parser.add_argument(
        "--tile-buffer",
        type=int,
        default=100,
        help="Buffer/overlap in pixels around each tile for registration stability (default: 100).",
    )
    parser.add_argument(
        "--parameter-map",
        default="rigid",
        help="Elastix default parameter map (e.g., rigid, affine, bspline).",
    )
    parser.add_argument(
        "--parameter-file",
        action="append",
        help="Optional elastix parameter file path. Repeat for multi-stage registration.",
    )
    parser.add_argument(
        "--moving-nodata",
        type=float,
        help="Optional override nodata value for moving image masking.",
    )
    parser.add_argument(
        "--fixed-nodata",
        type=float,
        help="Optional override nodata value for fixed image masking.",
    )
    parser.add_argument(
        "--output-nodata",
        type=float,
        help="Optional output nodata value. Defaults to moving nodata, then fixed nodata, else 0.",
    )
    parser.add_argument(
        "--min-valid-fraction",
        type=float,
        default=0.01,
        help="Minimum valid-mask fraction per tile required to run elastix (default: 0.01).",
    )
    parser.add_argument("--temp-dir", help="Optional parent directory for temporary tile files.")
    parser.add_argument(
        "--keep-temp-dir",
        action="store_true",
        help="Keep temporary tile directory for debugging.",
    )
    parser.add_argument(
        "--log-to-console",
        action="store_true",
        help="Enable verbose elastix logging.",
    )
    parser.add_argument(
        "--clip-fixed-to-moving",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Clip fixed image domain to moving-image bounds before alignment (default: enabled).",
    )
    parser.add_argument(
        "--output-on-moving-grid",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Write aligned output on the moving-image grid (default: true). "
            "Disable with --no-output-on-moving-grid to write on fixed-image grid."
        ),
    )
    parser.add_argument(
        "--enforce-mutual-valid-mask",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Use only pixels valid in both fixed and moving images for both elastix masks "
            "(default: true)."
        ),
    )
    parser.add_argument(
        "--registration-mode",
        choices=["default", "structural_wv3_lidar"],
        default="structural_wv3_lidar",
        help=(
            "Registration strategy. `default` uses raw bands; "
            "`structural_wv3_lidar` mirrors the elastix-wrapper flow "
            "(common-grid ROI + edge proxies + translation->rigid chain)."
        ),
    )
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    """CLI entrypoint for tile-aware pairwise image alignment."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if not os.path.isfile(args.moving_image):
        parser.error(f"--moving-image does not exist: {args.moving_image}")
    if not os.path.isfile(args.fixed_image):
        parser.error(f"--fixed-image does not exist: {args.fixed_image}")
    if args.band_index < 0:
        parser.error("--band-index must be >= 0.")
    if args.moving_band_index is not None and args.moving_band_index < 0:
        parser.error("--moving-band-index must be >= 0.")
    if args.fixed_band_index is not None and args.fixed_band_index < 0:
        parser.error("--fixed-band-index must be >= 0.")
    if args.tile_size <= 0:
        parser.error("--tile-size must be > 0.")
    if args.tile_buffer < 0:
        parser.error("--tile-buffer must be >= 0.")
    if args.min_valid_fraction <= 0 or args.min_valid_fraction > 1:
        parser.error("--min-valid-fraction must be in (0, 1].")

    result = align_image_pair(
        moving_image_path=args.moving_image,
        fixed_image_path=args.fixed_image,
        output_image_path=args.output_image,
        band_index=args.band_index,
        moving_band_index=args.moving_band_index,
        fixed_band_index=args.fixed_band_index,
        tiling=not args.no_tiling,
        tile_size=args.tile_size,
        tile_buffer=args.tile_buffer,
        parameter_map=args.parameter_map,
        parameter_file_paths=args.parameter_file,
        moving_nodata=args.moving_nodata,
        fixed_nodata=args.fixed_nodata,
        output_nodata=args.output_nodata,
        min_valid_fraction=args.min_valid_fraction,
        temp_dir=args.temp_dir,
        keep_temp_dir=args.keep_temp_dir,
        log_to_console=args.log_to_console,
        clip_fixed_to_moving=args.clip_fixed_to_moving,
        output_on_moving_grid=args.output_on_moving_grid,
        enforce_mutual_valid_mask=args.enforce_mutual_valid_mask,
        registration_mode=args.registration_mode,
    )

    print(
        json.dumps(
            {
                "output_image_path": result.output_image_path,
                "total_tiles": result.total_tiles,
                "successful_tiles": result.successful_tiles,
                "skipped_tiles": result.skipped_tiles,
                "temp_dir": result.temp_dir,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
