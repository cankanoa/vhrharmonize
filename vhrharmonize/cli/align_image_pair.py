#!/usr/bin/env python3
"""Align one image (moving) to another image (fixed) using coregix."""

import argparse
import json
import os
import sys
from typing import Optional

from vhrharmonize.preprocess.alignment import align_image_pair


def build_parser() -> argparse.ArgumentParser:
    """Build CLI parser for pairwise image alignment."""
    parser = argparse.ArgumentParser(
        description=(
            "Align a moving image (A) to a fixed image (B). "
            "Uses edge-proxy registration on full extent by default."
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
        "--use-edge-proxies",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use edge-proxy images rather than raw intensities for registration (default: true).",
    )
    parser.add_argument(
        "--split-factor",
        type=int,
        default=0,
        help="Split the solve and apply domains into 2^k chunks (default: 0).",
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
        help="Minimum valid-mask fraction required to run alignment (default: 0.01).",
    )
    parser.add_argument(
        "--solve-resolution",
        type=float,
        help="Optional target pixel size, in raster CRS units, for the registration solve.",
    )
    parser.add_argument("--temp-dir", help="Optional parent directory for temporary working files.")
    parser.add_argument(
        "--keep-temp-dir",
        action="store_true",
        help="Keep the temporary working directory for debugging.",
    )
    parser.add_argument(
        "--clip-fixed-to-moving",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Clip fixed image domain to moving-image bounds before alignment.",
    )
    parser.add_argument(
        "--output-on-moving-grid",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write aligned output on the moving-image grid (default: true).",
    )
    parser.add_argument(
        "--trim-edge-invalid",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="After alignment, set pixels adjacent to irregular exterior invalid boundaries to nodata.",
    )
    parser.add_argument(
        "--edge-trim-depth",
        type=int,
        default=8,
        help="Number of pixels to trim inward from each exterior invalid boundary (default: 8).",
    )
    parser.add_argument(
        "--edge-trim-detection-band-index",
        type=int,
        default=0,
        help="0-based band index used to detect edge artifacts (default: 0).",
    )
    parser.add_argument(
        "--edge-trim-invalid-below",
        type=float,
        help="For edge trimming, treat values <= this threshold as invalid.",
    )
    parser.add_argument(
        "--edge-trim-invalid-above",
        type=float,
        help="For edge trimming, treat values >= this threshold as invalid.",
    )
    parser.add_argument(
        "--enforce-mutual-valid-mask",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Use only pixels valid in both fixed and moving images for both masks.",
    )
    parser.add_argument(
        "--log-to-console",
        action="store_true",
        help="Enable alignment logging.",
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
    if args.min_valid_fraction <= 0 or args.min_valid_fraction > 1:
        parser.error("--min-valid-fraction must be in (0, 1].")
    if args.solve_resolution is not None and args.solve_resolution <= 0:
        parser.error("--solve-resolution must be > 0.")
    if args.split_factor < 0:
        parser.error("--split-factor must be >= 0.")
    if args.edge_trim_depth <= 0:
        parser.error("--edge-trim-depth must be > 0.")
    if args.edge_trim_detection_band_index < 0:
        parser.error("--edge-trim-detection-band-index must be >= 0.")

    result = align_image_pair(
        moving_image_path=args.moving_image,
        fixed_image_path=args.fixed_image,
        output_image_path=args.output_image,
        band_index=args.band_index,
        moving_band_index=args.moving_band_index,
        fixed_band_index=args.fixed_band_index,
        moving_nodata=args.moving_nodata,
        fixed_nodata=args.fixed_nodata,
        output_nodata=args.output_nodata,
        min_valid_fraction=args.min_valid_fraction,
        temp_dir=args.temp_dir,
        keep_temp_dir=args.keep_temp_dir,
        split_factor=args.split_factor,
        clip_fixed_to_moving=args.clip_fixed_to_moving,
        output_on_moving_grid=args.output_on_moving_grid,
        trim_edge_invalid=args.trim_edge_invalid,
        edge_trim_depth=args.edge_trim_depth,
        edge_trim_detection_band_index=args.edge_trim_detection_band_index,
        edge_trim_invalid_below=args.edge_trim_invalid_below,
        edge_trim_invalid_above=args.edge_trim_invalid_above,
        enforce_mutual_valid_mask=args.enforce_mutual_valid_mask,
        use_edge_proxies=args.use_edge_proxies,
        solve_resolution=args.solve_resolution,
        log_to_console=args.log_to_console,
    )

    print(
        json.dumps(
            {
                "output_image_path": result.output_image_path,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
