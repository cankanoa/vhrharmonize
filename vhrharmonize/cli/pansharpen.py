#!/usr/bin/env python3
"""Standalone pansharpen from existing orthorectified MS + PAN rasters."""

import argparse
import os
import sys
from typing import Optional


def main(argv: Optional[list[str]] = None) -> int:
    """Run the standalone pansharpen CLI.
    Args:
        argv: Optional command line arguments.
    Returns:
        Process exit code.
    """
    parser = argparse.ArgumentParser(
        description="Run pansharpening from existing orthorectified MS/PAN rasters.",
    )
    parser.add_argument("--mul-ortho", required=True, help="Orthorectified multispectral raster path.")
    parser.add_argument("--pan-ortho", required=True, help="Orthorectified panchromatic raster path.")
    parser.add_argument("--output", required=True, help="Output pansharpened raster path.")
    parser.add_argument(
        "--nodata-value",
        type=float,
        default=-9999,
        help="Output NoData value to apply after pansharpening (default: -9999).",
    )
    args = parser.parse_args(argv)

    if not os.path.isfile(args.mul_ortho):
        parser.error(f"--mul-ortho does not exist: {args.mul_ortho}")
    if not os.path.isfile(args.pan_ortho):
        parser.error(f"--pan-ortho does not exist: {args.pan_ortho}")

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    from vhrharmonize.preprocess.pansharpening import pansharpen_image

    pansharpen_image(
        args.mul_ortho,
        args.pan_ortho,
        args.output,
        change_nodata_value=args.nodata_value,
    )
    print(f"Wrote pansharpened raster: {args.output}")
    return 0


__all__ = ["main"]


if __name__ == "__main__":
    sys.exit(main())
