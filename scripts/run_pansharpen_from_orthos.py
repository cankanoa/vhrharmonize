#!/usr/bin/env python3
"""Standalone pansharpen from existing orthorectified MS + PAN rasters."""

import argparse
import importlib.util
import os
import sys

def _load_pansharpen_func():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    module_path = os.path.join(os.path.dirname(script_dir), "vhrharmonize", "pansharpen.py")
    spec = importlib.util.spec_from_file_location("vhr_pansharpen", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load pansharpen module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.pansharpen_image


def main(argv=None) -> int:
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

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    pansharpen_image = _load_pansharpen_func()
    pansharpen_image(
        args.mul_ortho,
        args.pan_ortho,
        args.output,
        change_nodata_value=args.nodata_value,
    )
    print(f"Wrote pansharpened raster: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
