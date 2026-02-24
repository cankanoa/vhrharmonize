#!/usr/bin/env python3
"""Run OmniCloudMask on an existing raster and write masked output."""

import argparse
import json
import os
import sys
from typing import Dict, Optional


def _parse_json_dict(raw_json: Optional[str]) -> Dict:
    if not raw_json:
        return {}
    parsed = json.loads(raw_json)
    if not isinstance(parsed, dict):
        raise ValueError("Expected a JSON object for omnicloud kwargs.")
    return parsed


def _build_output_path(input_path: str, suffix: str) -> str:
    base, ext = os.path.splitext(input_path)
    return f"{base}{suffix}{ext}"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Apply OmniCloudMask to an existing raster (standalone).",
    )
    parser.add_argument("--input-raster", required=True, help="Input raster path (e.g. pansharpened tif).")
    parser.add_argument(
        "--output-raster",
        help="Output masked raster path. Default: <input>_cloudmasked.tif",
    )
    parser.add_argument(
        "--output-mask",
        help="Output binary mask path. Default: <input>_cloudmask.tif",
    )
    parser.add_argument("--nodata-value", type=float, default=-9999, help="NoData for masked output.")
    parser.add_argument("--red-band-index", type=int, default=5, help="1-based red band index.")
    parser.add_argument("--green-band-index", type=int, default=3, help="1-based green band index.")
    parser.add_argument("--nir-band-index", type=int, default=7, help="1-based NIR band index.")
    parser.add_argument("--cloud-classes", default="1,2,3", help="Comma-separated cloud/cloud-shadow class ids.")
    parser.add_argument("--buffer-pixels", type=int, default=0, help="Cloud-mask dilation buffer size in pixels.")
    parser.add_argument(
        "--omnicloud-kwargs-json",
        help='Optional JSON kwargs passed to omnicloudmask predict_from_array, e.g. \'{"batch_size": 8}\'',
    )
    args = parser.parse_args(argv)

    if not os.path.isfile(args.input_raster):
        parser.error(f"--input-raster does not exist: {args.input_raster}")
    if args.buffer_pixels < 0:
        parser.error("--buffer-pixels must be >= 0.")

    cloud_classes = [int(x.strip()) for x in args.cloud_classes.split(",") if x.strip()]
    omnicloud_kwargs = _parse_json_dict(args.omnicloud_kwargs_json)

    output_raster = args.output_raster or _build_output_path(args.input_raster, "_cloudmasked")
    output_mask = args.output_mask or _build_output_path(args.input_raster, "_cloudmask")

    from vhrharmonize import create_cloud_mask_with_omnicloudmask, apply_binary_cloud_mask_to_image

    print(f"Running OmniCloudMask on: {args.input_raster}")
    create_cloud_mask_with_omnicloudmask(
        args.input_raster,
        output_mask,
        red_band_index=args.red_band_index,
        green_band_index=args.green_band_index,
        nir_band_index=args.nir_band_index,
        cloud_classes=cloud_classes,
        buffer_pixels=args.buffer_pixels,
        omnicloud_kwargs=omnicloud_kwargs,
    )
    apply_binary_cloud_mask_to_image(
        args.input_raster,
        output_mask,
        output_raster,
        output_nodata_value=args.nodata_value,
    )
    print(f"Wrote mask: {output_mask}")
    print(f"Wrote masked raster: {output_raster}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
