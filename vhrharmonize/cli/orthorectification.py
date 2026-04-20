"""CLI wrappers for orthorectification preprocess helpers."""

from __future__ import annotations

import argparse
import os
import sys

from vhrharmonize.preprocess.orthorectification import (
    gcp_refined_rpc_orthorectification,
    qgis_gcps_to_csv,
    qgis_gcps_to_geojson,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run orthorectification helper commands.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    orthorectify = subparsers.add_parser("orthorectify", help="Run RPC orthorectification.")
    orthorectify.add_argument("--input-image-path", required=True)
    orthorectify.add_argument("--output-image-path", required=True)
    orthorectify.add_argument("--dem-image-path", required=True)
    orthorectify.add_argument("--output-epsg", required=True, type=int)
    orthorectify.add_argument("--gcp-geojson-file-path")
    orthorectify.add_argument("--output-nodata-value", type=float)
    orthorectify.add_argument("--dtype")
    orthorectify.add_argument("--output-resolution", type=float)
    orthorectify.add_argument("--output-resolution-x", type=float)
    orthorectify.add_argument("--output-resolution-y", type=float)

    csv_parser = subparsers.add_parser("qgis-gcps-to-csv", help="Convert QGIS GCP text to CSV.")
    csv_parser.add_argument("--input-gcp-path", required=True)
    csv_parser.add_argument("--output-csv-path", required=True)
    csv_parser.add_argument("--output-epsg", type=int)

    geojson_parser = subparsers.add_parser(
        "qgis-gcps-to-geojson",
        help="Convert QGIS GCP text to Orthority-compatible GeoJSON.",
    )
    geojson_parser.add_argument("--input-image-path", required=True)
    geojson_parser.add_argument("--qgis-gcp-file-path", required=True)
    geojson_parser.add_argument("--file-name", required=True)
    geojson_parser.add_argument("--dem-file-path", required=True)
    geojson_parser.add_argument("--output-geojson-path", required=True)
    geojson_parser.add_argument("--force-positive-pixel-values", action="store_true")
    return parser


def _resolve_output_resolution(args):
    if args.output_resolution is not None:
        return float(args.output_resolution)
    if args.output_resolution_x is not None or args.output_resolution_y is not None:
        if args.output_resolution_x is None or args.output_resolution_y is None:
            raise ValueError(
                "Provide both --output-resolution-x and --output-resolution-y when using tuple resolution."
            )
        return (float(args.output_resolution_x), float(args.output_resolution_y))
    return None


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.command == "orthorectify":
        gcp_refined_rpc_orthorectification(
            input_image_path=args.input_image_path,
            output_image_path=args.output_image_path,
            dem_image_path=args.dem_image_path,
            output_epsg=args.output_epsg,
            gcp_geojson_file_path=args.gcp_geojson_file_path,
            output_nodata_value=args.output_nodata_value,
            dtype=args.dtype,
            output_resolution=_resolve_output_resolution(args),
        )
        return 0

    if args.command == "qgis-gcps-to-csv":
        csv_text = qgis_gcps_to_csv(
            input_gcp_path=args.input_gcp_path,
            output_epsg=args.output_epsg,
        )
        os.makedirs(os.path.dirname(args.output_csv_path) or ".", exist_ok=True)
        with open(args.output_csv_path, "w", encoding="utf-8") as handle:
            handle.write(csv_text)
        return 0

    if args.command == "qgis-gcps-to-geojson":
        qgis_gcps_to_geojson(
            input_image_path=args.input_image_path,
            qgis_gcp_file_path=args.qgis_gcp_file_path,
            file_name=args.file_name,
            dem_file_path=args.dem_file_path,
            output_geojson_path=args.output_geojson_path,
            force_positive_pixel_values=args.force_positive_pixel_values,
        )
        return 0

    raise ValueError(f"Unsupported orthorectification command: {args.command}")


__all__ = ["main"]


if __name__ == "__main__":
    sys.exit(main())
