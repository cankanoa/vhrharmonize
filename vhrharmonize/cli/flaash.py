"""Run shared FLAASH preprocessing directly from the CLI."""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Any

from vhrharmonize.preprocess.atmospheric_correction import (
    run_flaash,
)


def _load_params(args: argparse.Namespace) -> dict[str, Any]:
    """Load FLAASH parameters from CLI input.
    Args:
        args: Parsed CLI arguments.
    Returns:
        Deserialized FLAASH parameter dictionary.
    """
    if bool(args.params_json) == bool(args.params_json_file):
        raise ValueError("Provide exactly one of --params-json or --params-json-file.")

    if args.params_json_file:
        with open(args.params_json_file, "r", encoding="utf-8") as handle:
            params = json.load(handle)
    else:
        params = json.loads(args.params_json)

    if not isinstance(params, dict):
        raise ValueError("FLAASH params must deserialize to a JSON object.")
    return params


def _build_parser() -> argparse.ArgumentParser:
    """Build the FLAASH CLI parser.
    Args:
        None.
    Returns:
        Configured argument parser.
    """
    parser = argparse.ArgumentParser(description="Run ENVI FLAASH with a supplied parameter object.")
    parser.add_argument("--params-json", help="Inline JSON object of FLAASH parameters.")
    parser.add_argument("--params-json-file", help="Path to a JSON file containing FLAASH parameters.")
    parser.add_argument("--output-params-path", required=True, help="Path to write the executed params dump.")
    parser.add_argument("--envi-engine-path", required=True, help="ENVI engine path used by envipyengine.")
    parser.add_argument(
        "--convert-paths-for-windows",
        action="store_true",
        help="Convert WSL /mnt/<drive>/ paths into Windows paths before execution.",
    )
    parser.add_argument(
        "--delete-output-before-run",
        help="Optional output image path to delete before running FLAASH.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the FLAASH CLI.
    Args:
        argv: Optional command line arguments.
    Returns:
        Process exit code.
    """
    args = _build_parser().parse_args(argv)
    params = _load_params(args)
    output_raster = params.get("OUTPUT_RASTER_URI")
    input_raster = None
    if isinstance(params.get("INPUT_RASTER"), dict):
        input_raster = params["INPUT_RASTER"].get("url")
    if not input_raster or not output_raster:
        raise ValueError("FLAASH params must include INPUT_RASTER.url and OUTPUT_RASTER_URI.")

    os.makedirs(os.path.dirname(args.output_params_path) or ".", exist_ok=True)
    run_flaash(
        input_raster=input_raster,
        output_raster=output_raster,
        params=params,
        output_params_path=args.output_params_path,
        envi_engine_path=args.envi_engine_path,
        convert_paths_for_windows=args.convert_paths_for_windows,
        delete_output_before_run=args.delete_output_before_run,
    )
    return 0


__all__ = ["main"]


if __name__ == "__main__":
    sys.exit(main())
