"""CLI for radiometric normalization via SpectralMatch."""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any

from vhrharmonize.preprocess.radiometric_normalization import radiometric_normalization


def _coerce_unknown_arg_value(raw_value: str) -> Any:
    """Coerce an unknown CLI value into a Python scalar.
    Args:
        raw_value: Raw CLI token value.
    Returns:
        Parsed Python value.
    """
    lowered = raw_value.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    if lowered == "null":
        return None
    try:
        if "." in raw_value:
            return float(raw_value)
        return int(raw_value)
    except ValueError:
        return raw_value


def _apply_unknown_match_args(args: argparse.Namespace, unknown_args: list[str]) -> None:
    """Apply passthrough SpectralMatch CLI arguments.
    Args:
        args: Parsed CLI namespace to mutate.
        unknown_args: Unknown CLI tokens to interpret as match_* arguments.
    Returns:
        None.
    """
    idx = 0
    while idx < len(unknown_args):
        token = unknown_args[idx]
        if not token.startswith("--"):
            raise SystemExit(f"Unsupported extra argument syntax: {token}")

        key_token = token[2:]
        inline_value = None
        if "=" in key_token:
            key_token, inline_value = key_token.split("=", 1)
        normalized_key = key_token.replace("-", "_")
        if not normalized_key.startswith("match_"):
            raise SystemExit(f"Unrecognized argument: --{key_token}")

        if inline_value is not None:
            value = _coerce_unknown_arg_value(inline_value)
        elif idx + 1 < len(unknown_args) and not unknown_args[idx + 1].startswith("--"):
            idx += 1
            value = _coerce_unknown_arg_value(unknown_args[idx])
        else:
            value = True

        setattr(args, normalized_key, value)
        idx += 1


def _collect_prefixed_kwargs(namespace: argparse.Namespace, prefix: str) -> dict[str, Any]:
    """Collect namespace values with a shared prefix.
    Args:
        namespace: Parsed CLI namespace.
        prefix: Prefix to strip from matching keys.
    Returns:
        Collected keyword arguments.
    """
    collected: dict[str, Any] = {}
    for key, value in vars(namespace).items():
        if not key.startswith(prefix) or value is None:
            continue
        stripped_key = key[len(prefix):]
        if stripped_key:
            collected[stripped_key] = value
    return collected


def _json_dict(value: str | None) -> dict[str, Any]:
    """Parse an optional JSON object string.
    Args:
        value: Optional JSON object string.
    Returns:
        Parsed JSON dictionary.
    """
    if not value:
        return {}
    parsed = json.loads(value)
    if not isinstance(parsed, dict):
        raise ValueError("Expected a JSON object for extra normalization kwargs.")
    return parsed


def _build_parser() -> argparse.ArgumentParser:
    """Build the radiometric normalization CLI parser.
    Args:
        None.
    Returns:
        Configured argument parser.
    """
    parser = argparse.ArgumentParser(description="Run radiometric normalization via SpectralMatch.")
    parser.add_argument("--input-image", action="append", required=True, help="Input raster. Repeat for multiple images.")
    parser.add_argument(
        "--output-image",
        required=True,
        help="Merged normalized output image path.",
    )
    parser.add_argument("--temp-dir", help="Optional shared temp directory.")
    parser.add_argument("--delete-temp-dir", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--debug-logs", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--cache", default="auto")
    parser.add_argument("--nodata-value", type=float)
    parser.add_argument("--window-size", type=int, default=1024)
    parser.add_argument("--image-threads", default="auto")
    parser.add_argument("--io-threads", default="auto")
    parser.add_argument("--tile-threads", default="auto")
    parser.add_argument("--calculation-dtype", default="float32")
    parser.add_argument("--output-dtype")
    parser.add_argument("--save-as-cog", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument(
        "--extra-kwargs-json",
        help="Optional JSON object of additional SpectralMatch pipeline kwargs.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the radiometric normalization CLI.
    Args:
        argv: Optional command line arguments.
    Returns:
        Process exit code.
    """
    parser = _build_parser()
    args, unknown_args = parser.parse_known_args(argv)
    _apply_unknown_match_args(args, unknown_args)

    extra_kwargs = _json_dict(args.extra_kwargs_json)
    extra_kwargs.update(_collect_prefixed_kwargs(args, "match_"))

    output = radiometric_normalization(
        shared_input_images=args.input_image,
        shared_output_image_path=args.output_image,
        shared_temp_dir=args.temp_dir,
        delete_temp_dir=args.delete_temp_dir,
        shared_debug_logs=args.debug_logs,
        shared_cache=args.cache,
        shared_custom_nodata_value=args.nodata_value,
        shared_window_size=args.window_size,
        shared_image_threads=args.image_threads,
        shared_io_threads=args.io_threads,
        shared_tile_threads=args.tile_threads,
        shared_calculation_dtype=args.calculation_dtype,
        shared_output_dtype=args.output_dtype,
        shared_save_as_cog=args.save_as_cog,
        **extra_kwargs,
    )
    print(output)
    return 0


__all__ = ["main"]


if __name__ == "__main__":
    sys.exit(main())
