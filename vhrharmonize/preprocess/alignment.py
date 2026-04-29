"""Pairwise image alignment wrapper around coregix."""

from __future__ import annotations

from dataclasses import dataclass
import os

from vhrharmonize.preprocess.helpers import log


@dataclass
class AlignmentResult:
    """Summary of pairwise alignment execution."""

    output_image_path: str


def _extract_output_path(result) -> str:
    for attribute_name in ("output_image_path", "aligned_image_path", "path"):
        value = getattr(result, attribute_name, None)
        if isinstance(value, str) and value:
            return value
    raise ValueError("coregix align_image_pair result did not include an output image path.")


def align_image_pair(
    moving_image_path: str,
    fixed_image_path: str,
    output_image_path: str,
    *,
    band_index: int = 0,
    moving_band_index: int | None = None,
    fixed_band_index: int | None = None,
    moving_nodata: float | None = None,
    fixed_nodata: float | None = None,
    output_nodata: float | None = None,
    min_valid_fraction: float = 0.01,
    temp_dir: str | None = None,
    keep_temp_dir: bool = False,
    split_factor: int = 2,
    clip_fixed_to_moving: bool = False,
    output_on_moving_grid: bool = True,
    trim_edge_invalid: bool = False,
    edge_trim_depth: int = 8,
    edge_trim_detection_band_index: int = 0,
    edge_trim_invalid_below: float | None = None,
    edge_trim_invalid_above: float | None = None,
    enforce_mutual_valid_mask: bool = False,
    use_edge_proxies: bool = True,
    solve_resolution: float | None = None,
    log_to_console: bool = False,
) -> AlignmentResult:
    """Align a moving image to a fixed image using coregix."""
    if band_index < 0:
        raise ValueError("band_index must be >= 0.")
    if moving_band_index is not None and moving_band_index < 0:
        raise ValueError("moving_band_index must be >= 0.")
    if fixed_band_index is not None and fixed_band_index < 0:
        raise ValueError("fixed_band_index must be >= 0.")
    if min_valid_fraction <= 0 or min_valid_fraction > 1:
        raise ValueError("min_valid_fraction must be in (0, 1].")
    if split_factor < 0:
        raise ValueError("split_factor must be >= 0.")
    if edge_trim_depth <= 0:
        raise ValueError("edge_trim_depth must be > 0.")
    if edge_trim_detection_band_index < 0:
        raise ValueError("edge_trim_detection_band_index must be >= 0.")
    if solve_resolution is not None and solve_resolution <= 0:
        raise ValueError("solve_resolution must be > 0 when provided.")

    try:
        from coregix import align_image_pair as coregix_align_image_pair
    except ImportError as exc:
        raise RuntimeError(
            "coregix is not installed. Install alignment extras with `pip install -e \".[align]\"`."
        ) from exc

    log(
        f"Running alignment moving={os.path.basename(moving_image_path)} fixed={os.path.basename(fixed_image_path)} split_factor={split_factor}",
        enabled=log_to_console,
        step="alignment",
    )
    result = coregix_align_image_pair(
        moving_image_path=moving_image_path,
        fixed_image_path=fixed_image_path,
        output_image_path=output_image_path,
        band_index=band_index,
        moving_band_index=moving_band_index,
        fixed_band_index=fixed_band_index,
        moving_nodata=moving_nodata,
        fixed_nodata=fixed_nodata,
        output_nodata=output_nodata,
        min_valid_fraction=min_valid_fraction,
        temp_dir=temp_dir,
        keep_temp_dir=keep_temp_dir,
        clip_fixed_to_moving=clip_fixed_to_moving,
        output_on_moving_grid=output_on_moving_grid,
        split_factor=split_factor,
        trim_edge_invalid=trim_edge_invalid,
        edge_trim_depth=edge_trim_depth,
        edge_trim_detection_band_index=edge_trim_detection_band_index,
        edge_trim_invalid_below=edge_trim_invalid_below,
        edge_trim_invalid_above=edge_trim_invalid_above,
        enforce_mutual_valid_mask=enforce_mutual_valid_mask,
        use_edge_proxies=use_edge_proxies,
        solve_resolution=solve_resolution,
        log_to_console=log_to_console,
    )
    output_path = _extract_output_path(result)
    log(f"Wrote output {os.path.basename(output_path)}", enabled=log_to_console, step="alignment")
    return AlignmentResult(output_image_path=output_path)


__all__ = ["AlignmentResult", "align_image_pair"]
