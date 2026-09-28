"""Pairwise image alignment wrapper around coregix."""

from __future__ import annotations

from dataclasses import dataclass
import os
from typing import Any

from .base import FunctionPlugin
from vhrharmonize.io.logging import _log, _logged_operation


def coregix_align_image_pair(**kwargs):
    """Load the optional registration backend only when alignment runs."""
    from coregix import align_image_pair

    return align_image_pair(**kwargs)


@dataclass
class AlignmentResult:
    """Summary of pairwise alignment execution."""

    output_image_path: str


def _extract_output_path(result: Any) -> str:
    """Extract an output image path from a coregix result.
    Args:
        result: Returned object from coregix alignment.
    Returns:
        Extracted output image path.
    """
    for attribute_name in ("output_image_path", "aligned_image_path", "path"):
        value = getattr(result, attribute_name, None)
        if isinstance(value, str) and value:
            return value
    raise ValueError("coregix align_image_pair result did not include an output image path.")


@_logged_operation(
    "alignment", inputs=("moving_image_path", "fixed_image_path"), outputs=("output_image_path",)
)
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
    custom_nodata_value: float | None = None,
    min_valid_fraction: float = 0.01,
    temp_dir: str | None = None,
    delete_temp_dir: bool = True,
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
    scene_basename: str | None = None,
) -> AlignmentResult:
    """Align a moving image to a fixed image.
    Args:
        moving_image_path: Moving image path to transform.
        fixed_image_path: Fixed reference image path.
        output_image_path: Output aligned image path.
        band_index: Default band index used when band-specific overrides are not provided.
        moving_band_index: Optional moving-image band index override.
        fixed_band_index: Optional fixed-image band index override.
        moving_nodata: Optional moving-image nodata override.
        fixed_nodata: Optional fixed-image nodata override.
        custom_nodata_value: Optional output nodata override.
        min_valid_fraction: Minimum valid overlap fraction required for registration.
        temp_dir: Optional temp directory for coregix intermediates.
        delete_temp_dir: Whether to delete the coregix temp directory after alignment.
        split_factor: Chunking factor used by coregix.
        clip_fixed_to_moving: Whether to clip the fixed image to the moving bounds before registration.
        output_on_moving_grid: Whether to write the aligned result on the moving grid.
        trim_edge_invalid: Whether to trim invalid alignment edge artifacts.
        edge_trim_depth: Edge trim depth in pixels.
        edge_trim_detection_band_index: Detection band index used for edge trimming.
        edge_trim_invalid_below: Optional lower invalid-value threshold for edge trimming.
        edge_trim_invalid_above: Optional upper invalid-value threshold for edge trimming.
        enforce_mutual_valid_mask: Whether to enforce mutual valid masks during alignment.
        use_edge_proxies: Whether to use edge proxies during alignment.
        solve_resolution: Optional solve resolution override.
        log_to_console: Whether to emit console logs.
        scene_basename: Optional scene basename for log prefixes.
    Returns:
        Alignment result summary.
    """
    if not output_image_path:
        raise ValueError("output_image_path must be an explicit output path")
    output_image_path = str(output_image_path)
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

    _log(
        f"Processing moving={os.path.basename(moving_image_path)} fixed={os.path.basename(fixed_image_path)} split_factor={split_factor}",
        enabled=log_to_console,
        step="alignment",
        scene_basename=scene_basename,
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
        output_nodata=custom_nodata_value,
        min_valid_fraction=min_valid_fraction,
        temp_dir=temp_dir,
        keep_temp_dir=not delete_temp_dir,
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
    _log(
        f"Wrote output {os.path.basename(output_path)}",
        enabled=log_to_console,
        step="alignment",
        scene_basename=scene_basename,
    )
    return AlignmentResult(output_image_path=output_path)


__all__ = [
    "Alignment",
    "AlignmentResult",
    "align_image_pair",
]


class Alignment(FunctionPlugin):
    input_dependency_paths = frozenset({"fixed_image_path", "moving_image_path"})
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset({"output_image_path"})
    output_dependency_paths = output_path_resolution_paths
    output_target_paths = output_path_resolution_paths
    output_parent_creation_paths = output_path_resolution_paths
    output_collision_check_paths = output_path_resolution_paths
    output_reuse_paths = output_path_resolution_paths
    output_validation_paths = output_path_resolution_paths
    output_invalid_removal_paths = output_path_resolution_paths
    output_overview_calculation_paths = frozenset({"output_image_path"})
    output_temporary_cleanup_paths = output_path_resolution_paths
    output_hpc_staging_paths = output_path_resolution_paths
    output_hpc_download_paths = output_path_resolution_paths
    output_context_checkpoint_paths = frozenset({"output_image_path"})

    target = "vhrharmonize.plugins.alignment:align_image_pair"
