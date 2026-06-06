"""Radiometric normalization via the upstream SpectralMatch pipeline API."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

from spectralmatch.chain import pipeline as spectralmatch_pipeline
from vhrharmonize.preprocess.helpers import log


def _normalize_input_images(shared_input_images: Iterable[str]) -> list[str]:
    """Normalize shared input image paths.
    Args:
        shared_input_images: Iterable of candidate raster paths.
    Returns:
        Normalized non-empty raster path list.
    """
    normalized = [str(path) for path in shared_input_images if str(path)]
    if not normalized:
        raise ValueError("shared_input_images must contain at least one raster path.")
    return normalized


def radiometric_normalization(
    shared_input_images: Iterable[str],
    shared_output_image_path: str,
    *,
    method: str = "spectralmatch",
    shared_temp_dir: str | None = None,
    delete_temp_dir: bool = True,
    shared_debug_logs: bool = False,
    shared_cache: str = "auto",
    shared_custom_nodata_value: float | int | None = None,
    shared_window_size: int = 1024,
    shared_image_threads: int | str = "auto",
    shared_io_threads: int | str = "auto",
    shared_tile_threads: int | str = "auto",
    shared_calculation_dtype: str = "float32",
    shared_output_dtype: str | None = None,
    shared_save_as_cog: bool = False,
    log_to_console: bool = False,
    **kwargs: Any,
) -> str:
    """Run radiometric normalization with SpectralMatch.
    Args:
        shared_input_images: Input raster paths for normalization.
        shared_output_image_path: Output normalized raster path.
        method: Requested radiometric normalization method.
        shared_temp_dir: Optional SpectralMatch temp directory.
        delete_temp_dir: Whether SpectralMatch should delete its temp directory.
        shared_debug_logs: Whether SpectralMatch should emit debug logs.
        shared_cache: SpectralMatch cache setting to pass through.
        shared_custom_nodata_value: Optional custom nodata override.
        shared_window_size: SpectralMatch shared window size.
        shared_image_threads: SpectralMatch image thread setting.
        shared_io_threads: SpectralMatch IO thread setting.
        shared_tile_threads: SpectralMatch tile thread setting.
        shared_calculation_dtype: SpectralMatch calculation dtype.
        shared_output_dtype: Optional SpectralMatch output dtype override.
        shared_save_as_cog: Whether to save as a COG.
        log_to_console: Whether to emit console logs.
        kwargs: Additional SpectralMatch keyword arguments.
    Returns:
        Output normalized raster path.
    """
    method_norm = method.strip().lower()
    if method_norm != "spectralmatch":
        raise ValueError(f"Unsupported radiometric normalization method: {method}")

    output_path = Path(shared_output_image_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    input_images = _normalize_input_images(shared_input_images)
    log(
        f"Running radiometric normalization inputs={len(input_images)} output={output_path.name}",
        enabled=log_to_console,
        step="radiometric",
    )

    pipeline_kwargs = {
        "shared_input_images": input_images,
        "shared_output_image_path": str(output_path),
        "delete_temp_dir": delete_temp_dir,
        "shared_debug_logs": shared_debug_logs,
        "shared_cache": shared_cache,
        "shared_window_size": shared_window_size,
        "shared_image_threads": shared_image_threads,
        "shared_io_threads": shared_io_threads,
        "shared_tile_threads": shared_tile_threads,
        "shared_calculation_dtype": shared_calculation_dtype,
        "shared_save_as_cog": shared_save_as_cog,
    }
    if shared_temp_dir is not None:
        pipeline_kwargs["shared_temp_dir"] = shared_temp_dir
    if shared_custom_nodata_value is not None:
        pipeline_kwargs["shared_custom_nodata_value"] = shared_custom_nodata_value
    if shared_output_dtype is not None:
        pipeline_kwargs["shared_output_dtype"] = shared_output_dtype

    for key, value in kwargs.items():
        if value is not None:
            pipeline_kwargs[key] = value

    result = spectralmatch_pipeline(**pipeline_kwargs)

    if output_path.exists():
        log(f"Wrote output {output_path.name}", enabled=log_to_console, step="radiometric")
        return str(output_path)
    if isinstance(result, str) and result:
        return result
    if hasattr(result, "get") and callable(result.get):
        output_value = result.get("output")
        if output_value:
            return str(output_value)
    if hasattr(result, "shared_output_image_path") and getattr(result, "shared_output_image_path"):
        return str(getattr(result, "shared_output_image_path"))
    if hasattr(result, "output_image_path") and getattr(result, "output_image_path"):
        return str(getattr(result, "output_image_path"))

    raise RuntimeError("spectralmatch pipeline did not produce the requested output image.")


__all__ = ["radiometric_normalization"]
