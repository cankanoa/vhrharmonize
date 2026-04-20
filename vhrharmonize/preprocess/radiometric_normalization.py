"""Radiometric normalization via the upstream SpectralMatch pipeline API."""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any, Callable, Iterable

from vhrharmonize.preprocess.helpers import log


def _load_spectralmatch_pipeline() -> Callable[..., Any]:
    """Import the upstream ``spectralmatch.pipeline`` callable."""
    candidate_paths = [
        ("spectralmatch.pipeline", "pipeline"),
        ("spectralmatch", "pipeline"),
    ]
    errors = []
    for module_name, attr_name in candidate_paths:
        try:
            module = import_module(module_name)
            pipeline = getattr(module, attr_name, None)
            if callable(pipeline):
                return pipeline
        except Exception as exc:
            errors.append(f"{module_name}.{attr_name}: {exc}")
    raise RuntimeError(
        "spectralmatch pipeline function could not be imported. "
        "Install extras with `pip install -e \".[spectralmatch]\"` and ensure the package "
        "exposes `pipeline(...)`."
        + (f" Tried: {'; '.join(errors)}" if errors else "")
    )


def _normalize_input_images(shared_input_images: Iterable[str]) -> list[str]:
    normalized = [str(path) for path in shared_input_images if str(path)]
    if not normalized:
        raise ValueError("shared_input_images must contain at least one raster path.")
    return normalized


def _build_pipeline_kwargs(
    *,
    shared_input_images: Iterable[str],
    shared_output_image_path: str,
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
    **kwargs: Any,
) -> dict[str, Any]:
    pipeline_kwargs = {
        "shared_input_images": _normalize_input_images(shared_input_images),
        "shared_output_image_path": str(shared_output_image_path),
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

    return pipeline_kwargs


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
    """Run radiometric normalization using the upstream SpectralMatch pipeline.

    See: https://spectralmatch.github.io/spectralmatch/api/pipeline/
    """
    method_norm = method.strip().lower()
    if method_norm != "spectralmatch":
        raise ValueError(f"Unsupported radiometric normalization method: {method}")

    output_path = Path(shared_output_image_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log("Running radiometric normalization", enabled=log_to_console, step="radiometric")

    pipeline_kwargs = _build_pipeline_kwargs(
        shared_input_images=shared_input_images,
        shared_output_image_path=str(output_path),
        shared_temp_dir=shared_temp_dir,
        delete_temp_dir=delete_temp_dir,
        shared_debug_logs=shared_debug_logs,
        shared_cache=shared_cache,
        shared_custom_nodata_value=shared_custom_nodata_value,
        shared_window_size=shared_window_size,
        shared_image_threads=shared_image_threads,
        shared_io_threads=shared_io_threads,
        shared_tile_threads=shared_tile_threads,
        shared_calculation_dtype=shared_calculation_dtype,
        shared_output_dtype=shared_output_dtype,
        shared_save_as_cog=shared_save_as_cog,
        **kwargs,
    )

    pipeline = _load_spectralmatch_pipeline()
    result = pipeline(**pipeline_kwargs)

    if output_path.exists():
        log("Wrote output", enabled=log_to_console, step="radiometric")
        return str(output_path)
    if isinstance(result, str) and result:
        return result
    if hasattr(result, "shared_output_image_path") and getattr(result, "shared_output_image_path"):
        return str(getattr(result, "shared_output_image_path"))
    if hasattr(result, "output_image_path") and getattr(result, "output_image_path"):
        return str(getattr(result, "output_image_path"))

    raise RuntimeError("spectralmatch pipeline did not produce the requested output image.")


__all__ = ["radiometric_normalization"]
