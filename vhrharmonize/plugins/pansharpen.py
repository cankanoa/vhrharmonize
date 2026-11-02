import os
from typing import Optional
from types import FunctionType, MethodType, SimpleNamespace

import numpy as np
import rasterio

from .base import FunctionPlugin
from vhrharmonize.io.logging import _log, _logged_operation


def _configure_pan_threads(pan_sharp, tile_threads):
    """Limit Orthority's hard-coded pools without changing shared module globals.

    Orthority has no thread-count argument. Give just this instance's two pool-
    creating methods private globals for CPU count and Rasterio's GDAL settings.
    The original algorithms and other concurrently running instances are intact.
    """
    threads = (os.cpu_count() or 1) if tile_threads is None else tile_threads
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError("tile_threads must be a positive integer or None.")

    def env(*args, **kwargs):
        return rasterio.Env(*args, **{**kwargs, "GDAL_NUM_THREADS": threads})

    for name in ("_get_stats", "process"):
        method = getattr(type(pan_sharp), name)
        namespace = {
            **method.__globals__,
            "os": SimpleNamespace(**{**vars(os), "cpu_count": lambda: threads}),
            "rio": SimpleNamespace(**{**vars(rasterio), "Env": env}),
        }
        limited = FunctionType(method.__code__, namespace, method.__name__,
                               method.__defaults__, method.__closure__)
        limited.__kwdefaults__ = method.__kwdefaults__
        setattr(pan_sharp, name, MethodType(limited, pan_sharp))
    for profile in pan_sharp._profiles.values():
        profile["num_threads"] = threads
    return threads


@_logged_operation(
    "pansharpen",
    inputs=("input_low_resolution_path", "input_high_resolution_path"),
    outputs=("output_image_path",),
    worker_progress=True,
)
def pansharpen_image(
    input_low_resolution_path: str,
    input_high_resolution_path: str,
    output_image_path: str,
    custom_nodata_value: Optional[float] = None,
    log_to_console: bool = False,
    scene_basename: str | None = None,
    tile_threads: int | None = 1,
) -> None:
    """Pansharpen a multispectral raster.
    Args:
        input_low_resolution_path: Input multispectral raster path.
        input_high_resolution_path: Input panchromatic raster path.
        output_image_path: Output pansharpened raster path.
        custom_nodata_value: NoData value written during pansharpening;
            None uses Orthority's dtype-specific default.
        log_to_console: Whether to emit console logs.
        scene_basename: Optional scene basename for log prefixes.
        tile_threads: Threads per image for Orthority workers and GDAL; None uses all CPUs.
    Returns:
        None.
    """
    if not output_image_path:
        raise ValueError("output_image_path must be an explicit output path")
    output_image_path = str(output_image_path)
    import orthority as oty

    pan_sharp = oty.PanSharpen(input_high_resolution_path, input_low_resolution_path)
    threads = _configure_pan_threads(pan_sharp, tile_threads)
    with rasterio.Env(GDAL_NUM_THREADS=threads):
        if custom_nodata_value is None:
            pan_sharp.process(output_image_path, write_mask=False, overwrite=True, build_ovw=False)
        else:
            from orthority.common import create_profile

            # Reuse Orthority's calculated intersection grid and output defaults.
            # An open writer keeps its profile when passed to PanSharpen.process().
            grid = pan_sharp._profiles["pan_to_pan"]
            with rasterio.open(input_low_resolution_path) as src:
                dtype = np.dtype(src.dtypes[0])
                count = sum(color != rasterio.enums.ColorInterp.alpha for color in src.colorinterp)
            nodata = float(custom_nodata_value)
            if dtype.kind in "iu":
                limits = np.iinfo(dtype)
                valid_nodata = (
                    np.isfinite(nodata) and nodata.is_integer()
                    and float(limits.min) <= nodata <= float(limits.max)
                )
            else:
                limits = np.finfo(dtype)
                valid_nodata = np.isnan(nodata) or (
                    np.isfinite(nodata) and float(limits.min) <= nodata <= float(limits.max)
                )
            if not valid_nodata:
                raise ValueError(f"custom_nodata_value is not representable in {dtype.name}.")
            profile, _ = create_profile(
                driver="gtiff", shape=(count, grid["height"], grid["width"]),
                dtype=dtype.name, write_mask=False,
            )
            profile.update(crs=grid["crs"], transform=grid["transform"], nodata=nodata)
            with rasterio.open(output_image_path, "w", **profile) as dst:
                pan_sharp.process(dst, write_mask=False, overwrite=True, build_ovw=False)
    _log(
        f"Wrote output {os.path.basename(output_image_path)}",
        enabled=log_to_console,
        step="pansharpen",
        scene_basename=scene_basename,
    )


class Pansharpen(FunctionPlugin):
    tile_threads_parameter = "tile_threads"
    input_dependency_paths = frozenset({"input_high_resolution_path", "input_low_resolution_path"})
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset({"output_image_path"})
    output_dependency_paths = output_path_resolution_paths
    output_parent_creation_paths = output_path_resolution_paths
    output_collision_check_paths = output_path_resolution_paths
    output_reuse_paths = output_path_resolution_paths
    output_validation_paths = output_path_resolution_paths
    output_invalid_removal_paths = output_path_resolution_paths
    output_overview_calculation_paths = frozenset({"output_image_path"})
    output_temporary_cleanup_paths = output_path_resolution_paths
    output_hpc_staging_paths = output_path_resolution_paths
    output_hpc_download_paths = output_path_resolution_paths

    target = "vhrharmonize.plugins.pansharpen:pansharpen_image"
