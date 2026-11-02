import os
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import rasterio
from affine import Affine
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT

from .base import FunctionPlugin
from vhrharmonize.io.logging import _log, _logged_operation
from vhrharmonize.io.progress import raster_windows


@dataclass(frozen=True)
class CloudMaskResult:
    """Summary of cloud mask generation and application."""

    output_mask_path: str
    output_raster_path: str
    mask_pixel_count: int


@_logged_operation(
    "cloud_mask", inputs=("input_image_path", "inference_image_path"), outputs=("output_raster_path", "output_mask_path"),
    worker_progress=True,
)
def cloudmask_raster(
    input_image_path: str,
    output_raster_path: str,
    output_mask_path: str,
    *,
    red_band_index: int,
    green_band_index: int,
    nir_band_index: int,
    inference_image_path: Optional[str] = None,
    cloud_classes: Sequence[int] = (1, 2, 3),
    buffer_pixels: int = 0,
    omnicloud_kwargs: Optional[dict] = None,
    inference_resolution_m: Optional[float] = 10.0,
    cloud_mask_value: int = 1,
    custom_nodata_value: Optional[float] = None,
    allow_mask_reprojection: bool = True,
    log_to_console: bool = False,
    scene_basename: str | None = None,
) -> CloudMaskResult:
    """Create a cloud mask and apply it to ``input_image_path``.

    ``inference_image_path`` optionally supplies a separate raster for prediction,
    such as orthorectified MUL imagery when masking a pansharpened product. Band
    indices and inference resolution apply to that raster. By default, prediction
    uses ``input_image_path`` too. The saved mask stays on the inference grid and
    is resampled with nearest neighbor onto the application grid when needed.
    """
    if not output_raster_path:
        raise ValueError("output_raster_path must be an explicit output path")
    output_raster_path = str(output_raster_path)
    if not output_mask_path:
        raise ValueError("output_mask_path must be an explicit output path")
    output_mask_path = str(output_mask_path)
    _log(
        f"Processing {os.path.basename(input_image_path)}",
        enabled=log_to_console,
        step="cloudmask",
        scene_basename=scene_basename,
    )
    create_cloud_mask_with_omnicloudmask(
        inference_image_path if inference_image_path is not None else input_image_path,
        output_mask_path,
        red_band_index=red_band_index,
        green_band_index=green_band_index,
        nir_band_index=nir_band_index,
        cloud_classes=cloud_classes,
        buffer_pixels=buffer_pixels,
        omnicloud_kwargs=omnicloud_kwargs,
        inference_resolution_m=inference_resolution_m,
        log_to_console=log_to_console,
        scene_basename=scene_basename,
    )
    apply_binary_cloud_mask_to_image(
        input_image_path,
        output_mask_path,
        output_raster_path,
        cloud_mask_value=cloud_mask_value,
        custom_nodata_value=custom_nodata_value,
        allow_mask_reprojection=allow_mask_reprojection,
        log_to_console=log_to_console,
        scene_basename=scene_basename,
    )

    mask_pixel_count = 0
    with rasterio.open(output_mask_path) as src:
        for _, window in raster_windows(src, desc="Counting cloud pixels"):
            block = src.read(1, window=window)
            mask_pixel_count += int((block == cloud_mask_value).sum())

    return CloudMaskResult(
        output_mask_path=output_mask_path,
        output_raster_path=output_raster_path,
        mask_pixel_count=mask_pixel_count,
    )


def _ensure_parent_dir(path: str) -> None:
    """Create a parent directory when needed.
    Args:
        path: File path whose parent should exist.
    Returns:
        None.
    """
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def _normalize_omnicloud_output(prediction: np.ndarray) -> np.ndarray:
    """Normalize OmniCloudMask output to a 2D class mask.
    Args:
        prediction: Raw OmniCloudMask prediction array.
    Returns:
        Normalized 2D uint8 class mask.
    """
    pred = np.asarray(prediction)
    pred = np.squeeze(pred)
    if pred.ndim == 3:
        # Class-probability cubes can be class-first or class-last.
        # Keep singleton class dimensions as plain masks, otherwise reduce.
        if pred.shape[0] == 1:
            pred = pred[0]
        elif pred.shape[-1] == 1:
            pred = pred[..., 0]
        elif np.issubdtype(pred.dtype, np.floating):
            if pred.shape[0] <= pred.shape[-1]:
                pred = np.argmax(pred, axis=0)
            else:
                pred = np.argmax(pred, axis=-1)
        else:
            raise ValueError(f"Unexpected 3D omnicloudmask output shape: {pred.shape}")
    if pred.ndim != 2:
        raise ValueError(f"Unexpected omnicloudmask output shape: {pred.shape}")
    return pred.astype(np.uint8)


def _is_cuda_oom_error(exc: Exception) -> bool:
    """Return whether an exception is a CUDA OOM error.
    Args:
        exc: Exception to inspect.
    Returns:
        True when the exception indicates CUDA OOM.
    """
    msg = str(exc).lower()
    return "cuda out of memory" in msg or "torch.outofmemoryerror" in msg


def _build_cloud_binary_mask(
    class_mask: np.ndarray,
    cloud_classes: Sequence[int],
    buffer_pixels: int = 0,
) -> np.ndarray:
    """Build a binary cloud mask from class values.
    Args:
        class_mask: OmniCloudMask class mask.
        cloud_classes: Class ids treated as cloudy.
        buffer_pixels: Optional dilation buffer in pixels.
    Returns:
        Binary uint8 cloud mask.
    """
    binary_mask = np.isin(class_mask, np.array(list(cloud_classes), dtype=np.uint8))
    if buffer_pixels > 0:
        from scipy.ndimage import binary_dilation

        binary_mask = binary_dilation(binary_mask, iterations=buffer_pixels)
    return binary_mask.astype(np.uint8)


@_logged_operation("cloud_mask", inputs=("input_image_path",), outputs=("output_mask_path",), worker_progress=True)
def create_cloud_mask_with_omnicloudmask(
    input_image_path: str,
    output_mask_path: str,
    red_band_index: int,
    green_band_index: int,
    nir_band_index: int,
    *,
    cloud_classes: Sequence[int] = (1, 2, 3),
    buffer_pixels: int = 0,
    omnicloud_kwargs: Optional[dict] = None,
    inference_resolution_m: Optional[float] = 10.0,
    log_to_console: bool = False,
    scene_basename: str | None = None,
) -> str:
    """
    Run OmniCloudMask on an input raster and write a binary cloud mask.

    Output mask values:
    - 1: cloud/cloud-shadow (or buffered cloud zone)
    - 0: clear
    """

    if red_band_index < 1 or green_band_index < 1 or nir_band_index < 1:
        raise ValueError("Band indexes must be 1-based and >= 1.")
    if buffer_pixels < 0:
        raise ValueError("buffer_pixels must be >= 0.")
    from omnicloudmask import predict_from_array

    omnicloud_kwargs = dict(omnicloud_kwargs or {})
    # Keep large prediction mosaics off GPU by default to reduce VRAM pressure.
    omnicloud_kwargs.setdefault("mosaic_device", "cpu")
    _log(
        f"OmniCloudMask setup resolution={inference_resolution_m}m buffer={buffer_pixels} classes={list(cloud_classes)}",
        enabled=log_to_console,
        step="cloudmask",
        scene_basename=scene_basename,
    )

    with rasterio.open(input_image_path) as src:
        if max(red_band_index, green_band_index, nir_band_index) > src.count:
            raise ValueError(
                f"Requested band index exceeds band count ({src.count}) for {input_image_path}"
            )
        rgbn_profile = src.profile.copy()
        rgbn_profile.update(count=1, dtype="uint8", nodata=255)

        # OmniCloudMask is designed for ~10-50 m inputs. For VHR scenes, downsample
        # inference input by default and later reproject mask back when applying.
        rgbn = None
        if (
            inference_resolution_m is not None
            and inference_resolution_m > 0
            and src.crs is not None
            and src.res is not None
            and src.res[0] > 0
            and src.res[1] > 0
            and src.crs.is_projected
            and (src.res[0] < inference_resolution_m or src.res[1] < inference_resolution_m)
        ):
            # Resample the complete source extent and use the exact same grid
            # for the mask. Fractional source pixels must not shift its footprint.
            out_w = max(1, int(np.ceil(src.width * src.res[0] / inference_resolution_m)))
            out_h = max(1, int(np.ceil(src.height * src.res[1] / inference_resolution_m)))
            nodata = src.nodata if src.nodata is not None else 0.0
            with WarpedVRT(
                src,
                transform=src.transform * Affine.scale(src.width / out_w, src.height / out_h),
                width=out_w, height=out_h, dtype="float32", nodata=nodata,
                resampling=Resampling.average,
            ) as inference:
                rgbn = inference.read([red_band_index, green_band_index, nir_band_index])
                rgbn_profile.update(
                    height=out_h, width=out_w, transform=inference.transform,
                )
            _log(
                f"Downsampling inference input to {inference_resolution_m:.3f} map units",
                enabled=log_to_console,
                step="cloudmask",
                scene_basename=scene_basename,
            )
        else:
            rgbn = src.read([red_band_index, green_band_index, nir_band_index]).astype(np.float32)

        if "no_data_value" not in omnicloud_kwargs:
            src_nodata = src.nodata
            omnicloud_kwargs["no_data_value"] = float(src_nodata) if src_nodata is not None else 0.0
        try:
            raw_mask = predict_from_array(rgbn, **omnicloud_kwargs)
        except Exception as exc:
            if not _is_cuda_oom_error(exc):
                raise
            retry_kwargs = dict(omnicloud_kwargs)
            retry_kwargs["inference_device"] = "cpu"
            retry_kwargs["mosaic_device"] = "cpu"
            retry_kwargs.setdefault("batch_size", 1)
            retry_kwargs.setdefault("patch_size", 768)
            retry_kwargs.setdefault("patch_overlap", 192)
            _log(
                "Retrying on CPU after CUDA OOM",
                enabled=log_to_console,
                step="cloudmask",
                scene_basename=scene_basename,
            )
            raw_mask = predict_from_array(rgbn, **retry_kwargs)

        class_mask = _normalize_omnicloud_output(raw_mask)
        binary_mask = _build_cloud_binary_mask(
            class_mask, cloud_classes, buffer_pixels=buffer_pixels
        )

        _ensure_parent_dir(output_mask_path)
        with rasterio.open(output_mask_path, "w", **rgbn_profile) as dst:
            dst.write(binary_mask, 1)

    return output_mask_path


@_logged_operation(
    "apply_cloud_mask",
    inputs=("input_image_path", "cloud_mask_path"),
    outputs=("output_image_path",),
    worker_progress=True,
)
def apply_binary_cloud_mask_to_image(
    input_image_path: str,
    cloud_mask_path: str,
    output_image_path: str,
    *,
    cloud_mask_value: int = 1,
    custom_nodata_value: Optional[float] = None,
    allow_mask_reprojection: bool = True,
    log_to_console: bool = False,
    scene_basename: str | None = None,
) -> str:
    """
    Apply a binary cloud mask to an image by assigning NoData to masked pixels.
    """
    with rasterio.open(input_image_path) as src, rasterio.open(cloud_mask_path) as mask_src:
        mask_reader = mask_src
        if (
            (mask_src.width, mask_src.height) != (src.width, src.height)
            or mask_src.transform != src.transform
            or mask_src.crs != src.crs
        ):
            if not allow_mask_reprojection:
                raise ValueError(
                    "Cloud mask grid does not match image grid "
                    f"(mask shape {(mask_src.height, mask_src.width)}, "
                    f"image shape {(src.height, src.width)}); "
                    "allow_mask_reprojection is False."
                )
            if src.crs is None or mask_src.crs is None:
                raise ValueError(
                    "Cannot reproject cloud mask because source image or mask has no CRS."
                )
            mask_reader = WarpedVRT(
                mask_src,
                crs=src.crs,
                transform=src.transform,
                width=src.width,
                height=src.height,
                resampling=Resampling.nearest,
            )
            _log(
                "Reprojecting cloud mask onto image grid",
                enabled=log_to_console,
                step="cloudmask",
                scene_basename=scene_basename,
            )

        custom_nodata_value = custom_nodata_value if custom_nodata_value is not None else src.nodata
        if custom_nodata_value is None:
            raise ValueError(
                "No output nodata value was provided and input raster has no nodata set."
            )

        profile = src.profile.copy()
        profile.update(
            nodata=custom_nodata_value,
            BIGTIFF="IF_SAFER",
        )
        _ensure_parent_dir(output_image_path)
        with rasterio.open(output_image_path, "w", **profile) as dst:
            for _, window in raster_windows(src, desc="Applying cloud mask"):
                mask_block = mask_reader.read(1, window=window)
                cloud_pixels = mask_block == cloud_mask_value
                data_block = src.read(window=window)
                data_block[:, cloud_pixels] = custom_nodata_value
                dst.write(data_block, window=window)
        if mask_reader is not mask_src:
            mask_reader.close()

    _log(
        "Wrote mask and masked raster",
        enabled=log_to_console,
        step="cloudmask",
        scene_basename=scene_basename,
    )
    return output_image_path


__all__ = [
    "CloudMask",
    "CloudMaskResult",
    "cloudmask_raster",
    "create_cloud_mask_with_omnicloudmask",
    "apply_binary_cloud_mask_to_image",
]


class CloudMask(FunctionPlugin):
    input_dependency_paths = frozenset({"input_image_path", "inference_image_path"})
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset({"output_mask_path", "output_raster_path"})
    output_dependency_paths = output_path_resolution_paths
    output_parent_creation_paths = output_path_resolution_paths
    output_collision_check_paths = output_path_resolution_paths
    output_reuse_paths = output_path_resolution_paths
    output_validation_paths = output_path_resolution_paths
    output_invalid_removal_paths = output_path_resolution_paths
    output_overview_calculation_paths = frozenset({"output_mask_path", "output_raster_path"})
    output_temporary_cleanup_paths = output_path_resolution_paths
    output_hpc_staging_paths = output_path_resolution_paths
    output_hpc_download_paths = output_path_resolution_paths

    target = "vhrharmonize.plugins.cloud_mask:cloudmask_raster"
