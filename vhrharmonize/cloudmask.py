import os
from typing import Iterable, Optional, Sequence

import numpy as np
import rasterio
from scipy.ndimage import binary_dilation


def _ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def _normalize_omnicloud_output(prediction: np.ndarray) -> np.ndarray:
    pred = np.asarray(prediction)
    if pred.ndim == 3:
        # Support confidence/class-first output by taking argmax if needed.
        pred = np.argmax(pred, axis=0)
    if pred.ndim != 2:
        raise ValueError(f"Unexpected omnicloudmask output shape: {pred.shape}")
    return pred.astype(np.uint8)


def _build_cloud_binary_mask(
    class_mask: np.ndarray,
    cloud_classes: Sequence[int],
    buffer_pixels: int = 0,
) -> np.ndarray:
    binary_mask = np.isin(class_mask, np.array(list(cloud_classes), dtype=np.uint8))
    if buffer_pixels > 0:
        binary_mask = binary_dilation(binary_mask, iterations=buffer_pixels)
    return binary_mask.astype(np.uint8)


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
) -> str:
    """
    Run OmniCloudMask on an input raster and write a binary cloud mask.

    Output mask values:
    - 1: cloud/cloud-shadow (or buffered cloud zone)
    - 0: clear
    """
    from omnicloudmask import predict_from_array

    if red_band_index < 1 or green_band_index < 1 or nir_band_index < 1:
        raise ValueError("Band indexes must be 1-based and >= 1.")

    omnicloud_kwargs = omnicloud_kwargs or {}

    with rasterio.open(input_image_path) as src:
        if max(red_band_index, green_band_index, nir_band_index) > src.count:
            raise ValueError(
                f"Requested band index exceeds band count ({src.count}) for {input_image_path}"
            )
        rgbn = src.read([red_band_index, green_band_index, nir_band_index]).astype(np.float32)
        class_mask = _normalize_omnicloud_output(predict_from_array(rgbn, **omnicloud_kwargs))
        binary_mask = _build_cloud_binary_mask(class_mask, cloud_classes, buffer_pixels=buffer_pixels)

        profile = src.profile.copy()
        profile.update(count=1, dtype="uint8", nodata=255)

        _ensure_parent_dir(output_mask_path)
        with rasterio.open(output_mask_path, "w", **profile) as dst:
            dst.write(binary_mask, 1)

    return output_mask_path


def apply_binary_cloud_mask_to_image(
    input_image_path: str,
    cloud_mask_path: str,
    output_image_path: str,
    *,
    cloud_mask_value: int = 1,
    output_nodata_value: Optional[float] = None,
) -> str:
    """
    Apply a binary cloud mask to an image by assigning NoData to masked pixels.
    """
    with rasterio.open(input_image_path) as src:
        data = src.read()
        profile = src.profile.copy()
        nodata_value = output_nodata_value if output_nodata_value is not None else src.nodata
        if nodata_value is None:
            raise ValueError(
                "No output nodata value was provided and input raster has no nodata set."
            )

    with rasterio.open(cloud_mask_path) as mask_src:
        mask = mask_src.read(1)
        if mask.shape != data.shape[1:]:
            raise ValueError(
                f"Cloud mask shape {mask.shape} does not match image shape {data.shape[1:]}"
            )

    cloud_pixels = mask == cloud_mask_value
    data[:, cloud_pixels] = nodata_value

    profile.update(nodata=nodata_value)
    _ensure_parent_dir(output_image_path)
    with rasterio.open(output_image_path, "w", **profile) as dst:
        dst.write(data)

    return output_image_path
