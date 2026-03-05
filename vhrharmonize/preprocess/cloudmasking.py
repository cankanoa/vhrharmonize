import os
from typing import Iterable, Optional, Sequence

import numpy as np
import rasterio
from affine import Affine
from scipy.ndimage import binary_dilation
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT


def _ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def _normalize_omnicloud_output(prediction: np.ndarray) -> np.ndarray:
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
    msg = str(exc).lower()
    return "cuda out of memory" in msg or "torch.outofmemoryerror" in msg


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
    inference_resolution_m: Optional[float] = 10.0,
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

    omnicloud_kwargs = dict(omnicloud_kwargs or {})
    # Keep large prediction mosaics off GPU by default to reduce VRAM pressure.
    omnicloud_kwargs.setdefault("mosaic_device", "cpu")

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
            # Match SpectralMatch sizing logic and GDAL buffered read behavior.
            from osgeo import gdal

            def _read_band_gdal_buffer(path: str, band_index: int, out_h: int, out_w: int) -> np.ndarray:
                ds = gdal.Open(path, gdal.GA_ReadOnly)
                if ds is None:
                    raise RuntimeError(f"Unable to open raster with GDAL: {path}")
                band = ds.GetRasterBand(int(band_index))
                raw = band.ReadRaster(
                    xoff=0,
                    yoff=0,
                    xsize=band.XSize,
                    ysize=band.YSize,
                    buf_xsize=out_w,
                    buf_ysize=out_h,
                    buf_type=gdal.GDT_Float32,
                )
                if raw is None:
                    raise RuntimeError(f"GDAL ReadRaster failed for band {band_index}: {path}")
                arr = np.frombuffer(raw, dtype=np.float32).reshape(out_h, out_w)
                return arr

            left, top = src.transform.c, src.transform.f
            px_w, px_h = abs(src.transform.a), abs(src.transform.e)
            right = left + src.width * src.transform.a
            bottom = top + src.height * src.transform.e
            out_w = max(1, int((right - left) / float(inference_resolution_m)))
            out_h = max(1, int((top - bottom) / float(inference_resolution_m)))
            red = _read_band_gdal_buffer(input_image_path, red_band_index, out_h, out_w)
            green = _read_band_gdal_buffer(input_image_path, green_band_index, out_h, out_w)
            nir = _read_band_gdal_buffer(input_image_path, nir_band_index, out_h, out_w)
            rgbn = np.stack([red, green, nir], axis=0)
            out_gt = (left, float(inference_resolution_m), 0.0, top, 0.0, -float(inference_resolution_m))
            rgbn_profile.update(
                height=out_h,
                width=out_w,
                transform=Affine.from_gdal(*out_gt),
            )
            print(
                f"Downsampling cloudmask inference input from ~{px_w:.3f} to {inference_resolution_m:.3f} map units"
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
            print(
                "OmniCloudMask CUDA OOM encountered. Retrying on CPU with conservative tiling."
            )
            raw_mask = predict_from_array(rgbn, **retry_kwargs)

        class_mask = _normalize_omnicloud_output(raw_mask)
        binary_mask = _build_cloud_binary_mask(class_mask, cloud_classes, buffer_pixels=buffer_pixels)

        _ensure_parent_dir(output_mask_path)
        with rasterio.open(output_mask_path, "w", **rgbn_profile) as dst:
            dst.write(binary_mask, 1)

    return output_mask_path


def apply_binary_cloud_mask_to_image(
    input_image_path: str,
    cloud_mask_path: str,
    output_image_path: str,
    *,
    cloud_mask_value: int = 1,
    output_nodata_value: Optional[float] = None,
    allow_mask_reprojection: bool = True,
) -> str:
    """
    Apply a binary cloud mask to an image by assigning NoData to masked pixels.
    """
    with rasterio.open(input_image_path) as src, rasterio.open(cloud_mask_path) as mask_src:
        mask_reader = mask_src
        if (mask_src.width, mask_src.height) != (src.width, src.height):
            if not allow_mask_reprojection:
                raise ValueError(
                    "Cloud mask shape "
                    f"{(mask_src.height, mask_src.width)} does not match image shape {(src.height, src.width)}"
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

        nodata_value = output_nodata_value if output_nodata_value is not None else src.nodata
        if nodata_value is None:
            raise ValueError(
                "No output nodata value was provided and input raster has no nodata set."
            )

        profile = src.profile.copy()
        profile.update(
            nodata=nodata_value,
            BIGTIFF="IF_SAFER",
        )
        _ensure_parent_dir(output_image_path)
        with rasterio.open(output_image_path, "w", **profile) as dst:
            for _, window in src.block_windows(1):
                mask_block = mask_reader.read(1, window=window)
                cloud_pixels = mask_block == cloud_mask_value
                data_block = src.read(window=window)
                data_block[:, cloud_pixels] = nodata_value
                dst.write(data_block, window=window)
        if mask_reader is not mask_src:
            mask_reader.close()

    return output_image_path
