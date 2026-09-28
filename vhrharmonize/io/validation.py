"""Read-only validation of reusable raster and JSON outputs."""
from __future__ import annotations
import os
import math
import json
import struct
from typing import List, Dict
from osgeo import gdal
from tifffile import TiffFile, TiffFileError

def _is_gdal_raster_path(path: str) -> bool:
    """Return whether a path should be validated as a GDAL raster."""
    extension = os.path.splitext(path)[1].lower()
    return extension in {".tif", ".tiff", ".dat", ".img", ".vrt"}


def _tiff_truncation_reason(path: str) -> str | None:
    """Check strip/tile bounds in all TIFF image directories without decoding pixels."""
    try:
        with TiffFile(path) as tif:
            file_size = tif.filehandle.size
            pending = [tif.pages]
            while pending:
                for page in pending.pop():
                    for offset, byte_count in zip(page.dataoffsets, page.databytecounts):
                        if offset and byte_count and offset + byte_count > file_size:
                            return (
                                f"Truncated TIFF: offset {offset} + byte count {byte_count} "
                                f"exceeds file size {file_size} (IFD at {page.offset})"
                            )
                    if page.pages is not None:
                        pending.append(page.pages)
    except (OSError, ValueError, TiffFileError, NotImplementedError, struct.error):
        # Inconclusive TIFF inspection: let the existing GDAL checks decide.
        pass
    return None


def _gdal_raster_is_valid(path: str, *, validity_check_grid_size: int = 2048) -> tuple[bool, str | None]:
    """Check raster readability and TIFF strip/tile bounds.
    Args:
        path: Raster path to validate.
        validity_check_grid_size: Pixel sampling grid size; defaults to 2048. 0 disables pixel validity sampling.
    Returns:
        Tuple of validity and optional reason.
    """
    if not os.path.exists(path):
        return False, "missing"
    dataset = None
    try:
        dataset = gdal.OpenEx(path, gdal.OF_RASTER)
        if dataset is None:
            return False, "GDAL open failed"
        if dataset.GetDriver().ShortName == "GTiff":
            reason = _tiff_truncation_reason(path)
            if reason:
                return False, reason
        band_count = dataset.RasterCount
        width = dataset.RasterXSize
        height = dataset.RasterYSize
        if band_count < 1 or width < 1 or height < 1:
            return False, f"invalid raster shape bands={band_count} size={width}x{height}"

        read_w = min(width, 256)
        read_h = min(height, 256)
        offsets = [
            (0, 0),
            (max(0, width - read_w), max(0, height - read_h)),
        ]
        for band_index in range(1, band_count + 1):
            band = dataset.GetRasterBand(band_index)
            if band is None:
                return False, f"missing band {band_index}"
            for xoff, yoff in offsets:
                data = band.ReadRaster(xoff, yoff, read_w, read_h)
                if data is None:
                    return False, f"GDAL ReadRaster failed for band {band_index}"
        if validity_check_grid_size > 0 and not _gdal_raster_has_valid_sample(
            dataset,
            width,
            height,
            validity_check_grid_size,
        ):
            return False, "no finite valid pixels found in validity sample"
    except Exception as exc:
        return False, str(exc)
    finally:
        dataset = None
    return True, None


def _gdal_raster_has_valid_sample(dataset: gdal.Dataset, width: int, height: int, grid_size: int) -> bool:
    """Return whether any sampled raster pixel is finite and non-nodata."""
    if grid_size < 1:
        return True
    block_size = 512 if grid_size == 1 else grid_size
    for band_index in range(1, dataset.RasterCount + 1):
        band = dataset.GetRasterBand(band_index)
        if band is None:
            continue
        for yoff in range(0, height, block_size):
            read_h = min(block_size, height - yoff)
            sample_y = 0 if grid_size == 1 else min(grid_size // 2, read_h - 1)
            for xoff in range(0, width, block_size):
                read_w = min(block_size, width - xoff)
                sample_x = 0 if grid_size == 1 else min(grid_size // 2, read_w - 1)
                window_xoff = xoff if grid_size == 1 else xoff + sample_x
                window_yoff = yoff if grid_size == 1 else yoff + sample_y
                window_w = read_w if grid_size == 1 else 1
                window_h = read_h if grid_size == 1 else 1
                if _gdal_band_window_has_valid_pixel(
                    band,
                    window_xoff,
                    window_yoff,
                    window_w,
                    window_h,
                ):
                    return True
    return False


def _gdal_band_window_has_valid_pixel(
    band: gdal.Band,
    xoff: int,
    yoff: int,
    read_w: int,
    read_h: int,
) -> bool:
    """Return whether a sampled band window has a finite non-nodata pixel."""
    struct_format = _gdal_data_type_struct_format(band.DataType)
    if struct_format is None:
        return True
    data = band.ReadRaster(xoff, yoff, read_w, read_h)
    if data is None:
        return False
    mask = None
    mask_band = band.GetMaskBand()
    if mask_band is not None:
        mask = mask_band.ReadRaster(xoff, yoff, read_w, read_h, buf_type=gdal.GDT_Byte)

    nodata = band.GetNoDataValue()
    item_size = struct.calcsize(struct_format)
    if read_w == 1 and read_h == 1:
        return _sample_pixel_is_valid(
            struct.unpack_from(struct_format, data, 0)[0],
            nodata,
            mask[0] if mask else 255,
        )

    for pixel_index in range(0, read_w * read_h):
        mask_value = mask[pixel_index] if mask else 255
        value = struct.unpack_from(struct_format, data, pixel_index * item_size)[0]
        if _sample_pixel_is_valid(value, nodata, mask_value):
            return True
    return False


def _gdal_data_type_struct_format(data_type: int) -> str | None:
    """Return a struct format for scalar GDAL data types."""
    formats = {
        gdal.GDT_Byte: "=B",
        gdal.GDT_UInt16: "=H",
        gdal.GDT_Int16: "=h",
        gdal.GDT_UInt32: "=I",
        gdal.GDT_Int32: "=i",
        gdal.GDT_Float32: "=f",
        gdal.GDT_Float64: "=d",
    }
    if hasattr(gdal, "GDT_Int8"):
        formats[getattr(gdal, "GDT_Int8")] = "=b"
    if hasattr(gdal, "GDT_UInt64"):
        formats[getattr(gdal, "GDT_UInt64")] = "=Q"
    if hasattr(gdal, "GDT_Int64"):
        formats[getattr(gdal, "GDT_Int64")] = "=q"
    return formats.get(data_type)


def _sample_pixel_is_valid(value: object, nodata: float | int | None, mask_value: int = 255) -> bool:
    """Return whether a scalar sample is finite and not nodata."""
    if mask_value == 0:
        return False
    try:
        sample = float(value)
    except (TypeError, ValueError):
        return False
    if not math.isfinite(sample):
        return False
    if nodata is None:
        return True
    if isinstance(nodata, float) and math.isnan(nodata):
        return True
    return sample != nodata


def _json_file_is_valid(path: str) -> tuple[bool, str | None]:
    """Check that a file can be parsed as JSON without validating its contents."""
    try:
        with open(path, "r", encoding="utf-8") as handle:
            json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        return False, str(exc)
    return True, None


def _existing_output_failures(
    output_paths: List[str],
    *,
    check_validity: bool,
    validity_check_grid_size: int,
    log_to_console: bool,
    step: str,
    scene_basename: str | None = None,
) -> Dict[str, str]:
    """Inspect every expected output without changing files, including partial sets."""
    failures: Dict[str, str] = {}
    for output_path in output_paths:
        if not os.path.exists(output_path):
            failures[output_path] = "missing"
            continue
        if not check_validity:
            continue
        if os.path.splitext(output_path)[1].lower() == ".json":
            is_valid, reason = _json_file_is_valid(output_path)
        elif _is_gdal_raster_path(output_path):
            is_valid, reason = _gdal_raster_is_valid(
                output_path,
                validity_check_grid_size=validity_check_grid_size,
            )
        else:
            continue
        if not is_valid:
            failures[output_path] = reason or "invalid"
            if log_to_console:
                print(f"[{step}] Invalid output {output_path}: {reason}")
    return failures


def _existing_outputs_are_reusable(
    output_paths: List[str],
    *,
    check_validity: bool,
    validity_check_grid_size: int,
    log_to_console: bool,
    step: str,
    scene_basename: str | None = None,
) -> bool:
    """Return whether all outputs exist and pass enabled checks, without deleting files."""
    return bool(output_paths) and not _existing_output_failures(
        output_paths,
        check_validity=check_validity,
        validity_check_grid_size=validity_check_grid_size,
        log_to_console=log_to_console,
        step=step,
        scene_basename=scene_basename,
    )

