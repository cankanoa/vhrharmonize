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


def _gdal_raster_is_valid(
    path: str, *, validity_raster_grid: int | bool = 255,
    validity_raster_check_stop: str = "readable",
) -> tuple[bool, str | None]:
    """Inspect metadata and optionally read a regular pixel grid in every band.

    False disables data/mask reads. One reads every pixel in bounded windows;
    N > 1 reads pixels at (0, N, 2N, ...) along both axes. ``readable`` visits
    the entire grid, including masks, allowing NoData. ``has_valid_data`` stops
    at the first finite, unmasked, non-NoData sample and rejects an empty grid.
    TIFF byte ranges (including overview IFDs) are checked in either mode.
    """
    if not os.path.exists(path):
        return False, "missing"
    dataset = None
    try:
        with gdal.ExceptionMgr(useExceptions=True):
            dataset = gdal.OpenEx(path, gdal.OF_RASTER)
            if dataset is None:
                return False, "GDAL open failed"
            if dataset.GetDriver().ShortName == "GTiff":
                reason = _tiff_truncation_reason(path)
                if reason:
                    return False, reason
            width, height, bands = dataset.RasterXSize, dataset.RasterYSize, dataset.RasterCount
            if min(width, height, bands) < 1:
                return False, f"invalid raster shape bands={bands} size={width}x{height}"
            for index in range(1, bands + 1):
                if dataset.GetRasterBand(index) is None:
                    return False, f"missing band {index}"
            if validity_raster_grid is not False:
                found = _read_raster_grid(dataset, width, height, validity_raster_grid,
                                          validity_raster_check_stop)
                if validity_raster_check_stop == "has_valid_data" and not found:
                    return False, "no finite valid pixels found in validity sample"
    except Exception as exc:
        return False, str(exc)
    finally:
        dataset = None
    return True, None


def _read_raster_grid(dataset, width, height, spacing, stop):
    """Read all requested samples unless the caller explicitly permits early success."""
    stride = 512 if spacing == 1 else spacing
    for index in range(1, dataset.RasterCount + 1):
        band = dataset.GetRasterBand(index)
        for y in range(0, height, stride):
            for x in range(0, width, stride):
                w, h = (min(512, width - x), min(512, height - y)) if spacing == 1 else (1, 1)
                valid = _gdal_band_window_has_valid_pixel(band, x, y, w, h,
                                                         check_values=stop == "has_valid_data")
                if stop == "has_valid_data" and valid:
                    return True
    return False


def _gdal_band_window_has_valid_pixel(
    band: gdal.Band,
    xoff: int,
    yoff: int,
    read_w: int,
    read_h: int,
    *, check_values: bool = True,
) -> bool:
    """Return whether a sampled band window has a finite non-nodata pixel."""
    struct_format = "=d"
    data = band.ReadRaster(xoff, yoff, read_w, read_h, buf_type=gdal.GDT_Float64)
    if data is None or len(data) != read_w * read_h * 8:
        raise IOError(f"Incomplete raster read at ({xoff}, {yoff})")
    mask = None
    mask_band = band.GetMaskBand()
    if mask_band is not None:
        mask = mask_band.ReadRaster(xoff, yoff, read_w, read_h, buf_type=gdal.GDT_Byte)
        if mask is None or len(mask) != read_w * read_h:
            raise IOError(f"Incomplete mask read at ({xoff}, {yoff})")
    if not check_values:
        return True

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
    validity_raster_grid: int | bool = 255,
    validity_raster_check_stop: str = "readable",
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
                validity_raster_grid=validity_raster_grid,
                validity_raster_check_stop=validity_raster_check_stop,
            )
        elif os.path.splitext(output_path)[1].lower() == ".gpkg":
            is_valid, reason = _gpkg_is_valid(output_path)
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
    validity_raster_grid: int | bool = 255,
    validity_raster_check_stop: str = "readable",
    log_to_console: bool,
    step: str,
    scene_basename: str | None = None,
) -> bool:
    """Return whether all outputs exist and pass enabled checks, without deleting files."""
    return bool(output_paths) and not _existing_output_failures(
        output_paths,
        check_validity=check_validity,
        validity_raster_grid=validity_raster_grid,
        validity_raster_check_stop=validity_raster_check_stop,
        log_to_console=log_to_console,
        step=step,
        scene_basename=scene_basename,
    )



def _gpkg_is_valid(path: str) -> tuple[bool, str | None]:
    """Check SQLite integrity and the feature tables referenced by GeoPackage metadata."""
    import sqlite3
    from contextlib import closing
    from pathlib import Path

    try:
        with closing(sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro", uri=True)) as db:
            results = db.execute("PRAGMA quick_check").fetchall()
            if results != [("ok",)]:
                return False, f"GeoPackage integrity check failed: {results}"
            tables = {row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            for name, in db.execute("SELECT table_name FROM gpkg_contents"):
                if name not in tables:
                    return False, f"GeoPackage missing table: {name}"
        with gdal.ExceptionMgr(useExceptions=True):
            ds = gdal.OpenEx(path, gdal.OF_VECTOR)
            if ds is None or ds.GetLayerCount() < 1:
                return False, "GeoPackage has no readable vector layers"
            for layer in ds:
                layer.GetFeatureCount()
                layer.ResetReading()
                for feature in layer:
                    geometry = feature.GetGeometryRef()
                    if geometry is not None:
                        geometry.ExportToWkb()
            ds = None
    except Exception as exc:
        return False, str(exc)
    return True, None
