import numpy as np
import orthority as oty
import rasterio

from vhrharmonize.logging_utils import log


def pansharpen_image(
    input_low_resolution_path,
    input_high_resolution_path,
    output_image_path,
    change_nodata_value = None,
    log_to_console: bool = False,
    ):
    """Pansharpen a multispectral raster with a higher-resolution panchromatic raster."""
    log("Running pansharpen", enabled=log_to_console, step="pansharpen")
    pan_sharp = oty.PanSharpen(input_high_resolution_path, input_low_resolution_path)
    pan_sharp.process(output_image_path, write_mask=False, overwrite=True)

    if change_nodata_value is not None:
        _change_nodata_value(
            output_image_path,
            change_nodata_value,
            -32768,
            log_to_console=log_to_console,
        )
    log("Wrote output", enabled=log_to_console, step="pansharpen")


def _change_nodata_value(
    input_image_path,
    new_nodata_value,
    old_nodata,
    *,
    log_to_console: bool = False,
    ):

    """
    Replaces all pixels equal to `old_nodata` with `new_nodata_value` in-place.
    Uses block-wise IO to avoid loading very large rasters into memory.

    :param input_image_path: Path to the original raster file.
    :param new_nodata_value: Value to use as the new NoData.
    :param old_nodata: The pixel value currently used as NoData.
    """
    replaced_any = False
    with rasterio.open(input_image_path, "r+") as src:
        for band_idx in range(1, src.count + 1):
            for _, window in src.block_windows(band_idx):
                band_data = src.read(band_idx, window=window)
                if np.issubdtype(band_data.dtype, np.floating):
                    mask = np.isclose(band_data, old_nodata)
                else:
                    mask = band_data == old_nodata
                if np.any(mask):
                    band_data[mask] = new_nodata_value
                    src.write(band_data, band_idx, window=window)
                    replaced_any = True

        # GTiff stores a single dataset nodata (TIFFTAG_GDAL_NODATA).
        src.nodata = new_nodata_value

    if replaced_any:
        log("Updated nodata values", enabled=log_to_console, step="pansharpen")
