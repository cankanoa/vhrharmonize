import os
from typing import Optional

import numpy as np
import orthority as oty
import rasterio

from vhrharmonize.preprocess.helpers import log


def pansharpen_image(
    input_low_resolution_path: str,
    input_high_resolution_path: str,
    output_image_path: str,
    change_nodata_value: Optional[float] = None,
    log_to_console: bool = False,
    scene_basename: str | None = None,
    ) -> None:
    """Pansharpen a multispectral raster.
    Args:
        input_low_resolution_path: Input multispectral raster path.
        input_high_resolution_path: Input panchromatic raster path.
        output_image_path: Output pansharpened raster path.
        change_nodata_value: Optional nodata value to apply after pansharpening.
        log_to_console: Whether to emit console logs.
        scene_basename: Optional scene basename for log prefixes.
    Returns:
        None.
    """
    log(
        f"Running pansharpen mul={os.path.basename(input_low_resolution_path)} pan={os.path.basename(input_high_resolution_path)}",
        enabled=log_to_console,
        step="pansharpen",
        scene_basename=scene_basename,
    )
    pan_sharp = oty.PanSharpen(input_high_resolution_path, input_low_resolution_path)
    pan_sharp.process(output_image_path, write_mask=False, overwrite=True)

    if change_nodata_value is not None:
        _change_nodata_value(
            output_image_path,
            change_nodata_value,
            -32768,
            log_to_console=log_to_console,
            scene_basename=scene_basename,
        )
    log(
        f"Wrote output {os.path.basename(output_image_path)}",
        enabled=log_to_console,
        step="pansharpen",
        scene_basename=scene_basename,
    )


def _change_nodata_value(
    input_image_path: str,
    new_nodata_value: float,
    old_nodata: float,
    *,
    log_to_console: bool = False,
    scene_basename: str | None = None,
    ) -> None:
    """Replace nodata values in-place.
    Args:
        input_image_path: Raster path to update in-place.
        new_nodata_value: New nodata value to assign.
        old_nodata: Old nodata value to replace.
        log_to_console: Whether to emit console logs.
        scene_basename: Optional scene basename for log prefixes.
    Returns:
        None.
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
        log(
            "Updated nodata values",
            enabled=log_to_console,
            step="pansharpen",
            scene_basename=scene_basename,
        )
