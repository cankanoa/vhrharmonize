"""Raster tiling utilities."""

import os
from typing import List, Optional

import numpy as np
import rasterio
from rasterio.windows import Window


def tile_image(
    input_image_path: str,
    output_dir: str,
    tile_width: int,
    tile_height: Optional[int] = None,
    buffer: int = 0,
    output_prefix: Optional[str] = None,
    skip_empty: bool = False,
    empty_check_band_index: int = 1,
) -> List[str]:
    """Split a raster into tiles with optional pixel overlap between adjacent tiles.

    Args:
        input_image_path: Input raster path.
        output_dir: Directory where output tiles will be written.
        tile_width: Tile width in pixels.
        tile_height: Tile height in pixels. Defaults to `tile_width`.
        buffer: Pixel overlap between adjacent tiles. Must be >= 0 and < tile size.
        output_prefix: Optional filename prefix. Defaults to input basename.
        skip_empty: If True, skip tiles with no valid pixels in `empty_check_band_index`.
        empty_check_band_index: 1-based band index used for empty-tile detection.

    Returns:
        List of written tile file paths.
    """
    if tile_width <= 0:
        raise ValueError("tile_width must be > 0")
    if tile_height is None:
        tile_height = tile_width
    if tile_height <= 0:
        raise ValueError("tile_height must be > 0")
    if buffer < 0:
        raise ValueError("buffer must be >= 0")
    if buffer >= tile_width or buffer >= tile_height:
        raise ValueError("buffer must be smaller than both tile_width and tile_height")
    if empty_check_band_index < 1:
        raise ValueError("empty_check_band_index must be >= 1")

    os.makedirs(output_dir, exist_ok=True)
    step_x = tile_width - buffer
    step_y = tile_height - buffer
    written: List[str] = []

    with rasterio.open(input_image_path) as src:
        if empty_check_band_index > src.count:
            raise ValueError(
                f"empty_check_band_index={empty_check_band_index} exceeds band count={src.count}"
            )

        if not output_prefix:
            output_prefix = os.path.splitext(os.path.basename(input_image_path))[0]

        for y_off in range(0, src.height, step_y):
            row_idx = y_off // step_y
            for x_off in range(0, src.width, step_x):
                col_idx = x_off // step_x
                width = min(tile_width, src.width - x_off)
                height = min(tile_height, src.height - y_off)
                window = Window(col_off=x_off, row_off=y_off, width=width, height=height)

                if skip_empty:
                    mask = src.read_masks(empty_check_band_index, window=window)
                    if not np.any(mask > 0):
                        continue

                tile_profile = src.profile.copy()
                tile_profile.update(
                    width=width,
                    height=height,
                    transform=src.window_transform(window),
                    tiled=False,
                )

                tile_data = src.read(window=window)
                output_name = f"{output_prefix}_r{row_idx:04d}_c{col_idx:04d}.tif"
                output_path = os.path.join(output_dir, output_name)

                with rasterio.open(output_path, "w", **tile_profile) as dst:
                    dst.write(tile_data)

                written.append(output_path)

    return written


__all__ = ["tile_image"]
