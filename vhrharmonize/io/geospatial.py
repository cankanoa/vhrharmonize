import warnings
from typing import Iterable, Optional
import rasterio
import rasterio.mask
import numpy as np

from pyproj import Transformer
from shapely.ops import transform
from shapely.geometry.base import BaseGeometry
from shapely.geometry import mapping


from vhrharmonize.preprocess.helpers import _logged_operation


@_logged_operation("overviews", inputs=("input_image_path",), outputs=("input_image_path",), allow_nested=True)
def calculate_raster_overviews(
    input_image_path: str,
    overview_scales: Optional[Iterable[int]],
    *,
    resampling: str = "nearest",
    log_to_console: bool = False,
    scene_basename: str | None = None,
    scene_index: int = 1,
    scene_total: int = 1,
    ) -> str:
    """Build internal raster overviews.
    Args:
        input_image_path: Input raster path.
        overview_scales: Requested overview decimation factors.
        resampling: Rasterio overview resampling method name.
    Returns:
        Input raster path after overview creation.
    """
    factors = []
    seen = set()
    for value in overview_scales or []:
        factor = int(value)
        if factor <= 1 or factor in seen:
            continue
        seen.add(factor)
        factors.append(factor)

    if not factors:
        return input_image_path

    with rasterio.open(input_image_path, "r+") as dataset:
        dataset.build_overviews(factors, rasterio.enums.Resampling[resampling])
        dataset.update_tags(ns="rio_overview", resampling=resampling)
    return input_image_path


def get_image_percentile_value(
    input_image_path: str,
    percentile: float = 50.0,
    mask: Optional[BaseGeometry] = None,
    mask_crs: object = 4326,
    ) -> float:
    """Compute a raster percentile value.
    Args:
        input_image_path: Input raster path.
        percentile: Requested percentile in the inclusive range 0 to 100.
        mask: Optional Shapely geometry used to limit sampled pixels.
        mask_crs: CRS of the mask geometry; defaults to EPSG:4326.
    Returns:
        Percentile value from valid raster pixels.
    """
    if percentile < 0 or percentile > 100:
        raise ValueError("percentile must be between 0 and 100.")

    with rasterio.open(input_image_path) as src:
        nodata = src.nodata
        collected = []

        if mask is not None:
            if src.crs is None:
                raise ValueError("Raster CRS is required to sample using a mask geometry")
            geometry = transform(Transformer.from_crs(mask_crs, src.crs, always_xy=True).transform, mask)
            geometries = [mapping(geometry)]
        else:
            geometries = None

        for i in range(1, src.count + 1):
            if geometries:
                band_array, _ = rasterio.mask.mask(src, geometries, indexes=i, filled=False)
            else:
                band_array = src.read(i, masked=True)

            if nodata is not None:
                band_array = np.ma.masked_equal(band_array, nodata)

            if band_array.mask.all():
                continue

            valid = band_array.compressed() if np.ma.isMaskedArray(band_array) else band_array.ravel()
            if valid.size:
                collected.append(valid)

    if not collected:
        warnings.warn(
            f"Warning: Percentile value in file '{input_image_path}' could not be found. "
            "This could be because the mask is outside of the bounds of this image."
        )
        return float("nan")

    values = np.concatenate(collected)
    return float(np.percentile(values, percentile))
