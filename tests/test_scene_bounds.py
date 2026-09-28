"""IMD footprints are materialized only by consumers that need geometry."""
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from pyproj import Transformer
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import Polygon, box

from vhrharmonize.io.geospatial import get_image_percentile_value
from vhrharmonize.io.metadata import materialize_geometry


CORNERS = {
    "ULLon": -155.76202148, "ULLat": 20.03416250,
    "URLon": -155.60452034, "URLat": 20.01967459,
    "LRLon": -155.60207398, "LRLat": 19.88577036,
    "LLLon": -155.76204859, "LLLat": 19.90160157,
}


class SceneBoundsTests(unittest.TestCase):
    def test_geojson_footprint_and_reprojection(self):
        geometry = {"type": "Polygon", "coordinates": [[[CORNERS[f"{corner}Lon"], CORNERS[f"{corner}Lat"]] for corner in ("UL", "UR", "LR", "LL", "UL")]]}
        original = deepcopy(geometry)
        polygon = materialize_geometry(geometry)
        self.assertLess(polygon.area, box(*polygon.bounds).area)
        self.assertEqual(geometry, original)
        projected = materialize_geometry(geometry, epsg=32605)
        x, y = Transformer.from_crs(4326, 32605, always_xy=True).transform(CORNERS["ULLon"], CORNERS["ULLat"])
        self.assertAlmostEqual(projected.exterior.coords[0][0], x)
        self.assertAlmostEqual(projected.exterior.coords[0][1], y)

    def test_dem_sampling_reprojects_in_memory_geometry(self):
        with TemporaryDirectory() as directory:
            path = str(Path(directory) / "dem.tif")
            with rasterio.open(path, "w", driver="GTiff", width=4, height=4, count=1,
                               dtype="float32", crs=3857, transform=from_origin(0, 4000, 1000, 1000)) as dst:
                dst.write(np.tile([10, 10, 100, 100], (4, 1)).astype("float32"), 1)
            lon, lat = Transformer.from_crs(3857, 4326, always_xy=True).transform(2000, 4000)
            self.assertEqual(get_image_percentile_value(path, mask=box(0, 0, lon, lat)), 10)
            self.assertEqual(get_image_percentile_value(path), 55)


if __name__ == "__main__":
    unittest.main()
