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
from vhrharmonize.providers.standardized import StandardizedMetadata, materialize_scene_bounds
from vhrharmonize.providers.worldview import WorldViewMetadata


CORNERS = {
    "ULLon": -155.76202148, "ULLat": 20.03416250,
    "URLon": -155.60452034, "URLat": 20.01967459,
    "LRLon": -155.60207398, "LRLat": 19.88577036,
    "LLLon": -155.76204859, "LLLat": 19.90160157,
}


class SceneBoundsTests(unittest.TestCase):
    def test_corners_preserve_footprint_instead_of_rectangular_envelope(self):
        raw = {"BAND_C": dict(CORNERS), "BAND_B": dict(CORNERS)}
        original = deepcopy(raw)
        polygon = materialize_scene_bounds(raw)
        self.assertIsInstance(polygon, Polygon)
        self.assertEqual(list(polygon.exterior.coords)[:4], [
            (CORNERS[f"{corner}Lon"], CORNERS[f"{corner}Lat"])
            for corner in ("UL", "UR", "LR", "LL")
        ])
        self.assertLess(polygon.area, box(*polygon.bounds).area)
        self.assertEqual(raw, original)
        projected = materialize_scene_bounds(raw, epsg=32605)
        x, y = Transformer.from_crs(4326, 32605, always_xy=True).transform(CORNERS["ULLon"], CORNERS["ULLat"])
        self.assertAlmostEqual(projected.exterior.coords[0][0], x)
        self.assertAlmostEqual(projected.exterior.coords[0][1], y)

    def test_missing_bounds_do_not_prevent_metadata_loading(self):
        metadata = StandardizedMetadata.from_worldview_metadata(
            WorldViewMetadata(imd_file="scene.IMD", photo_basename=None, raw_metadata={"IMAGE_1": {"satId": "WV03"}})
        )
        self.assertEqual(metadata.sensor_id, "WV03")
        with self.assertRaisesRegex(ValueError, "no scene corner"):
            materialize_scene_bounds(metadata.source_metadata)

    def test_bad_and_inconsistent_corners_fail_when_requested(self):
        for values, message in [
            ({"ULLon": 0}, "Incomplete"),
            (dict(CORNERS, ULLat=float("nan")), "Invalid"),
            (dict(CORNERS, ULLat=100), "Invalid"),
            ({key: 0 for key in CORNERS}, "valid polygon"),
        ]:
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                materialize_scene_bounds({"BAND_C": values})
        with self.assertRaisesRegex(ValueError, "disagree"):
            materialize_scene_bounds({"BAND_C": CORNERS, "BAND_B": dict(CORNERS, ULLon=-150)})

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
