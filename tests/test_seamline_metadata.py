from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import geopandas as gpd
from shapely.geometry import box

from vhrharmonize.preprocess import seamline_metadata


class SeamlineMetadataTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.output_path = str(Path(self.temp.name) / "seamline_metadata.gpkg")
        self.kwargs = dict(
            layer="footprints", image_field_name="source_image",
            footprint_source="calculate_bounds",
            calculate_bounds_eight_connected=True, epsg=6635,
            run_from_existing_check_validity=True,
        )

    def state(self, image_path):
        return SimpleNamespace(
            current_files=[image_path],
            scene=SimpleNamespace(
                primary_basename=Path(image_path).stem, scene_id="scene",
                catalog_id="catalog",
                mul_image=SimpleNamespace(
                    shp_file="source.shp", imd_file="source.IMD", basename="source",
                    standardized_metadata=SimpleNamespace(to_dict=lambda: {"sun_elevation": 40.0}),
                ),
            ),
        )

    def write(self, states):
        return seamline_metadata.write_seamline_metadata_gpkg(
            states, self.output_path, **self.kwargs,
        )

    def read(self):
        return gpd.read_file(self.output_path, layer="footprints")

    def test_add_missing_and_preserve_existing_records_and_other_layers(self):
        with patch.object(seamline_metadata, "_valid_data_polygon_from_image", return_value=box(0, 0, 1, 1)):
            self.write([self.state("/old/a.tif"), self.state("/old/extra.tif")])
        existing = self.read()
        gpd.GeoDataFrame({"name": ["keep"]}, geometry=[box(0, 0, 2, 2)], crs=6635).to_file(
            self.output_path, layer="other", driver="GPKG",
        )

        with patch.object(seamline_metadata, "_valid_data_polygon_from_image", return_value=box(2, 2, 3, 3)) as calculate:
            self.write([self.state("/new/a.tif"), self.state("/new/b.tif"), self.state("/new/b.tif")])
        calculate.assert_called_once_with("/new/b.tif", eight_connected=True)
        result = self.read()
        self.assertEqual(list(result.image_basename), ["a.tif", "extra.tif", "b.tif"])
        self.assertTrue(result.iloc[:2].reset_index(drop=True).equals(existing))
        self.assertEqual(list(gpd.read_file(self.output_path, layer="other").name), ["keep"])

    def test_all_matching_does_not_calculate_or_write(self):
        with patch.object(seamline_metadata, "_valid_data_polygon_from_image", return_value=box(0, 0, 1, 1)):
            self.write([self.state("/old/a.tif")])
        with patch.object(seamline_metadata, "_valid_data_polygon_from_image") as calculate, patch.object(gpd.GeoDataFrame, "to_file") as save:
            self.write([self.state("/new/a.tif")])
        calculate.assert_not_called()
        save.assert_not_called()

    def test_disabled_validation_replaces_output(self):
        with patch.object(seamline_metadata, "_valid_data_polygon_from_image", return_value=box(0, 0, 1, 1)):
            self.write([self.state("a.tif")])
            self.kwargs["run_from_existing_check_validity"] = False
            self.write([self.state("b.tif")])
        self.assertEqual(list(self.read().image_basename), ["b.tif"])

    def test_calculation_failure_preserves_existing_file(self):
        with patch.object(seamline_metadata, "_valid_data_polygon_from_image", return_value=box(0, 0, 1, 1)):
            self.write([self.state("a.tif")])
        with patch.object(seamline_metadata, "_valid_data_polygon_from_image", side_effect=RuntimeError("read failed")):
            with self.assertRaisesRegex(RuntimeError, "read failed"):
                self.write([self.state("a.tif"), self.state("b.tif")])
        self.assertEqual(list(self.read().image_basename), ["a.tif"])


if __name__ == "__main__":
    unittest.main()
