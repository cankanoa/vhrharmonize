from pathlib import Path
from tempfile import TemporaryDirectory
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import patch
from concurrent.futures import ProcessPoolExecutor
import sys

import geopandas as gpd
from osgeo import gdal, osr
from shapely.geometry import box, MultiPolygon

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
                    imd_file="source.IMD", basename="source",
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
        with patch.object(seamline_metadata, "_iter_seamline_metadata_results") as calculate, patch.object(seamline_metadata, "_open_seamline_metadata_writer") as save:
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

    def test_results_are_committed_as_they_arrive_and_resume_after_failure(self):
        with patch.object(seamline_metadata, "_valid_data_polygon_from_image", return_value=box(0, 0, 1, 1)):
            self.write([self.state("a.tif")])

        def results(*args, **kwargs):
            yield "b.tif", box(1, 1, 2, 2)
            # A separate reader can see b before the next worker result arrives.
            self.assertEqual(list(self.read().image_basename), ["a.tif", "b.tif"])
            raise RuntimeError("worker failed")

        states = [self.state(name) for name in ("a.tif", "b.tif", "c.tif")]
        with patch.object(seamline_metadata, "_iter_seamline_metadata_results", results):
            with self.assertRaisesRegex(RuntimeError, "worker failed"):
                self.write(states)
        with patch.object(seamline_metadata, "_valid_data_polygon_from_image", return_value=box(2, 2, 3, 3)) as calculate:
            self.write(states)
        calculate.assert_called_once_with("c.tif", eight_connected=True)
        self.assertEqual(list(self.read().image_basename), ["a.tif", "b.tif", "c.tif"])

    def test_out_of_order_results_keep_numeric_fields_and_multipart_geometry(self):
        states = [self.state("a.tif"), self.state("b.tif")]
        states[0].scene.mul_image.standardized_metadata.to_dict = lambda: {"score": 1.5}
        states[1].scene.mul_image.standardized_metadata.to_dict = lambda: {"score": None}
        geometry = MultiPolygon([box(0, 0, 1, 1), box(2, 2, 3, 3)])

        def results(*args, **kwargs):
            yield "b.tif", geometry
            yield "a.tif", box(4, 4, 5, 5)

        with patch.object(seamline_metadata, "_iter_seamline_metadata_results", results):
            self.write(states)
        result = self.read()
        self.assertEqual(list(result.image_basename), ["b.tif", "a.tif"])
        self.assertEqual(result.score.iloc[1], 1.5)
        self.assertTrue(result.geometry.iloc[0].equals(geometry))

    def test_real_process_pool_polygonizes_images(self):
        states = []
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(6635)
        for name in ("a", "b", "c"):
            image_path = str(Path(self.temp.name) / f"{name}.tif")
            dataset = gdal.GetDriverByName("GTiff").Create(image_path, 4, 4, 1, gdal.GDT_Byte)
            dataset.SetGeoTransform((0, 1, 0, 4, 0, -1))
            dataset.SetProjection(srs.ExportToWkt())
            dataset.GetRasterBand(1).SetNoDataValue(0)
            dataset.GetRasterBand(1).Fill(1)
            dataset = None
            states.append(self.state(image_path))
        self.kwargs["concurrent_processing"] = 2
        with patch.object(seamline_metadata, "ProcessPoolExecutor", wraps=ProcessPoolExecutor) as executor:
            self.write(states)
        executor.assert_called_once_with(max_workers=2)
        result = self.read()
        self.assertEqual(set(result.image_basename), {"a.tif", "b.tif", "c.tif"})
        self.assertTrue(all(geometry.equals(box(0, 0, 4, 4)) for geometry in result.geometry))

    def test_package_bounds_workers_use_imd_and_reuse_existing_without_corners(self):
        from test_scene_bounds import CORNERS

        states = [self.state("a.tif"), self.state("b.tif")]
        for state in states:
            state.scene.mul_image.standardized_metadata.source_metadata = {"BAND_C": CORNERS}
        self.kwargs.update(footprint_source="package_bounds", epsg=4326, concurrent_processing=2)
        self.write(states)
        result = self.read()
        self.assertEqual(set(result.image_basename), {"a.tif", "b.tif"})
        expected = seamline_metadata.materialize_scene_bounds({"BAND_C": CORNERS})
        self.assertTrue(all(geometry.equals(expected) for geometry in result.geometry))
        # Reused records need neither an input raster nor valid IMD bounds.
        for state in states:
            state.scene.mul_image.standardized_metadata.source_metadata = {}
        with patch.object(seamline_metadata, "_iter_seamline_metadata_results") as workers:
            self.write(states)
        workers.assert_not_called()

    def test_dask_submits_only_missing_images_and_writes_in_completion_order(self):
        from concurrent.futures import Future

        with patch.object(seamline_metadata, "_valid_data_polygon_from_image", return_value=box(0, 0, 1, 1)):
            self.write([self.state("a.tif")])
        submitted = []

        class Client:
            closed = False

            def submit(self, function, *args):
                submitted.append(args[0])
                future = Future()
                future.set_result(function(*args))
                return future

            def close(self):
                self.closed = True

        client = Client()
        dask = ModuleType("dask")
        distributed = ModuleType("dask.distributed")
        distributed.as_completed = lambda futures: reversed(futures)
        dask.distributed = distributed
        self.kwargs.update(concurrent_processing_backend="dask", dask_scheduler_address="tcp://scheduler:8786")
        with patch.dict(sys.modules, {"dask": dask, "dask.distributed": distributed}), patch.object(
            seamline_metadata, "_make_dask_client", return_value=client,
        ) as make_client, patch.object(seamline_metadata, "_valid_data_polygon_from_image", return_value=box(1, 1, 2, 2)):
            self.write([self.state(name) for name in ("a.tif", "b.tif", "c.tif")])
        self.assertEqual(submitted, ["b.tif", "c.tif"])
        self.assertEqual(list(self.read().image_basename), ["a.tif", "c.tif", "b.tif"])
        self.assertEqual(make_client.call_args.args[0].dask_scheduler_address, "tcp://scheduler:8786")
        self.assertTrue(client.closed)



if __name__ == "__main__":
    unittest.main()
