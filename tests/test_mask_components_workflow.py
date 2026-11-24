"""Numbered mask products remain real workflow outputs and merge inputs."""
from pathlib import Path

import geopandas as gpd
import pytest
from osgeo import gdal
from shapely.geometry import MultiPolygon, box

from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import import_settings, stage, transfer


def recipe(tmp_path, make_test_raster):
    source = make_test_raster(tmp_path / "source.tif", width=32, height=32, crs="EPSG:32605")
    vector = tmp_path / "mask.gpkg"
    gpd.GeoDataFrame(geometry=[MultiPolygon([box(2, -10, 10, -2), box(20, -24, 28, -16)])],
                     crs=32605).to_file(vector, driver="GPKG")
    return {
        "shared": {"plugin": "shared", "core:run": True, "core:log_to_console": False,
                   "core:show_progress": False, "core:save_statistics_path": None,
                   "core:concurrent_processing": 2, "core:tile_threads": 1,
                   "core:delete_temp_steps_proactively": False,
                   "core:run_from_existing": False,
                   "param:window_scales": [2]},
        "files": import_settings(source, tmp_path),
        "mask": {"plugin": "mask_rasters", "core:run": True,
                 "core:require_outputs": "param:output_images",
                 "core:calculate_overviews": True,
                 "param:input_images": "collect:mul",
                 "param:output_images": [str(tmp_path / "output/masked.tif")],
                 "core:requires": str(vector),
                 "param:vector_mask": ["include", str(vector)],
                 "param:seperate_multipolygons": True,
                 "param:resume_from_outputs": "validate",
                 "const:masked_images": "returned:$"},
    }


@pytest.mark.parametrize("generated_input", [False, True])
def test_split_mask_then_merge_with_real_outputs(tmp_path, make_test_raster, generated_input):
    config = recipe(tmp_path, make_test_raster)
    if generated_input:
        mask = config.pop("mask")
        config["copy"] = {"plugin": "band_math", "core:run": True,
                          "param:input_images": "collect:mul", "param:threshold_math": "B1",
                          "var:mul": "expr:const.temp_dir & '/prepared.tif'",
                          "param:output_images": "collect:mul"}
        config["mask"] = mask
    config["merge"] = {"plugin": "merge_rasters", "core:run": True,
                       "core:require_outputs": "param:output_image_path",
                       "param:input_images": "const:masked_images",
                       "param:output_image_path": str(tmp_path / "output/mosaic.tif")}
    workflow = Workflow(config)
    mask_node = next(n for n in workflow.nodes if n.step["name"] == "mask")
    workflow.run()
    outputs = [str(tmp_path / f"output/masked_{i}.tif") for i in (1, 2)]
    assert workflow.context["const"]["masked_images"] == outputs
    assert not (tmp_path / "output/masked.tif").exists()
    assert mask_node.paths("output_validation_paths") == outputs
    for filename in outputs:
        with gdal.Open(filename) as ds:
            assert ds.RasterXSize == ds.RasterYSize == 8
            assert ds.GetRasterBand(1).GetOverviewCount() == 1
            assert (ds.ReadAsArray() == 1).all()
    with gdal.Open(str(tmp_path / "output/mosaic.tif")) as ds:
        assert (ds.ReadAsArray() == 1).sum() == 128
    mtimes = [Path(p).stat().st_mtime_ns for p in outputs]
    Workflow(config).run()
    assert [Path(p).stat().st_mtime_ns for p in outputs] == mtimes
    # Cleanup must retain the resolved file list after an input is removed.
    if generated_input:
        (tmp_path / "temp/prepared.tif").unlink()
        assert mask_node.paths("output_temporary_cleanup_paths") == outputs


def test_split_output_staging_downloads_numbered_files(tmp_path, make_test_raster):
    config = recipe(tmp_path, make_test_raster)
    staged, uploads, downloads = stage(config, tmp_path)
    outputs = [str(tmp_path / f"output/masked_{i}.tif") for i in (1, 2)]
    assert set(downloads) == set(outputs)
    transfer(uploads)
    Workflow(staged).run()
    for filename in outputs:
        with gdal.Open(downloads[filename]) as ds:
            assert (ds.ReadAsArray() == 1).all()


def test_split_outputs_participate_in_collision_checks(tmp_path, make_test_raster):
    config = recipe(tmp_path, make_test_raster)
    config["conflict"] = {"plugin": "band_math", "core:run": True,
                          "core:require_outputs": "param:output_images",
                          "param:input_images": "collect:mul", "param:threshold_math": "B1",
                          "param:output_images": [str(tmp_path / "output/masked_1.tif")]}
    with pytest.raises(ValueError, match="collision"):
        Workflow(config).run()


def test_shared_split_option_updates_outputs_and_scenes(tmp_path, make_test_raster):
    config = recipe(tmp_path, make_test_raster)
    config["shared"]["param:seperate_multipolygons"] = True
    config["mask"].pop("param:seperate_multipolygons")
    workflow = Workflow(config)
    assert not (tmp_path / "output/masked_1.tif").exists()
    workflow.run()
    assert {r["id"] for r in workflow.records} == {"masked_1", "masked_2"}
    assert not (tmp_path / "output/masked.tif").exists()


@pytest.mark.parametrize("cleanup", [False, True])
def test_split_components_become_independent_scenes(tmp_path, make_test_raster, cleanup):
    config = recipe(tmp_path, make_test_raster)
    if cleanup:
        config["shared"].update({"core:delete_temp_steps_proactively": True,
                                 "core:delete_temp_dir": True, "core:cleanup_dirs": ["const:temp_dir"]})
        config["mask"]["core:require_outputs"] = False
        config["mask"]["param:output_images"] = [str(tmp_path / "temp/masked.tif")]
        mask = config.pop("mask")
        config["prepare"] = {"plugin": "band_math", "core:run": True,
                             "param:input_images": "collect:mul", "param:threshold_math": "B1",
                             "var:mul": "expr:const.temp_dir & '/prepared.tif'",
                             "param:output_images": "collect:mul"}
        config["mask"] = mask
    config["files"]["var:footprint_metadata"] = {"view_zenith": 12}
    # Real downstream processing must see each numbered file, with distinct names.
    config["adjust"] = {
        "plugin": "band_math", "core:run": True,
        "param:input_images": "collect:current_image_paths",
        "param:threshold_math": "B1 + 2",
        "var:current_image_paths": "expr:const.temp_dir & '/' & var.basename & '_adjusted.tif'",
        "param:output_images": "collect:current_image_paths",
    }
    config["merge"] = {
        "plugin": "merge_rasters", "core:run": True,
        "core:require_outputs": "param:output_image_path",
        "param:input_images": "collect:current_image_paths",
        "param:output_image_path": str(tmp_path / "output/mosaic.tif"),
    }
    for _ in range(2):
        workflow = Workflow(config)
        workflow.run()
        assert {r["id"] for r in workflow.records} == {"masked_1", "masked_2"}
        for record in workflow.records:
            assert record["context"]["var"]["footprint_metadata"] == {"view_zenith": 12}
            assert str(tmp_path / "source.tif") in record["source_paths"]
            assert not any("masked_" in p for p in record["source_paths"])
        with gdal.Open(str(tmp_path / "output/mosaic.tif")) as ds:
            assert (ds.ReadAsArray() == 3).sum() == 128
        if cleanup:
            assert not list((tmp_path / "temp").glob("*.tif"))
            assert (tmp_path / "source.tif").exists()
