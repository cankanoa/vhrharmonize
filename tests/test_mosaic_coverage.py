"""Ranked initial coverage limits raster work without changing raw footprint defaults."""

from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
from osgeo import gdal, osr
from pyproj import Transformer
from shapely.geometry import box, mapping
from shapely.ops import transform

from vhrharmonize.plugins.create_footprints import create_footprints
from vhrharmonize.workflow.config import load_config
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.registry import load_plugin, plugin_names
from workflow_helpers import install_function


def test_removed_seamline_metadata_plugin_has_no_alternative_implementation():
    assert "seamline_metadata" not in plugin_names()
    with pytest.raises(ValueError, match="Unknown workflow plugin"):
        load_plugin("seamline_metadata")


def test_worldview_uses_metadata_coverage_then_crops_and_keeps_raw_footprints():
    recipe = load_config("configs/example.worldview.yml")
    names = list(recipe)
    assert [name for name, step in recipe.items() if step.get("plugin") == "create_footprints"] == [
        "standardize_mosaic_coverage", "create_footprints"
    ]
    assert names.index("standardize_mosaic_coverage") < names.index("mask_mosaic_coverage") < names.index("joint_coregistration")
    early, later = recipe["standardize_mosaic_coverage"], recipe["create_footprints"]
    assert early["param:footprint_source"] == "metadata"
    assert early["param:extrude_distance"] > 0
    assert later.get("param:footprint_source", "calculate_bounds") == "calculate_bounds"
    assert later.get("param:image_rank_function") is None
    assert later.get("param:outline_path") is None
    assert later.get("param:extrude_distance", 0) == 0
    assert recipe["mask_mosaic_coverage"]["param:vector_mask"] == [
        "include", "const:coverage_mask", "const:spectralmatch.image_field_name"
    ]


@pytest.mark.parametrize("views,suns,owner,score", [
    ([10, 20, 30], [50, 10, 30], "b", -30),
    ([15, 15, 15], [50, 10, 30], "b", -25),
    ([10], [30], "a", -40),
    ([15, 15, 15], [30, 30, 30], "a", -45),
    ([0, 10], [0, 10], "a", 0),
    ([0, 1, 2], [50, 49, 0], "c", -2),  # Standardizing would prefer a.
])
def test_worldview_prefers_lowest_combined_raw_angle_and_ignores_cloud_cover(tmp_path, views, suns, owner, score):
    recipe = load_config("configs/example.worldview.yml")
    expression = recipe["standardize_mosaic_coverage"]["param:image_rank_function"]
    assert "cloud_cover" not in expression
    records = [{"view_zenith": view, "solar_zenith": sun, "cloud_cover": 99 if index == 1 else 0,
                "geometry": mapping(box(0, 0, 10, 10))}
               for index, (view, sun) in enumerate(zip(views, suns))]
    output = create_footprints(
        [f"{name}.tif" for name in "abc"[:len(records)]], str(tmp_path / "angles.gpkg"),
        footprint_source="metadata", metadata_records=records, metadata_crs="EPSG:3857",
        image_rank_function=expression, image_rank_descending=True,
    )
    result = gpd.read_file(output)
    assert list(result.image) == [owner]
    assert result.weighted_score.iloc[0] == pytest.approx(score)


@pytest.mark.parametrize("prepare_images", [False, True])
def test_coverage_does_not_execute_or_filter_scenes_during_planning(tmp_path, monkeypatch, prepare_images):
    scene_records = []
    to_wgs84 = Transformer.from_crs(3857, 4326, always_xy=True).transform
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3857)
    for name, left, quality in (("partial", 5, 5), ("hidden", 0, 0), ("best", 0, 10), ("outside", 30, 20)):
        source = str(tmp_path / f"{name}.tif")
        dataset = gdal.GetDriverByName("GTiff").Create(source, 10, 10, 1, gdal.GDT_Int16)
        dataset.SetGeoTransform((left, 1, 0, 10, 0, -1))
        dataset.SetProjection(srs.ExportToWkt())
        dataset.GetRasterBand(1).SetNoDataValue(-9999)
        dataset.GetRasterBand(1).Fill(quality + 1)
        dataset = None
        scene_records.append({"scene_id": name, "basename": name, "suffix": "", "current_image_paths": source,
                              "kept_metadata": {"nested": name}, "footprint_metadata": {
                                  "scene_id": name, "quality": quality,
                                  "geometry": mapping(transform(to_wgs84, box(left, 0, left + 10, 10)))}})
    install_function(monkeypatch, "coverage_inputs", lambda: scene_records, scope="aggregate",
                     var_records_return="$", var_id_return="scene_id")
    outline = str(tmp_path / "outline.gpkg")
    gpd.GeoDataFrame(geometry=[box(2, 2, 13, 8)], crs=3857).to_file(outline)
    recipe = load_config("configs/example.worldview.yml")
    seen = []

    def inspect_images(input_images):
        seen.extend(input_images)
        assert len(input_images) == 2
        for path in input_images:
            with gdal.Open(path) as dataset:
                assert dataset.RasterXSize < 10
                assert dataset.RasterYSize < 10
                assert np.any(dataset.GetRasterBand(1).ReadAsArray() != -9999)

    install_function(monkeypatch, "inspect_coverage", inspect_images, scope="aggregate", input_paths={"input_images"})
    config = {
        "shared": {"plugin": "shared", "core:run": True, "core:log_to_console": False,
                   "core:delete_temp_steps_proactively": False, "param:epsg": 3857,
                   "param:custom_nodata_value": -9999},
        "input": {"plugin": "coverage_inputs", "core:run": True, "const:temp_dir": str(tmp_path / "temp")},
        **({"prepare": {
            "plugin": "file_source", "core:run": True,
            "param:input_path": "var:current_image_paths",
            "var:suffix": "_prepared",
            "var:current_image_paths": "expr:const.temp_dir & '/prepare/' & var.basename & var.suffix & '.tif'",
            "param:output_path": "var:current_image_paths",
        }} if prepare_images else {}),
        "setup": recipe["setup_spectralmatch"],
        "coverage": {**recipe["standardize_mosaic_coverage"], "param:outline_path": outline,
                     "param:image_rank_function": "{quality}", "param:extrude_distance": 1},
        "mask": recipe["mask_mosaic_coverage"],
        "footprints": recipe["create_footprints"],
        "inspect": {"plugin": "inspect_coverage", "core:run": True, "core:require_outputs": True,
                    "core:requires": "const:footprints",
                    "param:input_images": "collect:current_image_paths"},
    }
    workflow = Workflow(config)
    counts = workflow.counts()
    assert [record["id"] for record in workflow.records] == ["partial", "hidden", "best", "outside"]
    assert counts["coverage"]["processing"] == 1
    assert seen == []
    assert not (tmp_path / "temp/standardize_mosaic_coverage/coverage.gpkg").exists()
    assert not (tmp_path / "temp/mask_mosaic_coverage").exists()

    result = workflow.run()
    assert [record["id"] for record in result] == ["partial", "best"]
    suffix = "_prepared_coverage" if prepare_images else "_coverage"
    assert seen == [str(tmp_path / "temp/mask_mosaic_coverage" / f"{name}{suffix}.tif")
                    for name in ("partial", "best")]
    assert [record["context"]["var"]["kept_metadata"]["nested"] for record in result] == ["partial", "best"]
    assert not (tmp_path / "temp/mask_mosaic_coverage/hidden_coverage.tif").exists()
    assert workflow.errors == []


@pytest.mark.parametrize("minimum, expected", [(1.0, 1), (None, 2)])
def test_vhr_passes_min_polygon_area_to_spectralmatch(tmp_path, minimum, expected):
    output = create_footprints(
        ["large.tif", "small.tif"], str(tmp_path / "filtered.gpkg"),
        footprint_source="metadata", metadata_crs="EPSG:3857", epsg=3857,
        metadata_records=[
            {"quality": 2, "geometry": mapping(box(0, 0, 10, 10))},
            {"quality": 1, "geometry": mapping(box(9, 0, 10.01, 0.01))},
        ],
        image_rank_function="{quality}", extrude_distance=500,
        min_polygon_area=minimum,
    )
    assert len(gpd.read_file(output)) == expected
