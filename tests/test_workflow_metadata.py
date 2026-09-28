from pathlib import Path
import json
import pytest
from vhrharmonize.io.metadata import read_metadata, parse_imd_text
from vhrharmonize.workflow.config import load_config
from vhrharmonize.workflow.engine import Workflow


def test_imd_preserves_nested_groups_arrays_and_values(tmp_path):
    filename = tmp_path / "source.IMD"
    filename.write_text(
        'version = "1";\nBEGIN_GROUP = IMAGE\nsun = 70.5;\nvalues = (1, 2, (3, 4));\ntime = 2020-01-01T12:00:00Z;\nEND_GROUP = IMAGE\nEND;\n'
    )
    assert read_metadata(str(filename)) == {
        "version": "1",
        "IMAGE": {"sun": 70.5, "values": [1, 2, [3, 4]], "time": "2020-01-01T12:00:00Z"},
    }


@pytest.mark.parametrize(
    "extension, contents, expected",
    [
        (".json", '{"properties": {"angle": 42}}', {"properties": {"angle": 42}}),
        (
            ".xml",
            "<scene><angle>42</angle><band>1</band><band>2</band></scene>",
            {"scene": {"angle": 42, "band": [1, 2]}},
        ),
        (".yaml", "properties:\n  angle: 42\n", {"properties": {"angle": 42}}),
    ],
)
def test_metadata_format_conversion(tmp_path, extension, contents, expected):
    filename = tmp_path / ("source" + extension)
    filename.write_text(contents)
    assert read_metadata(str(filename)) == expected


def test_worldview_recipe_mappings_and_cached_cloudmask_plan(
    make_worldview_bundle, make_test_raster, tmp_path
):
    bundle = make_worldview_bundle()
    filename = bundle["mul_imd"]
    text = filename.read_text().replace("cloudCover = 0.814;", "cloudCover = 0.1;")
    extra = "".join(
        (
            f"BEGIN_GROUP = {band}\nabsCalFactor = 1;\neffectiveBandwidth = 2;\nEND_GROUP = {band}\n"
            for band in ["BAND_G", "BAND_Y", "BAND_R", "BAND_RE", "BAND_N", "BAND_N2"]
        )
    )
    filename.write_text(text.replace("END;", extra + "END;"))
    config = load_config("configs/example.worldview.yml")
    config["import_files"].update(
        {"param:output_dir": str(tmp_path / "output"), "const:temp_dir": str(tmp_path / "temp")}
    )
    config["shared"]["core:validity_check_grid_size"] = 0
    config["import_files"]["param:search_glob"] = str(bundle["mul_tif"])
    for suffix in ["cloudmasked", "cloudmasked_mask"]:
        make_test_raster(
            tmp_path / "output" / f"{bundle['basename']}_py6s_ortho_pansharpen_{suffix}.tif"
        )
    workflow = Workflow(config)
    metadata = workflow.records[0]["context"]["var"]
    constants = workflow.records[0]["context"]["const"]
    assert constants["band_wavelengths_um"] == config["import_files"]["const:band_wavelengths_um"]
    assert constants["calibration"]["WV03"]["BAND_C"] == [0.905, -8.604]
    assert (
        not {"band_order", "band_wavelengths_um", "calibration", "temp_dir"}
        & metadata.keys()
    )
    assert metadata["output_dir"] == str(tmp_path / "output")
    assert workflow.records[0]["id"] == metadata["scene_id"] == bundle["basename"]
    assert metadata["solar_zenith"] == 17
    assert metadata["cloud_cover"] == 10
    assert metadata["day"] == 5
    assert metadata["month"] == 7
    assert metadata["dn_to_radiance_factors"][:2] == pytest.approx([0.4525, 0.47])
    assert workflow.records[0]["context"]["var"]["pan"] == str(bundle["pan_tif"])
    counts = workflow.counts()
    for step in [
        "fetch_atmosphere",
        "atmospheric_correction",
        "orthorectify_mul",
        "orthorectify_pan",
        "pansharpen",
    ]:
        assert counts[step]["processing"] == 0
    assert counts["cloud_mask"]["loaded"] == 1
    assert counts["alignment"]["processing"] == 1
    assert counts["seamline_metadata"]["processing"] == 1


def test_planet_recipe_uses_same_runner(make_test_raster, tmp_path):
    source = make_test_raster(tmp_path / "planet_AnalyticMS_SR.tif", count=4)
    (tmp_path / "planet_metadata.json").write_text(
        json.dumps(
            {
                "id": "planet-scene",
                "properties": {"acquired": "2020-01-02T00:00:00Z", "cloud_cover": 0.2},
            }
        )
    )
    config = load_config("configs/example.planet.yml")
    config["import_files"]["param:search_glob"] = str(source)
    workflow = Workflow(config)
    assert workflow.records[0]["id"] == "planet-scene"
    assert workflow.records[0]["context"]["var"]["scene_id"] == "planet-scene"
    assert workflow.records[0]["context"]["var"]["cloud_cover"] == 20
    assert workflow.counts()["cloud_mask"]["processing"] == 1


def test_seamline_aggregate_receives_projected_metadata_and_resumes(tmp_path, make_test_raster):
    import geopandas as gpd

    for name in ["a", "b"]:
        make_test_raster(tmp_path / (name + ".tif"))
    config = {
        "shared": {"plugin": "shared", "core:run": True, "core:log_to_console": False},
        "import_files": {
            "plugin": "import_files",
            "core:run": True,
            "param:search_glob": str(tmp_path / "*.tif"),
            "var:image": "returned:file_path",
            "var:input": "returned:file_path",
            "var:input_name": "expr:$split(var.file_path, '/')[-1]",
            "var:input_stem": r"expr:$replace($split(var.file_path, '/')[-1], /\.[^.]*$/, '')",
            "var:input_dir": r"expr:$replace(var.file_path, /[^\/]+$/, '')",
            "const:output_dir": str(tmp_path / "output"),
            "const:temp_dir": str(tmp_path / "temp"),
            "var:fields": {"scene_id": "var:input_stem", "quality": 10},
        },
        "seamline_metadata": {
            "plugin": "seamline_metadata",
            "core:run": True,
            "param:image_paths": "collect:image",
            "param:metadata_records": "collect:fields",
            "param:epsg": 4326,
            "const:footprints": "expr:const.output_dir & '/footprints.gpkg'",
            "param:output_path": "const:footprints",
        },
    }
    Workflow(config).run()
    Workflow(config).run()
    result = gpd.read_file(tmp_path / "output" / "footprints.gpkg")
    assert len(result) == 2
    assert set(result.scene_id) == {"a", "b"}
    assert set(result.quality) == {10}


def test_imd_multiple_assignments_and_quoted_semicolons():
    assert parse_imd_text('BEGIN_GROUP = A\nx = 1; y = "a;b";\nEND_GROUP = A\nEND;') == {
        "A": {"x": 1, "y": "a;b"}
    }
    with pytest.raises(ValueError, match="Unclosed"):
        parse_imd_text("BEGIN_GROUP = A\nx = 1;")


def test_footprints_reproject_to_output_crs(tmp_path, make_test_raster):
    from vhrharmonize.plugins.seamline_metadata import _calculate_seamline_metadata_geometry
    from pyproj import Transformer

    source = make_test_raster(tmp_path / "image.tif", crs="EPSG:4326")
    _, geometry = _calculate_seamline_metadata_geometry(
        str(source), {}, "calculate_bounds", True, 3857, "image", "unused.gpkg", False
    )
    x, y = Transformer.from_crs(4326, 3857, always_xy=True).transform(4, 4)
    assert geometry.bounds == pytest.approx((0, 0, x, y))
