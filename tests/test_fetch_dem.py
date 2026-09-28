"""DEM coverage and workflow aggregation, without external API requests."""

from pathlib import Path
from types import SimpleNamespace

import pytest
import rasterio
from rasterio.control import GroundControlPoint
from rasterio.rpc import RPC
from rasterio.transform import from_bounds
from rasterio.warp import transform_bounds
import yaml

from vhrharmonize.plugins import fetch_dem
from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import import_settings


@pytest.fixture
def dem_request(monkeypatch, make_test_raster, tmp_path):
    content = make_test_raster(tmp_path / "response/dem.tif").read_bytes()
    calls = []

    def get(url, **kwargs):
        calls.append(kwargs["params"])
        return SimpleNamespace(
            headers={"content-type": "image/tiff"}, content=content,
            raise_for_status=lambda: None,
        )

    monkeypatch.setattr(fetch_dem.requests, "get", get)
    monkeypatch.setenv("OPENTOPOGRAPHY_API_KEY", "test-key")
    return calls


def test_union_of_images_in_different_crss(dem_request, make_test_raster, tmp_path):
    first = make_test_raster(tmp_path / "first.tif")  # 0, 0, 4, 4
    second = make_test_raster(tmp_path / "second.tif", crs="EPSG:3857")
    projected = transform_bounds("EPSG:4326", "EPSG:3857", 5, -2, 7, 3)
    with rasterio.open(second, "r+") as dst:
        dst.transform = from_bounds(*projected, dst.width, dst.height)
    result = fetch_dem.download_opentopography_dem_for_bbox(
        image_paths=[str(first), str(second)], output_tif_path=str(tmp_path / "dem.tif")
    )
    assert Path(result).is_file()
    assert len(dem_request) == 1
    assert [float(dem_request[0][key]) for key in ("west", "south", "east", "north")] == pytest.approx([0, -2, 7, 4])


def test_single_image(dem_request, make_test_raster, tmp_path):
    image = make_test_raster(tmp_path / "image.tif")
    fetch_dem.download_opentopography_dem_for_bbox(
        image_paths=str(image), output_tif_path=str(tmp_path / "dem.tif")
    )
    assert len(dem_request) == 1
    assert dem_request[0]["east"] == "4.00000000"


@pytest.mark.parametrize("mode", ["rpc", "gcp"])
def test_sensor_georeferencing_without_affine_crs(mode, make_test_raster, tmp_path):
    image = make_test_raster(tmp_path / "image.tif", crs=None)
    with rasterio.open(image, "r+") as dst:
        if mode == "rpc":
            denominator = [1] + [0] * 19
            line, sample = [0] * 20, [0] * 20
            line[2], sample[1] = 1, 1
            rpc = RPC(height_off=0, height_scale=100, lat_off=10, lat_scale=2,
                      long_off=20, long_scale=2, line_off=2, line_scale=2,
                      samp_off=2, samp_scale=2, line_num_coeff=line,
                      samp_num_coeff=sample, line_den_coeff=denominator,
                      samp_den_coeff=denominator)
            dst.update_tags(ns="RPC", **rpc.to_gdal())
        else:
            dst.gcps = ([GroundControlPoint(row=r, col=c, x=18+c, y=8+r)
                         for r, c in [(0, 0), (0, 4), (4, 0), (4, 4)]],
                        rasterio.crs.CRS.from_epsg(4326))
    west, south, east, north = fetch_dem._image_bounds(str(image))
    assert 17 <= west <= 18
    assert 7 <= south <= 8
    assert east - west == pytest.approx(4)
    assert north - south == pytest.approx(4)


@pytest.mark.parametrize("options", [
    {}, {"image_paths": []}, {"min_lon": 0},
    {"image_paths": ["unused.tif"], "min_lon": 0},
    {"min_lon": 2, "min_lat": 0, "max_lon": 1, "max_lat": 1},
])
def test_invalid_inputs_do_not_download(options, dem_request, tmp_path):
    with pytest.raises(ValueError):
        fetch_dem.download_opentopography_dem_for_bbox(
            output_tif_path=str(tmp_path / "dem.tif"), **options
        )
    assert not dem_request


@pytest.mark.parametrize("aggregate", [True, False])
def test_workflow_collects_once_or_fetches_per_scene(aggregate, dem_request, make_test_raster, tmp_path):
    for name in ("first", "second"):
        make_test_raster(tmp_path / f"source/{name}.tif")
    example = yaml.safe_load((Path(__file__).parents[1] / "configs/example.worldview.yml").read_text())
    step = {**example["fetch_dem"], "core:run": True}
    if not aggregate:
        step.pop("const:dem_path")
        step.pop("param:output_tif_path")
        step.update({
            "core:scope": "scene", "param:image_paths": "var:mul",
            "var:dem": "expr:const.temp_dir & '/fetch_dem/' & var.basename & '.tif'",
            "param:output_tif_path": "var:dem",
        })
    config = {
        "shared": {"plugin": "shared", "core:run": True, "core:log_to_console": False},
        "import": {**import_settings(tmp_path / "source/*.tif", tmp_path),
                   "var:current_image_paths": "returned:file_path"},
        "dem": step,
        "consume": {
            "plugin": "file_source", "core:run": True,
            "param:input_path": "const:dem_path" if aggregate else "var:dem",
            "param:output_path": "expr:const.output_dir & '/' & var.basename & '.tif'",
        },
    }
    workflow = Workflow(config)
    assert not dem_request  # Planning does not download.
    workflow.run()
    assert len(dem_request) == (1 if aggregate else 2)
    assert len(list((tmp_path / "output").glob("*.tif"))) == 2
    assert not list((tmp_path / "temp/fetch_dem").glob("*.tif"))  # Default cleanup after all consumers.
    assert len(list((tmp_path / "source").glob("*.tif"))) == 2
