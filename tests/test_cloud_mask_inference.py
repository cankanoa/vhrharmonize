"""Predict on MUL imagery while preserving the pansharpened application grid."""

import sys
from types import SimpleNamespace

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from vhrharmonize.plugins.cloud_mask import (
    apply_binary_cloud_mask_to_image,
    cloudmask_raster,
)
from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import copy_step, import_settings


def test_separate_inference_raster_preserves_target_bands_and_grid(
    monkeypatch, tmp_path, make_test_raster,
):
    mul = make_test_raster(tmp_path / "mul.tif", count=3, width=2, height=2, fill=17)
    with rasterio.open(mul, "r+") as ds:
        ds.transform = from_origin(0, 4, 2, 2)
    pan = make_test_raster(tmp_path / "pansharpened.tif", count=8, fill=42)

    def predict(array, **kwargs):
        np.testing.assert_array_equal(array, np.full((3, 2, 2), 17))
        return np.array([[0, 1], [0, 0]], dtype="uint8")

    monkeypatch.setitem(sys.modules, "omnicloudmask", SimpleNamespace(predict_from_array=predict))
    result = cloudmask_raster(
        str(pan), str(tmp_path / "masked.tif"), str(tmp_path / "mask.tif"),
        inference_image_path=str(mul), red_band_index=1, green_band_index=2,
        nir_band_index=3, inference_resolution_m=None,
    )
    assert result.mask_pixel_count == 1
    with rasterio.open(result.output_mask_path) as mask:
        assert mask.shape == (2, 2)
        assert mask.transform == from_origin(0, 4, 2, 2)
    with rasterio.open(result.output_raster_path) as output, rasterio.open(pan) as source:
        assert output.profile == source.profile
        expected = source.read()
        expected[:, :2, 2:] = 0
        np.testing.assert_array_equal(output.read(), expected)


@pytest.mark.parametrize("allow_reprojection", [True, False])
def test_mask_grid_checked_even_when_dimensions_match(
    tmp_path, make_test_raster, allow_reprojection,
):
    image = make_test_raster(tmp_path / "image.tif", fill=42)
    mask = make_test_raster(tmp_path / "mask.tif", fill=1)
    with rasterio.open(mask, "r+") as ds:
        ds.transform = from_origin(2, 4, 1, 1)
    output = tmp_path / "masked.tif"
    if not allow_reprojection:
        with pytest.raises(ValueError, match="Cloud mask grid"):
            apply_binary_cloud_mask_to_image(
                str(image), str(mask), str(output), allow_mask_reprojection=False,
            )
        assert not output.exists()
    else:
        apply_binary_cloud_mask_to_image(str(image), str(mask), str(output))
        with rasterio.open(output) as ds:
            expected = np.full((1, 4, 4), 42)
            expected[:, :, 2:] = 0
            np.testing.assert_array_equal(ds.read(), expected)


def test_planner_rebuilds_inference_input_even_when_pansharpened_image_is_reused(
    monkeypatch, tmp_path, make_test_raster,
):
    source = make_test_raster(tmp_path / "raw.tif", count=3, fill=17)
    pan = make_test_raster(tmp_path / "pansharpened.tif", count=8, fill=42)
    mul = tmp_path / "mul.tif"
    output = tmp_path / "masked.tif"
    inference_shapes = []

    def predict(array, **kwargs):
        inference_shapes.append(array.shape)
        assert (array == 17).all()
        mask = np.zeros(array.shape[1:], dtype="uint8")
        mask[0, 0] = 1
        return mask

    monkeypatch.setitem(sys.modules, "omnicloudmask", SimpleNamespace(predict_from_array=predict))
    recipe = {
        "shared": {
            "plugin": "shared", "core:run": True, "core:log_to_console": False,
            "core:show_progress": False, "core:save_statistics_path": None,
            "core:load_statistics_path": None,
        },
        "files": import_settings(source, tmp_path),
        "ortho": copy_step("current_image_paths", "mul", str(mul)),
        "pan": copy_step("current_image_paths", "current_image_paths", str(pan)),
        "cloud": {
            "plugin": "cloud_mask", "core:run": True,
            "core:require_outputs": ["param:output_raster_path", "param:output_mask_path"],
            "param:input_image_path": "var:current_image_paths",
            "param:inference_image_path": "var:orthorectified_mul_path",
            "param:output_raster_path": str(output),
            "param:output_mask_path": str(tmp_path / "mask.tif"),
            "param:red_band_index": 1, "param:green_band_index": 2,
            "param:nir_band_index": 3, "param:inference_resolution_m": None,
        },
    }
    recipe["ortho"]["var:orthorectified_mul_path"] = "var:current_image_paths"
    workflow = Workflow(recipe)
    assert workflow.counts()["ortho"]["processing"] == 1
    assert workflow.counts()["pan"]["loaded"] == 1
    workflow.run()
    assert inference_shapes == [(3, 4, 4)]
    with rasterio.open(output) as ds:
        assert ds.count == 8
        expected = np.full((8, 4, 4), 42)
        expected[:, 0, 0] = 0
        np.testing.assert_array_equal(ds.read(), expected)
    with rasterio.open(pan) as ds:
        assert (ds.read() == 42).all()
