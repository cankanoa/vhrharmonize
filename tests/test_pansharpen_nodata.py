"""Real Orthority coverage of direct output NoData and output-grid preservation."""

import importlib

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin


@pytest.fixture
def inputs(tmp_path):
    pytest.importorskip("orthority")
    rng = np.random.default_rng(17)
    ms = rng.integers(100, 1000, (3, 16, 16)).astype("float32")
    pan = np.repeat(np.repeat(ms.mean(axis=0), 2, axis=0), 2, axis=1)
    pan += rng.uniform(-10, 10, pan.shape)
    return ms, pan


def write_inputs(tmp_path, inputs, dtype, *, cropped=False):
    ms, pan = (data.copy().astype(dtype) for data in inputs)
    ms_nodata = float("nan") if dtype == "float32" else (-32768 if dtype == "int16" else 0)
    ms[:, :, :2] = ms_nodata
    pan[:, :4] = 0
    pan[:, -1] = 0  # Pan-only invalid region must also remain invalid.
    paths = tmp_path / "ms.tif", tmp_path / "pan.tif"
    if cropped:
        ms = ms[:, 2:-2, 2:-2]
    for path, data, transform, nodata in (
        (paths[0], ms, from_origin(2 if cropped else 0, 14 if cropped else 16, 1, 1), ms_nodata),
        (paths[1], pan[None], from_origin(0, 16, .5, .5), 0),
    ):
        with rasterio.open(path, "w", driver="GTiff", count=len(data), dtype=dtype,
                           height=data.shape[1], width=data.shape[2], transform=transform,
                           crs="EPSG:32605", nodata=nodata) as dst:
            dst.write(data)
    return paths


@pytest.mark.parametrize("dtype,nodata", [
    ("int16", -9999), ("int16", 0), ("uint16", 65535),
    ("float32", -9999), ("float32", float("nan")), ("int16", None),
])
def test_real_pansharpen_writes_valid_nodata_without_reopening_output(tmp_path, monkeypatch, inputs, dtype, nodata):
    module = importlib.import_module("vhrharmonize.plugins.pansharpen")
    ms, pan = write_inputs(tmp_path, inputs, dtype)
    baseline, output = tmp_path / "baseline.tif", tmp_path / "output.tif"
    module.pansharpen_image(str(ms), str(pan), str(baseline))
    original_open = rasterio.open
    output_opens = []

    def record_open(path, mode="r", *args, **kwargs):
        if str(path) == str(output):
            output_opens.append(mode)
        return original_open(path, mode, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(rasterio, "open", record_open)
        module.pansharpen_image(str(ms), str(pan), str(output), custom_nodata_value=nodata)
    assert output_opens == ["w"]
    with rasterio.open(baseline) as reference, rasterio.open(output) as result:
        expected_nodata = reference.nodata if nodata is None else nodata
        assert np.isnan(result.nodata) if np.isnan(expected_nodata) else result.nodata == expected_nodata
        np.testing.assert_array_equal(result.read_masks(), reference.read_masks())
        valid = result.read_masks() > 0
        assert valid.any() and not valid.all()
        assert not valid[:, :, :4].any() and not valid[:, :, -1].any()
        np.testing.assert_allclose(result.read()[valid], reference.read()[valid])
        invalid = result.read()[~valid]
        assert np.isnan(invalid).all() if np.isnan(expected_nodata) else (invalid == expected_nodata).all()
        for key in ("transform", "crs", "height", "width", "count", "dtype", "compress", "tiled"):
            assert result.profile[key] == reference.profile[key]
        assert result.colorinterp == reference.colorinterp


def test_custom_nodata_preserves_orthority_cropped_output_grid(tmp_path, inputs):
    module = importlib.import_module("vhrharmonize.plugins.pansharpen")
    ms, pan = write_inputs(tmp_path, inputs, "int16", cropped=True)
    baseline, output = tmp_path / "baseline.tif", tmp_path / "output.tif"
    module.pansharpen_image(str(ms), str(pan), str(baseline))
    module.pansharpen_image(str(ms), str(pan), str(output), custom_nodata_value=-9999)
    with rasterio.open(baseline) as reference, rasterio.open(output) as result:
        assert result.shape == reference.shape == (24, 24)
        assert result.transform == reference.transform
        np.testing.assert_array_equal(result.read(), reference.read())


@pytest.mark.parametrize("dtype,nodata", [
    ("int16", -99999), ("uint16", -1), ("int16", .5),
    ("int16", float("nan")), ("float32", float("inf")), ("float32", 1e40),
])
def test_unrepresentable_custom_nodata_is_rejected_before_creating_output(tmp_path, inputs, dtype, nodata):
    module = importlib.import_module("vhrharmonize.plugins.pansharpen")
    ms, pan = write_inputs(tmp_path, inputs, dtype)
    output = tmp_path / "output.tif"
    with pytest.raises(ValueError, match="not representable"):
        module.pansharpen_image(str(ms), str(pan), str(output), custom_nodata_value=nodata)
    assert not output.exists()
