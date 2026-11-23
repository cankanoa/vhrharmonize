"""Compression must survive the single-scene bypass and core finalization."""
import numpy as np
import pytest
import rasterio
from vhrharmonize.workflow.registry import load_plugin


@pytest.mark.parametrize('compression', [True, False])
def test_single_scene_bypass_honors_compression(tmp_path, make_test_raster, compression):
    source = make_test_raster(tmp_path / 'source.tif', width=32, height=32, fill=7)
    destination = tmp_path / 'result.tif'
    plugin = load_plugin('global_regression')
    plugin.run(params={
        'input_images': [str(source)], 'output_images': [str(destination)],
        'skip_single_image': True, 'compression': compression,
        'validity_writing_postfix': '_writing',
    }, shared={})
    with rasterio.open(source) as a, rasterio.open(destination) as b:
        assert (b.compression is not None) == compression
        assert a.transform == b.transform and a.crs == b.crs
        np.testing.assert_array_equal(a.read(), b.read())
    assert not (tmp_path / 'result_writing.tif').exists()
