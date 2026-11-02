"""Exercise real GDAL RPC warping while keeping original files immutable."""

import hashlib
import importlib
from types import SimpleNamespace

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
from osgeo import gdal

ortho = importlib.import_module('vhrharmonize.plugins.orthorectification')


@pytest.mark.parametrize('refined', [False, True])
def test_direct_rpc_warp_preserves_source_and_applies_output_options(tmp_path, monkeypatch, refined):
    source, dem, output = (tmp_path / name for name in ('source.tif', 'dem.tif', 'ortho.tif'))
    data = np.arange(32 * 32, dtype='int16').reshape(32, 32) + 100
    data[:3] = -32768
    with rasterio.open(source, 'w', driver='GTiff', count=1, width=32, height=32,
                       dtype='int16', nodata=-32768) as dst:
        dst.write(data, 1)
        den = [1.] + [0.] * 19
        line, samp = [0.] * 20, [0.] * 20
        line[2], samp[1] = -1., 1.
        rpc = dict(LINE_OFF='15.5', LINE_SCALE='15.5', SAMP_OFF='15.5', SAMP_SCALE='15.5',
                   LAT_OFF='20', LAT_SCALE='.001', LONG_OFF='-155', LONG_SCALE='.001',
                   HEIGHT_OFF='0', HEIGHT_SCALE='1000',
                   LINE_NUM_COEFF=' '.join(map(str, line)), LINE_DEN_COEFF=' '.join(map(str, den)),
                   SAMP_NUM_COEFF=' '.join(map(str, samp)), SAMP_DEN_COEFF=' '.join(map(str, den)))
        dst.update_tags(ns='RPC', **rpc)
    with rasterio.open(dem, 'w', driver='GTiff', count=1, width=100, height=100,
                       dtype='float32', crs='EPSG:4326',
                       transform=from_origin(-156, 21, .02, .02)) as dst:
        dst.write(np.zeros((100, 100), dtype='float32'), 1)
    digest = hashlib.sha256(source.read_bytes()).digest()
    translate = gdal.Translate
    formats = []
    def translated(*args, **kwargs):
        formats.append(kwargs.get('format'))
        assert refined and kwargs.get('format') == 'VRT'
        return translate(*args, **kwargs)
    monkeypatch.setattr(gdal, 'Translate', translated)
    if refined:
        import orthority
        with rasterio.open(source) as src:
            camera = src.rpcs
        monkeypatch.setattr(orthority, 'RpcCameras', SimpleNamespace(
            from_images=lambda paths: SimpleNamespace(
                refine=lambda path: None, get=lambda path: SimpleNamespace(_rpc=camera))))
    ortho.gcp_refined_rpc_orthorectification(
        str(source), str(output), str(dem), 4326,
        gcp_geojson_file_path='gcps.geojson' if refined else None,
        custom_nodata_value=-9999, output_dtype='Float32', output_resolution=.00005,
    )
    assert formats == (['VRT'] if refined else [])
    assert hashlib.sha256(source.read_bytes()).digest() == digest
    with rasterio.open(output) as result:
        assert result.dtypes == ('float32',)
        assert result.nodata == -9999
        assert result.crs.to_epsg() == 4326
        assert result.res == (.00005, .00005)
        values = result.read(1, masked=True)
        assert values.mask.any() and not values.mask.all()
        assert values.min() >= 100
