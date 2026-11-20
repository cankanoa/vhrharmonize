"""Seamline samples retain their values without whole-source warp allocations."""

from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from osgeo import gdal, osr
import pytest

from vhrharmonize.plugins import _markov_sampling as sampling
from vhrharmonize.plugins.markov_triangles import markov_triangles


class _ForbiddenSampler:
    def __init__(self, crs):
        raise AssertionError('Worker used the unbounded native sampler')


def forbid_native_sampler():
    sampling._native._RasterSampler = _ForbiddenSampler


def raster(path, *, width=520, height=280, transform=(0, 1, 0, 280, 0, -1), sparse=False):
    crs = osr.SpatialReference()
    crs.ImportFromEPSG(32605)
    with gdal.GetDriverByName('GTiff').Create(
        str(path), width, height, 8 if sparse else 2, gdal.GDT_Int16,
        options=['TILED=YES', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256', 'SPARSE_OK=YES', 'BIGTIFF=YES'],
    ) as ds:
        ds.SetProjection(crs.ExportToWkt())
        ds.SetGeoTransform(transform)
        rows, cols = np.indices((256, 256) if sparse else (height, width))
        for number in range(1, ds.RasterCount + 1):
            band = ds.GetRasterBand(number)
            band.SetNoDataValue(-9999)
            values = (rows + cols + number * 50).astype(np.int16)
            if not sparse:
                values[:20, :20] = -9999
                if number == 2:
                    values[30:40, 30:40] = -9999
            band.WriteArray(values)
    return str(path), crs.ExportToWkt()


@pytest.mark.parametrize('gsd', [None, 1, 2, 100])
@pytest.mark.parametrize('rotated', [False, True])
def test_samples_match_native_reader_with_masks_edges_and_band_order(tmp_path, gsd, rotated):
    transform = (0.3, 1, 0.2 if rotated else 0, 280.4, 0.1 if rotated else 0, -1)
    path, crs = raster(tmp_path / 'source.tif', transform=transform)
    pixels = [(0, 0), (10, 10), (35, 35), (25, 25), (255, 255),
              (256, 256), (510, 270), (520, 280), (-1, -1), (530, 300)]
    coordinates = np.array([gdal.ApplyGeoTransform(transform, x, y) for x, y in pixels])
    native = sampling._native._RasterSampler(crs)
    bounded = sampling.BoundedRasterSampler(crs)
    try:
        for fallback in (False, True):
            expected = native.sample(path, coordinates, (2, 1), gsd, native_fallback=fallback)
            actual = bounded.sample(path, coordinates, (2, 1), gsd, native_fallback=fallback)
            np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=1e-8, equal_nan=True)
        assert bounded.sample(path, np.empty((0, 2)), (1,), gsd).shape == (0, 1)
    finally:
        native.close()
        bounded.close()
    assert not bounded.sources and not bounded.grids and not bounded.tiles


def test_cache_is_bounded_and_repeated_samples_reuse_tiles(tmp_path, monkeypatch):
    path, crs = raster(tmp_path / 'source.tif')
    sampler = sampling.BoundedRasterSampler(crs)
    sampler.cache_bytes = 8  # One coarse sampled pixel, to exercise eviction.
    calls = []
    original = gdal.Warp

    def warp(*args, **kwargs):
        calls.append(kwargs['format'])
        return original(*args, **kwargs)

    monkeypatch.setattr(gdal, 'Warp', warp)
    try:
        for point in ((100, 180), (100, 180), (300, 180), (100, 180)):
            sampler.sample(path, np.array([point]), (1,), 1)
            assert sampler.tile_bytes <= sampler.cache_bytes
        assert calls == ['VRT', 'MEM', 'MEM', 'MEM']
    finally:
        sampler.close()


def test_coarse_reads_use_tight_windows_and_only_extend_uncached_strips(tmp_path, monkeypatch):
    path, crs = raster(tmp_path / 'source.tif')
    sampler = sampling.BoundedRasterSampler(crs)
    calls = []
    original = gdal.Warp

    def warp(*args, **kwargs):
        if kwargs.get('format') == 'MEM':
            calls.append((kwargs['width'], kwargs['height']))
        return original(*args, **kwargs)

    monkeypatch.setattr(gdal, 'Warp', warp)
    try:
        # The same coarse bin grows from a point to a 3x3 window. Repeated
        # requests must not average the original pixel again.
        sampler.sample(path, np.array([[50., 230.]]), (1,), 10)
        assert calls == [(1, 1)]
        coordinates = np.array([[40., 240.], [60., 220.], [50., 230.]])
        actual = sampler.sample(path, coordinates, (1,), 10)
        assert sum(w*h for w,h in calls) == 9
        assert all(w <= 3 and h <= 3 for w,h in calls)
        count = len(calls)
        np.testing.assert_array_equal(sampler.sample(path, coordinates, (1,), 10), actual)
        assert len(calls) == count
        native = sampling._native._RasterSampler(crs)
        try:
            np.testing.assert_allclose(actual, native.sample(path, coordinates, (1,), 10),
                                       rtol=1e-8, atol=1e-8, equal_nan=True)
        finally:
            native.close()
    finally:
        sampler.close()


def test_overrides_restore_after_errors(monkeypatch):
    original_solver = sampling._native._solve_node
    original_sampler = sampling._native._RasterSampler

    def fail():
        assert sampling._native._RasterSampler is sampling.BoundedRasterSampler
        raise RuntimeError('read failed')

    monkeypatch.setattr(sampling, '_native_solve_node', fail)
    with pytest.raises(RuntimeError, match='read failed'):
        with sampling.bounded_markov_sampling():
            assert sampling._native._solve_node is sampling._solve_node
            sampling._native._solve_node()
    assert sampling._native._solve_node is original_solver
    assert sampling._native._RasterSampler is original_sampler


@pytest.mark.parametrize('execution', ['serial', 'native_pool', 'spawn_executor'])
def test_real_solver_preserves_seamlines_in_worker_backends(tmp_path, execution):
    import geopandas as gpd
    from shapely import symmetric_difference, union_all
    from spectralmatch import markov_triangles as native_markov

    a, _ = raster(tmp_path / 'a.tif', width=32, height=32, transform=(0, 1, 0, 32, 0, -1))
    b, _ = raster(tmp_path / 'b.tif', width=32, height=32, transform=(16, 1, 0, 32, 0, -1))
    params = dict(input_images=[a, b], mesh_max_triangle_area=16,
                  mesh_min_triangle_angle=0, image_quality_weight=0)
    expected_path = str(tmp_path / 'expected.gpkg')
    actual_path = str(tmp_path / 'actual.gpkg')
    native_markov(**params, output_mask=expected_path)
    if execution == 'spawn_executor':
        with ProcessPoolExecutor(max_workers=1, mp_context=get_context('spawn'),
                                 initializer=forbid_native_sampler) as executor:
            markov_triangles(**params, output_mask=actual_path, executor=executor)
    else:
        markov_triangles(**params, output_mask=actual_path,
                         image_threads=1 if execution == 'native_pool' else None)
    expected = gpd.read_file(expected_path).set_index('image')
    actual = gpd.read_file(actual_path).set_index('image')
    assert set(actual.index) == set(expected.index)
    for name in expected.index:
        assert symmetric_difference(actual.loc[name].geometry, expected.loc[name].geometry).area < 1e-8
    assert union_all(actual.geometry).area > 0


@pytest.mark.skipif(sys.platform != 'linux', reason='Linux address-space and RSS regression check')
def test_large_sparse_raster_samples_within_memory_limit(tmp_path):
    path, _ = raster(tmp_path / 'large.tif', width=8192, height=8192,
                     transform=(0, 0.5, 0, 4096, 0, -0.5), sparse=True)
    script = '''
import resource
import sys
import os
from pathlib import Path
import numpy as np
from osgeo import gdal
from vhrharmonize.plugins._markov_sampling import BoundedRasterSampler, _native
gdal.UseExceptions()
gdal.SetCacheMax(32 * 1024**2)
with gdal.Open(sys.argv[1]) as ds:
    crs = ds.GetProjection()
# Imported scientific libraries reserve substantial address space. Limit the
# additional allocation instead of assuming a particular import-time footprint.
baseline_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
baseline_vms = int(Path('/proc/self/statm').read_text().split()[0]) * os.sysconf('SC_PAGE_SIZE')
limit = baseline_vms + 1024**3
resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
point = np.array([[50., 4050.]])
native = _native._RasterSampler(crs)
try:
    native.sample(sys.argv[1], point, (1,), 100)
except RuntimeError as exc:
    assert 'allocate' in str(exc), str(exc)
else:
    raise AssertionError('Original sampler no longer reproduces the oversized allocation')
finally:
    native.close()
bounded = BoundedRasterSampler(crs)
try:
    result = bounded.sample(sys.argv[1], point, (1,), 100)
    assert np.isfinite(result).all(), result
finally:
    bounded.close()
added_mib = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss - baseline_rss) / 1024
assert added_mib < 512, added_mib
print('added_mib', added_mib)
'''
    env = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'}
    result = subprocess.run([sys.executable, '-c', script, path], env=env,
                            cwd=Path(__file__).parents[1], capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr


def test_adjacent_coarse_bins_share_one_warp(tmp_path, monkeypatch):
    path, crs = raster(tmp_path / 'adjacent.tif')
    sampler = sampling.BoundedRasterSampler(crs)
    calls = []
    original = gdal.Warp
    def warp(*args, **kwargs):
        if kwargs.get('format') == 'MEM':
            calls.append((kwargs['width'], kwargs['height']))
        return original(*args, **kwargs)
    monkeypatch.setattr(gdal, 'Warp', warp)
    coordinates = np.array([[159., 180.], [160., 180.]])
    native = sampling._native._RasterSampler(crs)
    try:
        actual = sampler.sample(path, coordinates, (1,2), 10)
        assert calls == [(2,1)]
        np.testing.assert_allclose(actual, native.sample(path, coordinates, (1,2), 10), equal_nan=True)
    finally:
        native.close(); sampler.close()


def test_worker_cache_reuses_arrays_and_invalidates_changed_sources(tmp_path, monkeypatch):
    path, crs = raster(tmp_path / 'cache.tif')
    calls = []
    original = sampling.BoundedRasterSampler._read_window
    def read(self, *args, **kwargs):
        calls.append(1)
        return original(self,*args,**kwargs)
    monkeypatch.setattr(sampling.BoundedRasterSampler,'_read_window',read)
    point=np.array([[100.,180.]])
    for _ in range(2):
        sampler=sampling.BoundedRasterSampler(crs)
        sampler.sample(path,point,(1,),10);sampler.close()
    assert len(calls)==1
    with gdal.Open(path,gdal.GA_Update) as ds:
        ds.GetRasterBand(1).WriteArray(np.full((10,10),3000,dtype=np.int16),100,100)
    stamp=os.stat(path).st_mtime_ns
    os.utime(path,ns=(stamp,stamp+1000000000))  # Deterministic on coarse-timestamp filesystems.
    sampler=sampling.BoundedRasterSampler(crs)
    sampler.sample(path,point,(1,),10);sampler.close()
    assert len(calls)==2
    assert sampling._worker_tile_bytes <= sampling._worker_cache_limit
