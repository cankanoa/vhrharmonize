"""Planning validation must be bounded, visible, and scoped to a single pass."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from threading import Barrier, Lock, get_ident
from types import SimpleNamespace

import pytest

from vhrharmonize.io import validation
from vhrharmonize.workflow.config import DEFAULT_SHARED
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.validation import ValidationPass
from vhrharmonize.workflow import validation as planning
from workflow_helpers import install_function
from test_validation_policy import SampleBand


def test_planning_checks_each_output_once_and_discards_cache(tmp_path, monkeypatch):
    output = tmp_path / 'output.json'
    output.write_text('{}')
    install_function(monkeypatch, 'cached_writer', lambda output_path: output_path,
                     output_paths=('output_path',))
    calls = []
    original = validation._json_file_is_valid
    def check(path):
        calls.append(path)
        return original(path)
    monkeypatch.setattr(validation, '_json_file_is_valid', check)
    workflow = Workflow({
        'settings': {'plugin': 'shared', 'core:run': True,
                     'core:save_statistics_path': None, 'core:load_statistics_path': None},
        'write': {'plugin': 'cached_writer', 'core:run': True, 'core:require_outputs': True,
                  'param:output_path': str(output)},
    })
    workflow.plan()
    assert calls == [str(output)]
    assert workflow.nodes[0].status == 'loaded'
    assert workflow._planning_validation is None
    output.write_text('{broken')
    assert workflow._failures(workflow.nodes[0], [str(output)])
    assert calls == [str(output), str(output)]


def test_file_concurrency_limits_progress_and_tile_setting(tmp_path, monkeypatch, capsys):
    files = [tmp_path / f'{i}.json' for i in range(4)]
    for file in files:
        file.write_text('{}')
    active = 0
    peak = 0
    lock, barrier = Lock(), Barrier(2)
    original = planning._existing_output_failures
    def check(paths, **options):
        nonlocal active, peak
        assert options['tile_threads'] == 3
        with lock:
            active += 1
            peak = max(peak, active)
        try:
            barrier.wait(timeout=10)
            return original(paths, **options)
        finally:
            with lock:
                active -= 1
    @contextmanager
    def executor(workers, backend, scheduler):
        assert workers == 2 and backend == 'process_pool'
        with ThreadPoolExecutor(workers) as pool:
            yield pool
    monkeypatch.setattr(planning, 'image_executor', executor)
    monkeypatch.setattr(planning, '_existing_output_failures', check)
    batch = ValidationPass({**DEFAULT_SHARED, 'concurrent_processing': 2, 'tile_threads': 3})
    assert batch.failures([str(f) for f in files] * 2, True) == {}
    assert peak == 2
    lines = [line for line in capsys.readouterr().out.splitlines() if 'core:validation' in line]
    assert len(lines) == 5
    assert '0%' in lines[0] and '25%' in lines[1] and '100%' in lines[-1]
    assert batch.failures([str(f) for f in files], True) == {}
    assert capsys.readouterr().out == ''


def test_validation_uses_real_process_pool(tmp_path):
    good, bad = tmp_path / 'good.json', tmp_path / 'bad.json'
    good.write_text('{}')
    bad.write_text('{broken')
    batch = ValidationPass({**DEFAULT_SHARED, 'concurrent_processing': 2, 'log_to_console': False})
    assert set(batch.failures([str(good), str(bad)], True)) == {str(bad)}


def test_disabled_validation_never_opens_files_or_starts_workers(tmp_path, monkeypatch):
    file = tmp_path / 'broken.tif'
    file.write_bytes(b'incomplete')
    def unexpected(*args, **kwargs):
        pytest.fail('Disabled validation must only check existence')
    monkeypatch.setattr(planning, 'image_executor', unexpected)
    monkeypatch.setattr(planning, '_existing_output_failures', unexpected)
    batch = ValidationPass({**DEFAULT_SHARED, 'concurrent_processing': 4})
    assert batch.failures([str(file)], False) == {}
    assert batch.failures([str(tmp_path / 'missing.tif')], False)


@pytest.mark.parametrize('spacing', [1, 2])
def test_tile_threads_read_same_grid_with_private_handles(monkeypatch, spacing):
    reads, owners = [], []
    barrier = Barrier(2)
    class Band(SampleBand):
        def __init__(self, owner, index):
            super().__init__()
            self.owner, self.index = owner, index
        def ReadRaster(self, x, y, w, h, **kwargs):
            assert self.owner == get_ident()
            reads.append((self.index, x, y, w, h))
            return super().ReadRaster(x, y, w, h, **kwargs)
    def open_dataset(*args):
        owner = get_ident()
        owners.append(owner)
        barrier.wait(timeout=10)
        return SimpleNamespace(RasterCount=2, GetRasterBand=lambda i: Band(owner, i))
    monkeypatch.setattr(validation.gdal, 'OpenEx', open_dataset)
    validation._read_raster_grid_parallel('unused', SimpleNamespace(RasterCount=2),
                                          9, 7, spacing, 'readable', 2)
    assert len(set(owners)) == 2
    stride = 512 if spacing == 1 else spacing
    expected = [(i, x, y, min(512, 9-x) if spacing == 1 else 1,
                 min(512, 7-y) if spacing == 1 else 1)
                for i in (1, 2) for y in range(0, 7, stride) for x in range(0, 9, stride)]
    assert sorted(reads) == sorted(expected)


@pytest.mark.parametrize('stop', ['readable', 'has_valid_data'])
def test_parallel_real_raster_validation(tmp_path, make_test_raster, stop):
    file = str(make_test_raster(tmp_path / 'raster.tif', width=64, height=64))
    assert validation._gdal_raster_is_valid(file, validity_raster_grid=2,
        validity_raster_check_stop=stop, tile_threads=3) == (True, None)


def test_parallel_raster_read_error_is_not_ignored(tmp_path, make_test_raster, monkeypatch):
    file = str(make_test_raster(tmp_path / 'raster.tif', width=64, height=64))
    original = validation._gdal_band_window_has_valid_pixel
    def check(band, x, y, *args, **kwargs):
        if (x, y) == (62, 62):
            raise IOError('bad tile at end of raster')
        return original(band, x, y, *args, **kwargs)
    monkeypatch.setattr(validation, '_gdal_band_window_has_valid_pixel', check)
    valid, reason = validation._gdal_raster_is_valid(file, validity_raster_grid=2, tile_threads=3)
    assert not valid and 'bad tile at end' in reason
