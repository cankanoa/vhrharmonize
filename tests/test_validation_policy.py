"""Validation timing, complete grid scans, and publication of individual outputs."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Event
from types import SimpleNamespace
import shutil
import struct

import pytest
from osgeo import gdal

from vhrharmonize.io import validation
from vhrharmonize.workflow.config import DEFAULT_SHARED, validate_config
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.output_files import OutputPolicy, OutputExecutor, writing_path
from vhrharmonize.plugins.spectralmatch_base import SpectralMatchPlugin
from workflow_helpers import install_function


class SampleBand:
    def __init__(self, fail=None):
        self.reads = []
        self.fail = fail
    def ReadRaster(self, x, y, w, h, **kwargs):
        self.reads.append((x, y, w, h))
        if self.fail == (x, y):
            return None
        return struct.pack('=d', 1) * (w * h)
    def GetMaskBand(self):
        return None
    def GetNoDataValue(self):
        return 0


def test_grid_readable_checks_all_bands_and_does_not_stop_at_valid_data():
    bands = [SampleBand(), SampleBand(fail=(4, 2))]
    ds = SimpleNamespace(RasterCount=2, GetRasterBand=lambda i: bands[i-1])
    with pytest.raises(IOError, match="Incomplete raster read"):
        validation._read_raster_grid(ds, 5, 3, 2, 'readable')
    assert bands[0].reads == [(x, y, 1, 1) for y in (0, 2) for x in (0, 2, 4)]
    assert bands[1].reads[-1] == (4, 2, 1, 1)


def test_grid_has_valid_data_explicitly_stops_early():
    band = SampleBand(fail=(2, 0))
    ds = SimpleNamespace(RasterCount=1, GetRasterBand=lambda i: band)
    assert validation._read_raster_grid(ds, 5, 3, 2, 'has_valid_data')
    assert band.reads == [(0, 0, 1, 1)]


def test_grid_one_reads_every_pixel_in_bounded_windows():
    band = SampleBand()
    ds = SimpleNamespace(RasterCount=1, GetRasterBand=lambda i: band)
    validation._read_raster_grid(ds, 1030, 517, 1, 'readable')
    assert sum(w * h for x, y, w, h in band.reads) == 1030 * 517
    assert len(band.reads) == 6


def test_grid_false_does_not_read_pixels(tmp_path, make_test_raster, monkeypatch):
    filename = str(make_test_raster(tmp_path / 'input.tif'))
    def unexpected(*args, **kwargs):
        pytest.fail('Pixel reading is disabled')
    monkeypatch.setattr(validation, '_read_raster_grid', unexpected)
    assert validation._gdal_raster_is_valid(filename, validity_raster_grid=False) == (True, None)
    with open(filename, 'r+b') as stream:
        stream.truncate(100)
    assert not validation._gdal_raster_is_valid(filename, validity_raster_grid=False)[0]


@pytest.mark.parametrize('setting,value', [('check_validity', True), ('check_validity', 'always'),
    ('validity_raster_grid', 0), ('validity_raster_grid', True), ('validity_raster_grid', -1),
    ('validity_raster_check_stop', 'first'), ('validity_writing_postfix', '../writing')])
def test_invalid_validation_settings_are_rejected(setting, value):
    with pytest.raises(ValueError):
        validate_config({'settings': {'plugin': 'shared', 'core:run': True, 'core:' + setting: value}})


def test_validation_defaults():
    assert {k: DEFAULT_SHARED[k] for k in ('check_validity', 'validity_raster_grid',
            'validity_raster_check_stop', 'validity_writing_postfix')} == {
        'check_validity': 'before_reuse', 'validity_raster_grid': 255,
        'validity_raster_check_stop': 'readable', 'validity_writing_postfix': '_writing'}


@pytest.mark.parametrize('mode,accepts', [('around', False), ('before_reuse', True), (False, True)])
def test_new_output_validation_happens_before_rename_only_in_around(tmp_path, mode, accepts):
    output = str(tmp_path / 'result.json')
    policy = OutputPolicy(check_validity=mode, validation_paths=(output,))
    def write(path):
        assert path.endswith('_writing.json')
        assert not Path(output).exists()
        Path(path).write_text('{incomplete')
        return path
    if accepts:
        assert policy.run(write, (output,), outputs=(output,)) == output
        assert Path(output).exists()
    else:
        with pytest.raises(RuntimeError, match='valid declared outputs'):
            policy.run(write, (output,), outputs=(output,))
        assert not Path(output).exists()
        assert Path(writing_path(output, '_writing')).exists()


def test_failed_output_never_replaces_final_and_next_attempt_starts_clean(tmp_path):
    output = str(tmp_path / 'image.txt')
    Path(output).write_text('previous complete image')
    policy = OutputPolicy()
    def fail(path):
        Path(path).write_text('partial')
        raise RuntimeError('interrupted')
    with pytest.raises(RuntimeError, match='interrupted'):
        policy.run(fail, (output,), outputs=(output,))
    assert Path(output).read_text() == 'previous complete image'
    def success(path):
        assert not Path(path).exists()
        Path(path).write_text('complete')
        return path
    assert policy.run(success, (output,), outputs=(output,)) == output
    assert Path(output).read_text() == 'complete'
    assert not Path(writing_path(output, '_writing')).exists()


@pytest.mark.parametrize('postfix', [None, ''])
def test_postfix_can_be_disabled(tmp_path, postfix):
    output = str(tmp_path / 'out.txt')
    def write(path):
        assert path == output
        Path(path).write_text('done')
    OutputPolicy(postfix=postfix).run(write, (output,), outputs=(output,))
    assert Path(output).exists()


def test_each_image_publishes_before_other_workers_finish(tmp_path):
    from spectralmatch.utils_logging import _run_logged_image
    first, second = (str(tmp_path / name) for name in ('first.json', 'second.json'))
    release, started = Event(), Event()
    def write(path):
        Path(path).write_text('{}')
        return path
    def block_then_fail(path):
        Path(path).write_text('{')
        started.set()
        assert release.wait(10)
        raise RuntimeError('second image failed')
    with ThreadPoolExecutor(2) as pool:
        executor = OutputExecutor(pool, SpectralMatchPlugin(),
            OutputPolicy(check_validity='around', validation_paths=(first, second)), (first, second))
        a = executor.submit(_run_logged_image, write, 'input_a', first, (first,))
        b = executor.submit(_run_logged_image, block_then_fail, 'input_b', second, (second,))
        try:
            assert started.wait(10)
            assert a.result(timeout=10)[0] == first
            assert Path(first).exists()
            assert not Path(second).exists()
            assert Path(writing_path(second, '_writing')).exists()
        finally:
            release.set()
        with pytest.raises(RuntimeError, match='second image failed'):
            b.result(timeout=10)


def test_overviews_are_created_before_publication(tmp_path, make_test_raster):
    source = str(make_test_raster(tmp_path / 'source.tif', width=64, height=64))
    output = str(tmp_path / 'out.tif')
    policy = OutputPolicy(check_validity='around', validation_paths=(output,),
                          overview_paths=(output,), overview_scales=[2, 4])
    policy.run(shutil.copyfile, (source, output), outputs=(output,))
    with gdal.Open(output) as ds:
        assert ds.GetRasterBand(1).GetOverviewCount() == 2
    assert not Path(writing_path(output, '_writing')).exists()


def test_gpkg_missing_table_is_not_reusable(tmp_path):
    import sqlite3
    p = tmp_path / 'bad.gpkg'
    with sqlite3.connect(p) as db:
        db.execute('CREATE TABLE gpkg_contents (table_name TEXT)')
        db.execute("INSERT INTO gpkg_contents VALUES ('footprints')")
    valid, reason = validation._gpkg_is_valid(str(p))
    assert not valid and 'missing table: footprints' in reason


@pytest.mark.parametrize('mode', ['around', 'before_reuse'])
def test_workflow_validation_timing_and_reuse(tmp_path, monkeypatch, mode):
    source, output = tmp_path / 'input.json', tmp_path / 'output.json'
    source.write_text('{bad')
    calls = []
    def write(input_path, output_path):
        calls.append(output_path)
        Path(output_path).write_text('{}')
        return output_path
    install_function(monkeypatch, 'validated_writer', write, input_paths=('input_path',), output_paths=('output_path',))
    config = {
        'settings': {'plugin': 'shared', 'core:run': True, 'core:check_validity': mode,
                     'core:show_progress': False, 'core:save_statistics_path': None, 'core:load_statistics_path': None},
        'write': {'plugin': 'validated_writer', 'core:run': True, 'core:require_outputs': True,
                  'param:input_path': str(source), 'param:output_path': str(output)},
    }
    if mode == 'around':
        with pytest.raises(RuntimeError, match='invalid inputs'):
            Workflow(config).run()
        assert not calls
        source.write_text('{}')
    Workflow(config).run()
    assert calls == [str(tmp_path / 'output_writing.json')]
    Workflow(config).run()
    assert len(calls) == 1
    output.write_text('{partial')
    Workflow(config).run()
    assert len(calls) == 2
    assert output.read_text() == '{}'


@pytest.mark.parametrize("workers", [1, 2])
def test_native_aggregate_workers_publish_images_with_overviews(tmp_path, make_test_raster, workers):
    sources = [str(make_test_raster(tmp_path / f'in{i}.tif', width=32, height=32)) for i in range(2)]
    outputs = [str(tmp_path / f'out{i}.tif') for i in range(2)]
    config = {
        'shared': {'plugin': 'shared', 'core:run': True, 'core:check_validity': 'around',
                   'core:concurrent_processing': workers, 'core:show_progress': True,
                   'core:save_statistics_path': None, 'core:load_statistics_path': None,
                   'param:window_scales': [2, 4]},
        'math': {'plugin': 'band_math', 'core:run': True, 'core:require_outputs': True,
                 'core:calculate_overviews': True, 'param:input_images': sources,
                 'param:output_images': outputs, 'param:threshold_math': 'B1 > 0'},
    }
    workflow = Workflow(config)
    workflow.run()
    for output in outputs:
        assert Path(output).exists() and not Path(writing_path(output, '_writing')).exists()
        with gdal.Open(output) as ds:
            assert ds.GetRasterBand(1).GetOverviewCount() == 2
    snapshot = workflow.get_progress()
    overview = snapshot['rows'][0]['children'][0]
    assert overview['done'] == 2
