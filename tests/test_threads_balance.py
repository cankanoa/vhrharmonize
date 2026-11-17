"""Idle image slots become tile workers without changing image pool ownership."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier
import time

import pytest

from vhrharmonize.workflow.config import DEFAULT_SHARED, validate_config
from vhrharmonize.workflow.concurrency import balanced_tile_threads, image_executor
from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import import_settings, install_function


@pytest.mark.parametrize('base,remaining,expected', [(2, 8, 2), (2, 20, 2), (2, 2, 5),
    (2, 1, 9), (None, 1, None), ('cpu', 1, 'cpu'), ('num_cpu', 1, 'num_cpu')])
def test_spare_slots_are_shared_without_expanding_all_cpu_requests(base, remaining, expected):
    assert balanced_tile_threads(base, 8, remaining) == expected


@pytest.mark.parametrize('invalid', [None, 1, 'true'])
def test_balance_requires_boolean(invalid):
    with pytest.raises(ValueError, match='threads_balance'):
        validate_config({'shared': {'plugin': 'shared', 'core:run': True, 'core:threads_balance': invalid}})
    assert DEFAULT_SHARED['threads_balance'] is False


@pytest.mark.parametrize('direction', ['horizontal', 'vertical'])
@pytest.mark.parametrize('enabled', [False, True])
def test_scene_jobs_receive_spare_capacity(tmp_path, monkeypatch, direction, enabled):
    from vhrharmonize.workflow import engine
    monkeypatch.setattr(engine, 'ProcessPoolExecutor', lambda max_workers, **kw: ThreadPoolExecutor(max_workers))
    for name in ('a', 'b'):
        (tmp_path / f'{name}.txt').write_text(name)
    seen, barrier = [], Barrier(2)
    def process(input_path, output_path, tile_threads=1):
        seen.append(tile_threads)
        barrier.wait(timeout=10)
        Path(output_path).write_text(str(tile_threads))
    install_function(monkeypatch, 'balanced_scene', process, input_paths=('input_path',),
                     output_paths=('output_path',), tile_threads_parameter='tile_threads')
    config = {
        'shared': {'plugin': 'shared', 'core:run': True, 'core:threads_balance': enabled,
                   'core:concurrent_processing': 4, 'core:tile_threads': 2,
                   'core:processing_direction': direction, 'core:show_progress': False,
                   'core:save_statistics_path': None, 'core:load_statistics_path': None},
        'files': import_settings(tmp_path/'*.txt', tmp_path),
        'process': {'plugin': 'balanced_scene', 'core:run': True, 'core:require_outputs': True,
                    'param:input_path': 'var:mul',
                    'param:output_path': "expr:const.output_dir & '/' & var.basename & '.txt'"},
    }
    Workflow(config, config_dir=tmp_path).run()
    assert seen == [3, 3] if enabled else seen == [2, 2]


def _remaining_worker(index, started_path, tile_threads=2):
    from spectralmatch.utils_thread_balance import current_tile_threads
    before = current_tile_threads(tile_threads)
    if index == 0:
        Path(started_path).write_text('ready')
        deadline = time.monotonic() + 20
        while time.monotonic() < deadline:
            after = current_tile_threads(tile_threads)
            if after > before:
                from spectralmatch.utils_multiprocessing import _get_executor
                with _get_executor('thread', tile_threads) as pool:
                    assert pool._max_workers == after
                    assert pool.submit(abs, -1).result(timeout=10) == 1
                return before, after
            time.sleep(.02)
        raise AssertionError('Remaining worker did not receive released capacity')
    deadline = time.monotonic() + 20
    while not Path(started_path).exists():
        assert time.monotonic() < deadline
        time.sleep(.02)
    return before, before


def test_spectralmatch_running_worker_refreshes_capacity_between_phases(tmp_path):
    from spectralmatch.utils_multiprocessing import _run_image_tasks
    from spectralmatch.utils_thread_balance import current_tile_threads
    with image_executor(4, 'process_pool', threads_balance=True, tile_threads=2) as pool:
        result = _run_image_tasks(_remaining_worker, [(0, str(tmp_path/'started'), 2), (1, str(tmp_path/'started'), 2)],
                                  input_paths=['a', 'b'], output_paths=[None, None], executor=pool)
        assert result[0] == (3, 5)
        assert pool.submit(abs, -1).result(timeout=10) == 1
    assert current_tile_threads(2) == 2


def test_single_mosaic_gets_all_unused_slots(tmp_path, monkeypatch, make_test_raster):
    from spectralmatch import utils
    seen = []
    original = utils._resolve_parallel_config
    def record(value, *args):
        seen.append(value)
        return original(value, *args)
    monkeypatch.setattr(utils, '_resolve_parallel_config', record)
    source = str(make_test_raster(tmp_path/'input.tif', width=32, height=32))
    output = str(tmp_path/'mosaic.tif')
    config = {
        'shared': {'plugin': 'shared', 'core:run': True, 'core:threads_balance': True,
                   'core:concurrent_processing': 4, 'core:tile_threads': 2,
                   'core:show_progress': False, 'core:save_statistics_path': None,
                   'core:load_statistics_path': None},
        'merge': {'plugin': 'merge_rasters', 'core:run': True, 'core:require_outputs': True,
                  'param:input_images': [source], 'param:output_image_path': output},
    }
    Workflow(config, config_dir=tmp_path).run()
    assert 5 in seen and Path(output).is_file()
