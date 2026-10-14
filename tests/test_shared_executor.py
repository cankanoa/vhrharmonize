"""Workflow-owned aggregate pools and independent per-image threading."""
import os
from pathlib import Path

import pytest

from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.config import validate_config
from workflow_helpers import install_function, import_settings


@pytest.mark.parametrize('workers', [1, 2])
def test_aggregate_uses_parent_coordinator_and_core_pool(tmp_path, monkeypatch, workers):
    seen = {}
    def aggregate(values, *, worker_pool=None, window_threads=None, image_threads=None):
        seen.update(pid=os.getpid(), executor=worker_pool, tiles=window_threads, images=image_threads)
        return [worker_pool.submit(abs, value).result(timeout=20) for value in values]
    install_function(monkeypatch, 'aggregate_pool', aggregate, scope='aggregate',
                     executor_parameter='worker_pool', tile_threads_parameter='window_threads')
    source = tmp_path / 'source.txt'
    source.write_text('input')
    config = {
        'shared': {'plugin': 'shared', 'core:run': True, 'core:concurrent_processing': workers,
                   'core:tile_threads': 3, 'core:show_progress': False, 'core:save_statistics_path': None},
        'files': import_settings(source, tmp_path),
        'aggregate': {'plugin': 'aggregate_pool', 'core:run': True, 'core:require_outputs': True,
                      'param:values': [-2, -3], 'param:image_threads': 8, 'param:window_threads': 7,
                      'const:result': 'returned:$'},
    }
    workflow = Workflow(config, config_dir=tmp_path)
    workflow.run()
    assert seen['pid'] == os.getpid()
    assert seen['executor']._max_workers == workers
    assert seen['tiles'] == 3 and seen['images'] == 8
    assert all('worker_pool' not in node.params for node in workflow.nodes)
    with pytest.raises(RuntimeError):
        seen['executor'].submit(abs, -1)  # The workflow, not the function, owns shutdown.


@pytest.mark.parametrize('value', [0, -1, True, 'bad', [], {}])
def test_tile_threads_rejects_invalid_values(value):
    with pytest.raises(ValueError, match='tile_threads'):
        validate_config({'shared': {'plugin': 'shared', 'core:run': True, 'core:tile_threads': value}})


def test_real_spectralmatch_function_uses_borrowed_workflow_pool(tmp_path, monkeypatch, make_test_raster):
    from vhrharmonize.plugins.compute_overviews import compute_overviews
    from spectralmatch import utils_multiprocessing as parallel
    from vhrharmonize.workflow import engine
    from contextlib import contextmanager
    submitted = []
    original = engine.image_executor
    @contextmanager
    def observe(*args):
        with original(*args) as pool:
            submit = pool.submit
            def record(function, *args, **kwargs):
                submitted.append(function.__name__)
                return submit(function, *args, **kwargs)
            pool.submit = record
            yield pool
    monkeypatch.setattr(engine, 'image_executor', observe)
    monkeypatch.setattr(parallel, '_get_executor', lambda *a, **kw: pytest.fail('nested pool'))
    source = make_test_raster(tmp_path / 'source.tif', width=16, height=16)
    config = {
        'shared': {'plugin': 'shared', 'core:run': True, 'core:concurrent_processing': 2,
                   'core:tile_threads': 1, 'core:show_progress': False, 'core:save_statistics_path': None},
        'files': import_settings(source, tmp_path),
        'overview': {'plugin': 'compute_overviews', 'core:run': True, 'core:require_outputs': True,
                     'param:input_images_paths': [str(source)], 'param:window_scales': [2, 4]},
    }
    Workflow(config, config_dir=tmp_path).run()
    assert '_run_logged_image' in submitted
    from osgeo import gdal
    dataset = gdal.Open(str(source))
    assert dataset.GetRasterBand(1).GetOverviewCount() == 2


def test_dask_executor_runs_spectral_jobs_and_remains_caller_owned():
    from distributed import Client, LocalCluster
    from spectralmatch.utils_multiprocessing import _run_image_tasks
    from vhrharmonize.workflow.concurrency import DaskExecutor
    with LocalCluster(n_workers=2, threads_per_worker=1, processes=False, dashboard_address=None) as cluster:
        with Client(cluster) as client:
            executor = DaskExecutor(client)
            result = _run_image_tasks(abs, [(-2,), (-3,)], input_paths=['a', 'b'],
                                      output_paths=[None, None], executor=executor)
            assert result == [2, 3]
            assert executor.submit(abs, -4).result(timeout=10) == 4
            with pytest.raises(TypeError):
                _run_image_tasks(abs, [('bad',)], input_paths=['c'], output_paths=[None], executor=executor)
            assert client.submit(abs, -5).result(timeout=10) == 5


def test_worldview_uses_only_core_scheduling_controls():
    from vhrharmonize.workflow.config import load_config
    config = load_config(Path(__file__).parents[1] / 'configs/example.worldview.yml')
    assert config['shared']['core:tile_threads'] == 1
    assert config['shared']['core:concurrent_processing'] == 1
    old = {'param:' + key for key in ('image_threads', 'tile_threads', 'concurrent_processing',
                                      'concurrent_processing_backend', 'dask_scheduler')}
    assert all(not old.intersection(settings) for settings in config.values())


def test_concurrency_injection_requires_explicit_declarations():
    from vhrharmonize.plugins.base import FunctionPlugin
    plugin = FunctionPlugin()
    plugin.function = lambda: (lambda executor=None, tile_threads=9: None)
    assert plugin.runtime_parameters(tile_threads=3) == {}
    plugin.tile_threads_parameter = "tile_threads"
    assert plugin.runtime_parameters(tile_threads=3) == {"tile_threads": 3}
    for value in ("cpu", "num_cpu"):
        assert plugin.runtime_parameters(tile_threads=value)["tile_threads"] == (os.cpu_count() or 1)


@pytest.mark.parametrize("declaration", ["missing", "bad-name", 4])
def test_invalid_concurrency_declarations_fail_clearly(declaration):
    from vhrharmonize.plugins.base import FunctionPlugin
    plugin = FunctionPlugin()
    plugin.function = lambda: (lambda executor=None: None)
    plugin.executor_parameter = declaration
    with pytest.raises(ValueError, match="executor_parameter"):
        plugin.runtime_parameters(tile_threads=1)


def test_spectralmatch_and_seamline_adapters_disable_internal_pools():
    from vhrharmonize.workflow.registry import load_plugin
    from vhrharmonize.parameters import function_parameters
    names = "joint_coregistration global_regression local_block_adjustment align_rasters create_footprints postprocess_footprints markov_triangles mask_rasters merge_rasters band_math create_cloud_mask_with_omnicloudmask process_raster_values_to_vector_polygons compute_overviews seamline_metadata".split()
    for name in names:
        plugin = load_plugin(name)
        values = plugin.runtime_parameters(tile_threads=3)
        assert values[plugin.executor_parameter] is None
        assert values['concurrent_processing_backend'] == 'process_pool'
        assert values['dask_scheduler'] is None
        assert set(values) <= set(function_parameters(plugin.function()))
        if name == 'seamline_metadata':
            assert values['concurrent_processing'] == 1
        else:
            assert values['image_threads'] is None


def test_undeclared_function_keeps_its_own_concurrency_arguments(tmp_path, monkeypatch):
    from vhrharmonize.workflow import engine
    seen = {}
    def aggregate(*, executor='native', tile_threads=9):
        seen.update(executor=executor, tiles=tile_threads)
    install_function(monkeypatch, 'unmanaged_pool', aggregate, scope='aggregate')
    monkeypatch.setattr(engine, 'image_executor', lambda *args: pytest.fail('undeclared pool'))
    source = tmp_path / 'source.txt'
    source.write_text('input')
    config = {
        'shared': {'plugin': 'shared', 'core:run': True, 'core:concurrent_processing': 1,
                   'core:tile_threads': 3, 'core:show_progress': False, 'core:save_statistics_path': None},
        'files': import_settings(source, tmp_path),
        'aggregate': {'plugin': 'unmanaged_pool', 'core:run': True, 'core:require_outputs': True,
                      'param:tile_threads': 7},
    }
    Workflow(config, config_dir=tmp_path).run()
    assert seen == {'executor': 'native', 'tiles': 7}
