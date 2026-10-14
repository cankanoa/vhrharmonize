"""Nested counters, worker isolation, transport and terminal hierarchy."""
from io import StringIO
from time import monotonic

import pytest
from prompt_toolkit.formatted_text import to_plain_text

from vhrharmonize.workflow.progress import ProgressState, StepProgress, TaskProgress, WorkflowProgress
from vhrharmonize.workflow.progress_terminal import TerminalProgressDisplay
from vhrharmonize.progress import validate_progress_snapshot


def collector():
    progress = WorkflowProgress.__new__(WorkflowProgress)
    progress.state = ProgressState(4)
    progress.state.rows['joint'] = StepProgress('joint', run=3, all=3)
    progress.state.tasks['task'] = TaskProgress('joint', 'all scenes', weight=3, started=monotonic())
    progress.run_id = 'test'
    progress.snapshot_status = 'running'
    return progress


def emit(progress, identity, name, *, parent=None, kind='jobs', event='start', total=5, done=0, **extra):
    stats = dict(n=done, total=total, prefix=name, unit='pairs', elapsed=1, rate=None,
                 substep=dict(id=identity, name=name, parent=parent, kind=kind, event=event,
                              total=total, done=done, unit='pairs', **extra))
    progress.state.handle(dict(kind='progress', task='task', stats=stats))


def test_pair_counts_active_workers_and_recursive_rows_do_not_change_scene_totals():
    p = collector()
    emit(p, 'root', 'joint_coregistration', kind='group', total=None)
    emit(p, 'group', 'Matching', parent='root', kind='group', total=None)
    emit(p, 'pairs', 'Tie points', parent='group', reused=2)
    for job in ('0', '1'):
        emit(p, 'pairs', 'Tie points', parent='group', event='job_start', job=job, scene='pair ' + job)
        p.state.handle(dict(kind='progress', task='task', stats=dict(n=1,total=10,prefix='Tiles',unit='tiles',elapsed=1,rate=1,
                            job_id='pairs:' + job, parent_substep='pairs', scene='pair ' + job)))
    data = validate_progress_snapshot(p.snapshot())
    row = data['rows'][0]
    leaf = row['children'][0]['children'][0]
    assert (row['all'], row['active'], data['total']['all'], data['total']['active']) == (3, 2, 3, 2)
    assert (leaf['all'], leaf['run'], leaf['reused'], leaf['active']) == (7, 5, 2, 2)
    assert len(data['active']) == 2
    assert {r['scene'] for r in data['active']} == {'pair 0', 'pair 1'}
    emit(p, 'pairs', 'Tie points', parent='group', event='advance', done=1, job='0')
    # Late starts/updates cannot resurrect an already committed job.
    emit(p, 'pairs', 'Tie points', parent='group', event='job_start', job='0')
    data = p.snapshot()
    leaf = data['rows'][0]['children'][0]['children'][0]
    assert (leaf['done'], leaf['active'], data['rows'][0]['active']) == (1, 1, 1)
    assert len(data['active']) == 1
    display = TerminalProgressDisplay(stream=StringIO())
    lines = display.table(data, width=180)
    text = [to_plain_text(line) for line in lines]
    assert text[2].startswith('joint')
    assert text[3].startswith('  Matching')
    assert text[4].startswith('    Tie points')
    assert '4(57%)' in text[4]  # Five overlap jobs minus one completed; seven including loaded.
    assert all('italic' in style for style, _ in lines[3] + lines[4])
    assert all('italic' not in style for style, _ in lines[2])
    p.state.complete('task')
    final = p.snapshot()
    assert final['rows'][0]['active'] == 0 and not final['active']
    assert final['rows'][0]['children'][0]['children'][0]['active'] == 0


def test_failed_phase_preserves_unfinished_work_and_child_validation():
    p = collector()
    emit(p, 'phase', 'Outputs', total=4)
    emit(p, 'phase', 'Outputs', total=4, event='advance', done=1, job='0')
    emit(p, 'phase', 'Outputs', total=4, event='failed', done=1)
    p.state.handle(dict(kind='failed', task='task'))
    data = p.snapshot()
    row = data['rows'][0]['children'][0]
    assert (row['status'], row['done'], row['run'], row['active']) == ('failed', 1, 4, 0)
    data['rows'][0]['children'][0]['all'] = -1
    with pytest.raises(ValueError):
        validate_progress_snapshot(data)


def test_compact_view_keeps_ancestors_and_does_not_merge_repeated_names():
    p = collector()
    emit(p, 'group', 'Repeated', kind='group')
    emit(p, 'nested', 'Repeated', parent='group', kind='group')
    emit(p, 'phase', 'Repeated', parent='nested')
    emit(p, 'phase', 'Repeated', parent='nested', event='job_start', job='0')
    data = p.snapshot()
    display = TerminalProgressDisplay(stream=StringIO())
    for count in (2, 3, 4):
        lines = [to_plain_text(line) for line in display.table(data, width=90, max_rows=count)]
        if any(line.startswith('    ') and 'Repeated' in line for line in lines):
            assert any(line.startswith('  Repeated') for line in lines)
            assert any(line.startswith('joint') for line in lines)


def _overlap_work(ready, value):
    from pathlib import Path
    from time import sleep, monotonic
    from spectralmatch.utils_progress import current_callback
    deadline = monotonic() + 20
    while not Path(ready).exists():
        if monotonic() >= deadline:
            raise RuntimeError('Both worker-start events were not delivered')
        sleep(.01)
    current_callback()(n=1, total=2, prefix='Matching tiles', unit='tiles', elapsed=1, rate=1)
    return value * 2


@pytest.mark.parametrize('backend', ['process_pool', 'dask'])
def test_live_nested_pair_progress_survives_worker_transport(tmp_path, backend):
    from contextlib import nullcontext
    from pathlib import Path
    from spectralmatch.utils_progress import operation
    from spectralmatch.utils_multiprocessing import _run_image_tasks
    from vhrharmonize.workflow.concurrency import image_executor
    from distributed import LocalCluster
    p = collector()
    ready = tmp_path / 'both_started'
    observed = []
    def receive(**stats):
        p.state.handle(dict(kind='progress', task='task', stats=stats))
        snapshot = p.snapshot()
        observed.append(snapshot)
        if snapshot['rows'][0]['active'] == 2:
            ready.touch()
    cluster_context = LocalCluster(n_workers=2, threads_per_worker=1, processes=False, dashboard_address=None) if backend == 'dask' else nullcontext(None)
    with cluster_context as cluster:
        scheduler = ('address', cluster.scheduler_address) if cluster else None
        with image_executor(1 if cluster else 2, backend, scheduler) as pool:
            with operation('joint_coregistration', callback=receive):
                with operation('Matching overlaps', callback=receive):
                    # operation(callback=...) does not replace current_callback;
                    # explicitly bind it for child phases as a public API does.
                    from spectralmatch.utils_progress import progress_context
                    with progress_context(receive):
                        results = _run_image_tasks(_overlap_work, [(str(ready), 1), (str(ready), 2)],
                            input_paths=[('a.tif','b.tif'), ('b.tif','c.tif')], output_paths=[None,None],
                            executor=pool, reused=1)
    assert results == [2,4]
    assert ready.exists()
    assert any(len(s['active']) == 2 for s in observed)
    row = p.snapshot()['rows'][0]['children'][0]['children'][0]
    assert (row['all'], row['reused'], row['done'], row['active'], row['unit']) == (3,1,2,0,'pairs')
    assert any(s['rows'][0]['children'][0]['children'][0]['done'] == 1 for s in observed
               if s['rows'][0].get('children') and s['rows'][0]['children'][0].get('children'))
    validate_progress_snapshot(p.snapshot())


def test_nested_jobs_keep_identity_without_double_counting_outer_workers():
    from vhrharmonize.io.substeps import JobReporter
    p = collector()
    emit(p, 'outer', 'Images', total=2)
    emit(p, 'outer', 'Images', total=2, event='job_start', job='0')
    emit(p, 'inner', 'Windows', parent='outer', total=3)
    emit(p, 'inner', 'Windows', parent='outer', total=3, event='job_start', job='0')
    emit(p, 'inner', 'Windows', parent='outer', total=3, event='job_start', job='1')
    task = p.state.tasks['task']
    assert task.substeps.active_workers() == 1
    events = []
    outer = JobReporter(lambda **stats: events.append(stats), {'id': 'outer'}, '0', 'image')
    inner = JobReporter(outer, {'id': 'inner'}, '1', 'window')
    inner(n=2, total=3, prefix='Matching', unit='tiles', elapsed=1, rate=2)
    assert events[0]['job_id'] == 'inner:1'
    assert events[0]['parent_substep'] == 'inner'
    assert events[0]['scene'] == 'window'
    emit(p, 'outer', 'Images', total=2, event='job_finished', job='0')
    assert task.substeps.active_workers() == 2
