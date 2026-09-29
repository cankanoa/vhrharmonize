"""Progress callbacks and dashboard accounting across execution backends."""

from concurrent.futures import Future, ThreadPoolExecutor
from io import StringIO
from pathlib import Path
from threading import Event
from types import ModuleType
import inspect
import json
import os
import sys
from time import monotonic, sleep
import yaml

import pytest
from rich.console import Console

from vhrharmonize.io.progress import current_callback, progress, reports_progress
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.progress import ProgressState, StepProgress, TaskProgress, WorkflowProgress
from vhrharmonize.workflow.progress_rich import RichProgressDisplay, _bar
from vhrharmonize.progress import read_progress_snapshot, render_progress, validate_progress_snapshot
from workflow_helpers import copy_step, import_settings, install_function


@pytest.fixture
def dashboards(monkeypatch):
    result, displays = [], []
    initialize = WorkflowProgress.__init__
    initialize_display = RichProgressDisplay.__init__

    def capture_display(self, **kwargs):
        initialize_display(self, console=Console(file=StringIO(), width=120))
        displays.append(self)

    def capture(self, workflow, **kwargs):
        initialize(self, workflow, **kwargs)
        self.test_display = displays[-1] if displays else None
        result.append(self)

    monkeypatch.setattr(RichProgressDisplay, "__init__", capture_display)
    monkeypatch.setattr(WorkflowProgress, "__init__", capture)
    return result


def recipe(tmp_path, *, direction="vertical", workers=1):
    for name in ("a", "b", "c"):
        source = tmp_path / "source" / f"{name}.txt"
        source.parent.mkdir(exist_ok=True)
        source.write_text(name)
    (tmp_path / "temp").mkdir()
    (tmp_path / "output").mkdir()
    (tmp_path / "temp/a_first.txt").write_text("a")
    (tmp_path / "output/c_second.txt").write_text("c")
    return {
        "shared": {"plugin": "shared", "core:run": True, "core:show_progress": True,
                   "core:log_to_console": True, "core:processing_direction": direction,
                   "core:concurrent_processing": workers},
        "files": import_settings(tmp_path / "source/*.txt", tmp_path),
        "prepare": copy_step("first", "mul", suffix="_first"),
        "deliver": copy_step("second", "first", 'expr:const.output_dir & "/" & var.basename & "_second.txt"'),
    }


@pytest.mark.parametrize("direction", ["horizontal", "vertical"])
@pytest.mark.parametrize("workers", [1, 2])
def test_counts_reuse_and_totals_with_real_workers(tmp_path, dashboards, direction, workers):
    config = recipe(tmp_path, direction=direction, workers=workers)
    config["shared"].pop("core:show_progress")  # Exercise the enabled default.
    Workflow(config).run()
    dashboard = dashboards[-1]
    state = dashboard.state
    assert list(state.rows) == ["prepare", "deliver"]
    assert (state.rows["prepare"].done, state.rows["prepare"].run, state.rows["prepare"].all) == (1, 1, 3)
    assert (state.rows["deliver"].done, state.rows["deliver"].run, state.rows["deliver"].all) == (2, 2, 3)
    assert state.rows["prepare"].reused == state.rows["deliver"].reused == 1
    assert (state.total().done, state.total().run, state.total().all, state.total().reused) == (3, 3, 6, 2)
    assert state.active() == 0
    assert all(t.duration is not None and t.done for t in state.tasks.values())
    output = dashboard.test_display.console.file.getvalue()
    assert "file_source" not in output and "core:" not in output
    assert "Unused" in output and "Done" in output and "Run" in output and "All" in output
    assert "3(50%)" in output and "6(100%)" in output
    assert state.rows["prepare"].unused == 1
    assert "gray: reused" not in output
    assert "\x1b[" not in output  # Redirected output gets one plain final summary.
    assert not dashboard.reader.is_alive()
    assert [(tmp_path / "output" / f"{name}_second.txt").read_text() for name in ("a", "b", "c")] == ["a", "b", "c"]
    Workflow(config).run()
    assert dashboards[-1].state.total().run == 0
    assert not dashboards[-1].state.tasks


def test_disable_does_not_start_dashboard(tmp_path, dashboards):
    config = recipe(tmp_path)
    config["shared"].update({"core:show_progress": False, "core:log_to_console": False})
    Workflow(config).run()
    assert dashboards == []


def test_forced_reprocessing_is_not_counted_as_reused(tmp_path, dashboards):
    config = recipe(tmp_path)
    config["shared"]["core:run_from_existing"] = False
    Workflow(config).run()
    state = dashboards[-1].state
    assert state.total().reused == 0
    assert state.total().done == state.total().run == state.total().all == 6


def test_callback_uses_tqdm_fields_and_restores_context_on_error():
    snapshots = []

    @reports_progress
    def calculate(fail=False):
        for i in progress(range(3), desc="matching tiles", unit="tiles", mininterval=0):
            if fail and i == 1:
                raise RuntimeError("tile failed")
        return 42

    assert "progress_callback" in inspect.signature(calculate).parameters
    assert calculate(progress_callback=lambda **stats: snapshots.append(stats)) == 42
    tiles = [s for s in snapshots if s["prefix"] == "matching tiles"]
    assert [s["n"] for s in tiles][0] == 0
    assert tiles[-1]["n"] == tiles[-1]["total"] == 3
    assert {"n", "total", "prefix", "unit", "elapsed", "rate"} <= tiles[-1].keys()
    assert len({s["operation"] for s in tiles}) == 1
    snapshots.clear()
    with pytest.raises(RuntimeError, match="tile failed"):
        calculate(True, progress_callback=lambda **stats: snapshots.append(stats))
    assert snapshots[-1]["status"] == "failed"
    assert current_callback() is None


def test_worker_callbacks_and_messages_share_parent_display(tmp_path, monkeypatch, dashboards):
    from vhrharmonize.workflow import engine

    monkeypatch.setattr(engine, "ProcessPoolExecutor", ThreadPoolExecutor)
    config = recipe(tmp_path, workers=2)

    def call(value, progress_callback=None):
        assert callable(progress_callback)
        print(f"processing {value}")
        for n in range(3):
            progress_callback(n=n, total=2, prefix="tiles", unit="tiles", elapsed=n, rate=1)
        return value

    install_function(monkeypatch, "callback_plugin", call)
    config["measure"] = {"plugin": "callback_plugin", "core:run": True, "param:value": "var:basename"}
    Workflow(config).run()
    state = dashboards[-1].state
    assert state.rows["measure"].done == 3
    messages = [str(m) for m in state.messages]
    assert all(f"processing {value}" in messages for value in ("a", "b", "c"))
    assert state.active() == 0


@pytest.mark.parametrize("direction", ["horizontal", "vertical"])
def test_dask_reports_events_and_removes_subscription(tmp_path, monkeypatch, dashboards, direction):
    handlers, subscriptions = {}, []
    clients = []

    class Client:
        def __init__(self, address):
            clients.append(self)

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def subscribe_topic(self, topic, handler):
            subscriptions.append(topic)
            handlers[topic] = handler

        def unsubscribe_topic(self, topic):
            handlers.pop(topic)

        def nthreads(self):
            return {"worker": 2}

        def submit(self, function, payload, **kwargs):
            future = Future()
            future.set_result(function(payload))
            return future

        def cancel(self, futures):
            pass

        def log_event(self, topic, event):
            handlers[topic]((0, event))

    distributed = ModuleType("dask.distributed")
    distributed.Client = Client
    distributed.as_completed = lambda futures: iter(reversed(list(futures)))
    distributed.get_worker = lambda: clients[-1]
    monkeypatch.setitem(sys.modules, "dask.distributed", distributed)
    monkeypatch.setitem(sys.modules, "distributed", distributed)
    config = recipe(tmp_path, direction=direction)
    config["shared"].update({"core:concurrent_processing_backend": "dask",
                              "core:dask_scheduler_address": "tcp://scheduler:8786"})
    Workflow(config).run()
    state = dashboards[-1].state
    assert state.total().done == state.total().run == 3
    assert state.active() == 0 and state.workers == 2
    assert subscriptions and not handlers


def test_failed_call_does_not_advance_and_restores_streams(tmp_path, monkeypatch, dashboards):
    def fail(value):
        print("before failure")
        raise RuntimeError("backend failed")

    install_function(monkeypatch, "fail", fail)
    config = recipe(tmp_path)
    config["failure"] = {"plugin": "fail", "core:run": True, "param:value": "var:basename"}
    stdout, stderr = sys.stdout, sys.stderr
    with pytest.raises(RuntimeError, match="backend failed"):
        Workflow(config).run()
    state = dashboards[-1].state
    assert state.rows["failure"].done == 0
    assert state.failed and state.active() == 0
    assert (sys.stdout, sys.stderr) == (stdout, stderr)
    assert "backend failed" in dashboards[-1].test_display.console.file.getvalue()


def test_overviews_have_separate_counts_after_named_step(tmp_path, make_test_raster, dashboards):
    source = make_test_raster(tmp_path / "source.tif")
    config = {
        "shared": {"plugin": "shared", "core:run": True, "core:show_progress": True, "param:window_scales": [2]},
        "files": import_settings(source, tmp_path),
        "deliver": {**copy_step("result", "mul", 'expr:const.output_dir & "/result.tif"'), "core:calculate_overviews": True},
    }
    Workflow(config).run()
    state = dashboards[-1].state
    assert list(state.rows) == ["deliver", "deliver / overviews"]
    assert state.rows["deliver / overviews"].done == 1
    assert state.total().done == state.total().run == 2


def test_reused_segment_and_eta_exclude_cached_work():
    state = ProgressState(workers=2)
    row = StepProgress("alignment", run=20, all=264, reused=244, done=3)
    state.rows[row.name] = row
    state.tasks["completed"] = TaskProgress(row.name, "scene", weight=3, done=True, duration=30)
    assert state.remaining_work(row) == 170
    bar = _bar({"all": row.all, "reused": row.reused, "done": row.done}, width=264)
    assert len(bar.plain.split(" ")[0]) == 264
    assert "%" not in bar.plain
    assert any(s.style == "grey50" and s.end - s.start == 244 for s in bar.spans)
    assert any(s.style == "green" and s.end - s.start == 3 for s in bar.spans)


def test_opaque_operation_elapsed_time_keeps_updating(tmp_path, monkeypatch):
    dashboard = WorkflowProgress(Workflow(recipe(tmp_path)))
    console = Console(file=StringIO(), width=120)
    monkeypatch.setattr("vhrharmonize.workflow.progress.monotonic", lambda: 100)
    dashboard.state.tasks["scene"] = TaskProgress("prepare", "scene", started=100)
    dashboard.state.handle({"task": "scene", "kind": "progress", "stats": {
        "n": 0, "total": None, "prefix": "opaque backend", "elapsed": 0, "rate": None,
    }})
    monkeypatch.setattr("vhrharmonize.workflow.progress.monotonic", lambda: 110)
    console.print(render_progress(dashboard.snapshot(), console=console))
    output = console.file.getvalue()
    assert "opaque backend" in output and "working · 10s · ETA estimating" in output


def test_runtime_scene_discovery_replans_dashboard(tmp_path, monkeypatch, dashboards):
    install_function(monkeypatch, "seed", lambda: [{"n": 1}, {"n": 2}], scene_records_return="$")
    install_function(monkeypatch, "double", lambda n: n * 2)
    install_function(monkeypatch, "reset", lambda values: [{"total": sum(values)}], scene_records_return="$")
    install_function(monkeypatch, "consume", lambda value: value)
    config = {
        "shared": {"plugin": "shared", "core:run": True, "core:show_progress": True},
        "source": {"plugin": "seed", "core:run": True},
        "first": {"plugin": "double", "core:run": True, "param:n": "var:n", "var:doubled": "returned:$"},
        "reset": {"plugin": "reset", "core:run": True, "param:values": "collect:doubled"},
        "last": {"plugin": "consume", "core:run": True, "param:value": "var:total"},
    }
    workflow = Workflow(config)
    assert workflow.counts()["last"]["pending"]
    workflow.run()
    state = dashboards[-1].state
    assert not state.rows["last"].pending
    assert (state.rows["last"].done, state.rows["last"].run, state.rows["last"].all) == (1, 1, 1)
    assert state.rows["first"].done == state.rows["reset"].done == 2


def test_real_dask_cluster(tmp_path, monkeypatch, dashboards):
    distributed = pytest.importorskip("distributed")
    delivered = Event()
    handle = ProgressState.handle

    def receive(self, event):
        handle(self, event)
        if event.get("stats", {}).get("prefix") == "remote tiles":
            delivered.set()

    monkeypatch.setattr(ProgressState, "handle", receive)

    def measure(value, progress_callback=None):
        progress_callback(n=1, total=2, prefix="remote tiles", unit="tiles", elapsed=1, rate=1)
        assert delivered.wait(10), "Dask did not deliver progress until after the call finished"
        return value

    install_function(monkeypatch, "measure_remote", measure)
    config = recipe(tmp_path)
    config["measure"] = {"plugin": "measure_remote", "core:run": True, "param:value": "var:basename"}
    with distributed.LocalCluster(n_workers=2, threads_per_worker=1, processes=False,
                                  dashboard_address=None) as cluster:
        config["shared"].update({"core:concurrent_processing_backend": "dask",
                                 "core:dask_scheduler_address": cluster.scheduler_address})
        Workflow(config).run()
    state = dashboards[-1].state
    assert state.total().done == state.total().run == 6
    assert state.active() == 0 and not state.failed


def test_progress_snapshot_is_readable_during_work_and_after_completion(tmp_path, monkeypatch, dashboards):
    from vhrharmonize.workflow.api import run_workflow

    filename = tmp_path / "recipe.yml"
    snapshot = Path(str(filename) + ".progress.json")
    observed = []

    def measure(value, progress_callback=None):
        progress_callback(n=1, total=2, prefix="processing tiles", unit="tiles", elapsed=1, rate=1)
        deadline = monotonic() + 5
        while monotonic() < deadline:
            data = json.loads(snapshot.read_text())
            if any(task["stats"].get("prefix") == "processing tiles" for task in data["active"]):
                observed.append(data)
                break
            sleep(0.05)
        else:
            pytest.fail("No active operation was saved while the function was running")
        return value

    install_function(monkeypatch, "measure_snapshot", measure)
    config = recipe(tmp_path)
    config["measure"] = {"plugin": "measure_snapshot", "core:run": True, "param:value": "var:basename"}
    filename.write_text(yaml.safe_dump(config, sort_keys=False))
    run_workflow(filename)
    final = json.loads(snapshot.read_text())
    assert observed and observed[0]["status"] == "running"
    assert final["status"] == "completed" and not final["active"]
    restored = read_progress_snapshot(snapshot)
    assert restored["total"]["done"] == restored["total"]["run"] == 6
    assert restored["rows"][0]["unused"] == 1
    assert not list(tmp_path.glob("*.progress.json.*.tmp"))


def test_phase_only_plugin_is_labeled_but_explicit_callbacks_are_supported(tmp_path, monkeypatch, dashboards):
    install_function(monkeypatch, "opaque", lambda value: value)
    install_function(monkeypatch, "reporting", lambda value, progress_callback=None: value)
    config = recipe(tmp_path)
    config["settings"] = {"core:run": True, "const:example": 1}
    for name in ("opaque", "reporting"):
        config[name] = {"plugin": name, "core:run": True, "param:value": "var:basename"}
    Workflow(config).run()
    panel = dashboards[-1]
    assert not panel.state.rows["opaque"].worker_progress
    assert panel.state.rows["reporting"].worker_progress
    console = Console(file=StringIO(), width=180)
    console.print(render_progress(panel.snapshot(), console=console))
    assert "opaque (no cb)" in console.file.getvalue()


def test_spectralmatch_aggregate_workers_feed_core_display(tmp_path, make_test_raster, monkeypatch, dashboards):
    import spectralmatch

    if not getattr(spectralmatch.align_rasters, "__worker_progress__", False):
        pytest.skip("Requires SpectralMatch with worker callbacks")
    for name in ("a", "b"):
        make_test_raster(tmp_path / "source" / f"{name}.tif")
    snapshots = []
    handle = ProgressState.handle

    def receive(self, event):
        if event["kind"] == "progress":
            snapshots.append(event["stats"])
        handle(self, event)

    monkeypatch.setattr(ProgressState, "handle", receive)
    config = {
        "shared": {"plugin": "shared", "core:run": True, "core:show_progress": True},
        "files": import_settings(tmp_path / "source/*.tif", tmp_path),
        "resize": {"plugin": "align_rasters", "core:run": True,
                   "param:input_images": "collect:mul",
                   "var:resized": 'expr:const.output_dir & "/" & var.basename & "_resized.tif"',
                   "param:output_images": "collect:resized", "param:image_threads": 2},
    }
    Workflow(config).run()
    row = dashboards[-1].state.rows["resize"]
    assert row.worker_progress and row.done == row.run == row.all == 2
    assert any(s.get("unit") == "images" and s["n"] == s["total"] == 2 for s in snapshots)
    assert any(s.get("scene") in {"a.tif", "b.tif"} and s["prefix"] == "Writing raster" for s in snapshots)


@pytest.mark.parametrize("workers", [1, 2])
def test_public_callback_and_polling_work_without_rich(tmp_path, monkeypatch, workers):
    import builtins

    original_import = builtins.__import__

    def no_rich(name, *args, **kwargs):
        if name == "rich" or name.startswith("rich.") or name.endswith("progress_rich"):
            raise AssertionError("Headless progress must not import the frontend")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_rich)
    config = recipe(tmp_path, workers=workers)
    config["shared"]["core:show_progress"] = False
    workflow = Workflow(config)
    assert workflow.get_progress() is None
    snapshots, parent = [], os.getpid()

    def receive(data):
        assert os.getpid() == parent
        snapshots.append(validate_progress_snapshot(data))
        assert workflow.get_progress()["run_id"] == data["run_id"]
        data["rows"].clear()  # Consumers own their copy, including nested values.

    workflow.run(progress_callback=receive)
    assert snapshots[0]["status"] == "running"
    assert snapshots[0]["total"]["done"] == 0
    assert snapshots[-1]["status"] == "completed"
    final = workflow.get_progress()
    assert len(final["rows"]) == 2
    assert final["total"]["done"] == final["total"]["run"] == 3
    assert final["total"]["all"] == 6 and final["total"]["percentages"]["done"] == 50
    assert final["total"]["fraction_done"] == 1 and final["total"]["eta_seconds"] == 0
    assert not final["active"]
    assert not any("consumer disabled" in message for message in final["messages"])
    final["total"]["done"] = -1
    assert workflow.get_progress()["total"]["done"] == 3


def test_public_callback_receives_live_dask_operations_without_a_display(tmp_path, monkeypatch):
    distributed = pytest.importorskip("distributed")
    delivered, snapshots = Event(), []

    def receive(data):
        snapshots.append(data)
        if any(t["stats"]["prefix"] == "remote tiles" for t in data["active"]):
            delivered.set()

    def measure(value, progress_callback=None):
        progress_callback(n=1, total=2, prefix="remote tiles", unit="tiles", elapsed=1, rate=1)
        assert delivered.wait(10), "Application callback did not receive progress during processing"
        return value

    install_function(monkeypatch, "measure_remote_api", measure)
    config = recipe(tmp_path)
    config["shared"]["core:show_progress"] = False
    config["measure"] = {"plugin": "measure_remote_api", "core:run": True, "param:value": "var:basename"}
    with distributed.LocalCluster(n_workers=2, threads_per_worker=1, processes=False,
                                  dashboard_address=None) as cluster:
        config["shared"].update({"core:concurrent_processing_backend": "dask",
                                 "core:dask_scheduler_address": cluster.scheduler_address})
        Workflow(config).run(progress_callback=receive)
    assert delivered.is_set()
    assert snapshots[-1]["total"]["done"] == snapshots[-1]["total"]["run"] == 6
    active = next(t for s in snapshots for t in s["active"] if t["stats"]["prefix"] == "remote tiles")
    assert active["eta_seconds"] == 1
    assert active["stats"]["n"] == 1 and active["stats"]["total"] == 2


def test_headless_hpc_reporting_and_failed_run_publish_final_snapshot(tmp_path, monkeypatch, capsys):
    from vhrharmonize import run_workflow

    def fail(value):
        raise RuntimeError("processing failed")

    install_function(monkeypatch, "fail_api", fail)
    config = recipe(tmp_path)
    config["shared"].update({"core:show_progress": False, "core:report_progress": True})
    config["failure"] = {"plugin": "fail_api", "core:run": True, "param:value": "var:basename"}
    filename = tmp_path / "headless.yml"
    filename.write_text(yaml.safe_dump(config, sort_keys=False))
    with pytest.raises(RuntimeError, match="processing failed"):
        run_workflow(filename)
    snapshot = read_progress_snapshot(str(filename) + ".progress.json")
    assert snapshot["status"] == snapshot["total"]["status"] == "failed"
    assert snapshot["rows"][-1]["done"] == 0
    assert snapshot["rows"][-1]["status"] == "failed"
    assert not snapshot["active"]
    assert "Workflow progress" not in capsys.readouterr().out
    assert any("processing failed" in m for m in snapshot["messages"])


def test_callback_failure_is_isolated_and_explicit_path_enables_reporting(tmp_path):
    from vhrharmonize import run_workflow

    config = recipe(tmp_path)
    config["shared"]["core:show_progress"] = False
    calls = []

    def broken(snapshot):
        calls.append(snapshot)
        raise RuntimeError("app disconnected")

    filename = tmp_path / "progress.json"
    run_workflow(config, progress_callback=broken, progress_path=filename)
    assert len(calls) == 1
    snapshot = read_progress_snapshot(filename)
    assert snapshot["status"] == "completed" and snapshot["total"]["done"] == 3
    assert any("app disconnected" in m for m in snapshot["messages"])
    with pytest.raises(ValueError, match="Unsupported"):
        validate_progress_snapshot({**snapshot, "version": 999})
    with pytest.raises(ValueError, match="Malformed"):
        validate_progress_snapshot({**snapshot, "total": {}})
    with pytest.raises(ValueError, match="Malformed"):
        validate_progress_snapshot({**snapshot, "total": {**snapshot["total"], "eta_seconds": float("nan")}})


def test_static_renderer_uses_numeric_snapshot_without_recomputing_eta(tmp_path):
    workflow = Workflow(recipe(tmp_path)).plan()
    reporter = WorkflowProgress(workflow)
    reporter.state.rows["prepare"].peak_active = 2
    reporter.state.tasks["sample"] = TaskProgress("prepare", "sample", done=True, weight=1, duration=10)
    data = reporter.snapshot()
    assert data["rows"][0]["eta_seconds"] == 5
    data["rows"][0]["eta_seconds"] = 83
    console = Console(file=StringIO(), width=140)
    console.print(render_progress(data, console=console))
    assert "~1m 23s" in console.file.getvalue()


def test_run_plugin_and_rich_consume_the_same_public_data(tmp_path, dashboards):
    from vhrharmonize import run_plugin

    snapshots = []
    run_plugin("file_source", recipe(tmp_path), progress_callback=snapshots.append)
    assert snapshots[-1]["status"] == "completed"
    assert snapshots[-1] == dashboards[-1].test_display.snapshot
    assert snapshots[-1]["total"]["done"] == 3
