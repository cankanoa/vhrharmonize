"""Core measurements, append-only OpenTelemetry records and pandas reports."""

from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
import json
import os
from pathlib import Path

import pandas as pd
import pytest

from vhrharmonize import StatisticsRecorder, load_statistics, run_workflow, summarize_statistics, validate_statistics_record
from vhrharmonize.workflow.engine import Workflow
from test_workflow_progress import recipe
from workflow_helpers import install_function, stage, transfer


def quiet_recipe(tmp_path, **kwargs):
    config = recipe(tmp_path, **kwargs)
    config["shared"].update({"core:show_progress": False, "core:log_to_console": False,
                              "core:save_statistics_path": "history.jsonl", "core:load_statistics_path": "history.jsonl"})
    return config


def spans(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def event(run="one", kind="task", name="alignment", seconds=2, status="completed", **attributes):
    return {"version": 1, "run_id": run, "run_started_ns": 1_700_000_000_000_000_000,
            "kind": kind, "name": name, "start_time_ns": 1_700_000_001_000_000_000,
            "duration_seconds": seconds, "status": status,
            "attributes": {"vhr.backend": "process_pool", **attributes}}


@pytest.mark.parametrize("direction", ["horizontal", "vertical"])
@pytest.mark.parametrize("workers", [1, 2])
def test_repeated_workflows_append_real_measurements_and_reuse_counts(tmp_path, workers, direction, capsys):
    config = quiet_recipe(tmp_path, workers=workers, direction=direction)
    workflow = Workflow(config, config_dir=tmp_path)
    events, parent = [], os.getpid()

    def consume(data):
        assert os.getpid() == parent
        assert workflow.get_progress() is None or workflow.get_progress()["run_id"] == data["run_id"]
        events.append(deepcopy(data))
        data["attributes"].clear()  # Other consumers receive independent copies.

    workflow.run(event_callback=consume)
    raw = tmp_path / "history.jsonl"
    first_bytes, first = raw.read_bytes(), spans(raw)
    assert capsys.readouterr().out == ""
    assert {s["attributes"]["vhr.run_id"] for s in first} == {workflow.get_progress()["run_id"]}
    assert len(first) == len(events)
    tasks = [s for s in first if s["attributes"]["vhr.kind"] == "task"]
    assert sorted(s["name"] for s in tasks) == ["deliver", "deliver", "files", "prepare"]
    assert all(s["attributes"]["vhr.measured"] and s["attributes"]["vhr.duration_seconds"] > 0 for s in tasks)
    root = next(s for s in first if s["attributes"]["vhr.kind"] == "workflow")
    assert root["status"]["status_code"] == "OK" and root["parent_id"] is None
    assert all(s["parent_id"] == root["context"]["span_id"] for s in first if s is not root)
    assert {s["context"]["trace_id"] for s in first} == {root["context"]["trace_id"]}
    assert {s["name"] for s in first if s["attributes"]["vhr.kind"] == "core"} >= {
        "initialization", "discovery", "build", "planning", "execution"}
    second = Workflow(config, config_dir=tmp_path)
    second.run()
    assert raw.read_bytes().startswith(first_bytes)
    assert len({s["attributes"]["vhr.run_id"] for s in spans(raw)}) == 2
    appended = spans(raw)[len(first):]
    assert [s["name"] for s in appended if s["attributes"]["vhr.kind"] == "task"] == ["files"]
    counts = [s["attributes"] for s in appended if s["attributes"]["vhr.kind"] == "step_summary"]
    assert sum(c["vhr.done"] for c in counts) == sum(c["vhr.run"] for c in counts) == 0
    report = summarize_statistics(raw, tmp_path / "summary.json")
    assert sum(r["samples"] for r in report["data"] if r["kind"] == "task") == 5


@pytest.mark.parametrize("bad_output", [False, True])
def test_failed_calls_and_output_validation_preserve_raw_failure(tmp_path, monkeypatch, bad_output):
    def calculate(input_path, output_path):
        if not bad_output:
            raise RuntimeError("calculation failed")
        return output_path  # A success return without the required output file.

    install_function(monkeypatch, "bad_statistics", calculate, input_paths=("input_path",),
                     output_paths=("output_path",))
    config = quiet_recipe(tmp_path)
    config["deliver"]["plugin"] = "bad_statistics"
    with pytest.raises(Exception):
        Workflow(config, config_dir=tmp_path).run()
    records = spans(tmp_path / "history.jsonl")
    task = next(s for s in records if s["name"] == "deliver" and s["attributes"]["vhr.kind"] == "task")
    assert task["attributes"]["vhr.status"] == "failed"
    assert task["attributes"]["vhr.duration_seconds"] >= 0
    assert records[-1]["attributes"]["vhr.kind"] == "workflow"
    assert records[-1]["status"]["status_code"] == "ERROR"
    assert records[-1]["attributes"]["vhr.status"] == "failed"


def test_event_api_works_independently_and_bad_consumer_does_not_break_recording(tmp_path):
    config = quiet_recipe(tmp_path)
    calls = []

    def broken(data):
        calls.append(data)
        raise RuntimeError("consumer disconnected")

    with pytest.warns(RuntimeWarning, match="consumer disconnected"):
        run_workflow(config, config_dir=str(tmp_path), event_callback=broken)
    assert len(calls) == 1
    assert spans(tmp_path / "history.jsonl")[-1]["attributes"]["vhr.status"] == "completed"
    config["shared"]["core:save_statistics_path"] = None
    original = (tmp_path / "history.jsonl").read_bytes()
    events = []
    run_workflow(config, config_dir=str(tmp_path), event_callback=events.append)
    assert events[-1]["kind"] == "workflow"
    assert (tmp_path / "history.jsonl").read_bytes() == original


def test_dry_run_and_hpc_history_input_never_overwrite_remote_output(tmp_path):
    config = quiet_recipe(tmp_path)
    run_workflow(config, config_dir=str(tmp_path), dry_run=True)
    raw = tmp_path / "history.jsonl"
    assert not raw.exists()
    with StatisticsRecorder(raw) as recorder:
        recorder(event(kind="workflow", name="workflow"))
    before = raw.read_bytes()
    staged, uploads, downloads = stage(config, tmp_path)
    assert raw.read_bytes() == before
    remote = staged["shared"]["core:save_statistics_path"]
    assert remote == downloads[str(raw)]
    assert remote.startswith(str(tmp_path / "remote/output/statistics"))
    loaded = staged["shared"]["core:load_statistics_path"]
    assert loaded == uploads[str(raw)] and loaded != remote
    assert loaded.startswith(str(tmp_path / "remote/reference/statistics"))
    transfer(uploads)
    Workflow(staged, config_dir=tmp_path).run()
    assert spans(Path(remote))[-1]["attributes"]["vhr.status"] == "completed"
    assert raw.read_bytes() == before


def test_statistics_cli_pandas_schema_and_filters(tmp_path, capsys):
    from vhrharmonize.cli.main import main

    raw, output = tmp_path / "history.jsonl", tmp_path / "summary.json"
    with StatisticsRecorder(raw) as recorder:
        for run, duration in (("one", 2), ("two", 4)):
            recorder(event(run=run, seconds=duration, **{"vhr.plugin": "alignment", "vhr.scene_units": 3}))
            recorder(event(run=run, kind="step_summary", seconds=0,
                           **{"vhr.done": 3, "vhr.run": 3, "vhr.all": 10, "vhr.reused": 7}))
            recorder(event(run=run, kind="workflow", name="workflow", seconds=5))
    original = raw.read_bytes()
    assert main(["statistics", "--input-path", str(raw), "--output-path", str(output)]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result == json.loads(output.read_text())
    frame = pd.read_json(output, orient="table")
    row = frame[frame["kind"] == "task"].iloc[0]
    assert (row["samples"], row["runs"], row["unfinished_runs"], row["scene_units"]) == (2, 2, 0, 6)
    assert (row["total_seconds"], row["mean_seconds"], row["min_seconds"], row["max_seconds"]) == (6, 3, 2, 4)
    assert row["p50_seconds"] == 3 and row["p95_seconds"] == pytest.approx(3.9)
    counts = frame[frame["kind"] == "step_summary"].iloc[0]
    assert counts["samples"] == 0 and pd.isna(counts["total_seconds"])
    assert (counts["done"], counts["run"], counts["all"], counts["reused"]) == (6, 6, 20, 14)
    separate = summarize_statistics(raw, output, per_run=True)
    assert {r["run_id"] for r in separate["data"]} == {"one", "two"}
    assert all(r["samples"] == 1 for r in separate["data"] if r["kind"] == "task")
    selected = summarize_statistics(raw, output, run_id="one")
    assert next(r for r in selected["data"] if r["kind"] == "task")["mean_seconds"] == 2
    assert raw.read_bytes() == original
    with pytest.raises(ValueError, match="No statistics"):
        summarize_statistics(raw, output, run_id="missing")


def test_analysis_handles_incomplete_runs_and_duplicate_spans(tmp_path):
    raw = tmp_path / "history.jsonl"
    with StatisticsRecorder(raw) as recorder:
        recorder(event())  # Simulate an application closing without workflow completion.
    original = raw.read_text()
    raw.write_text(original + original)
    data = summarize_statistics(raw, tmp_path / "summary.json")["data"]
    task = next(r for r in data if r["kind"] == "task")
    assert task["samples"] == task["unfinished_runs"] == 1
    root = next(r for r in data if r["kind"] == "workflow")
    assert root["status"] == "incomplete" and root["samples"] == 0 and root["total_seconds"] is None


def test_analysis_freezes_boundary_without_holding_up_appenders(tmp_path):
    from vhrharmonize.statistics import _raw_lines

    raw = tmp_path / "history.jsonl"
    with StatisticsRecorder(raw) as recorder:
        recorder(event(run="one", kind="workflow", name="workflow"))
    reader = _raw_lines(raw)
    assert json.loads(next(reader))["attributes"]["vhr.run_id"] == "one"
    with StatisticsRecorder(raw) as recorder:
        recorder(event(run="two", kind="workflow", name="workflow"))
    assert list(reader) == []
    assert len(spans(raw)) == 2


def test_aggregate_and_overviews_are_measured_by_core(tmp_path, monkeypatch, make_test_raster):
    from workflow_helpers import copy_step, import_settings

    for name in ("a", "b"):
        make_test_raster(tmp_path / "source" / f"{name}.tif")
    seen = []

    def aggregate(values):
        seen.extend(values)
        return values

    install_function(monkeypatch, "statistics_aggregate", aggregate, scope="aggregate")
    config = {
        "shared": {"plugin": "shared", "core:run": True, "core:show_progress": False,
                   "core:log_to_console": False, "core:save_statistics_path": "history.jsonl",
                   "param:window_scales": [2]},
        "files": import_settings(tmp_path / "source/*.tif", tmp_path),
        "copy": {**copy_step("result", "mul", 'expr:const.output_dir & "/" & var.basename & ".tif"', require_outputs=True),
                 "core:calculate_overviews": True},
        "batch": {"plugin": "statistics_aggregate", "core:run": True, "core:require_outputs": True, "param:values": "collect:result"},
    }
    Workflow(config, config_dir=tmp_path).run()
    assert len(seen) == 2
    tasks = [s for s in spans(tmp_path / "history.jsonl") if s["attributes"]["vhr.kind"] == "task"]
    assert sorted(s["name"] for s in tasks) == ["batch", "copy", "copy", "copy / overviews", "copy / overviews", "files"]
    batch = next(s for s in tasks if s["name"] == "batch")
    assert batch["attributes"]["vhr.scene_units"] == 2
    assert all(s["attributes"]["vhr.duration_seconds"] > 0 for s in tasks)


def test_analysis_never_overwrites_raw_data_and_identifies_malformed_line(tmp_path):
    raw, alias, output = tmp_path / "history.jsonl", tmp_path / "alias.jsonl", tmp_path / "summary.json"
    raw.write_text("")
    assert summarize_statistics(raw, output)["data"] == []
    os.link(raw, alias)
    for destination in (raw, alias):
        with pytest.raises(ValueError, match="differ"):
            summarize_statistics(raw, destination)
    raw.write_text('{"partial":')
    with StatisticsRecorder(raw) as recorder:
        recorder(event(kind="workflow", name="workflow"))
    assert raw.read_text().startswith('{"partial":\n{')
    before = output.read_bytes()
    with pytest.raises(ValueError, match=r"history.jsonl:1:"):
        summarize_statistics(raw, output)
    assert output.read_bytes() == before


def _append_from_process(arguments):
    raw, run = arguments
    with StatisticsRecorder(raw) as recorder:
        for _ in range(15):
            recorder(event(run=run))
        recorder(event(run=run, kind="workflow", name="workflow"))


def test_parallel_append_keeps_whole_records(tmp_path):
    raw = tmp_path / "history.jsonl"
    with ProcessPoolExecutor(max_workers=2) as executor:
        list(executor.map(_append_from_process, [(raw, str(i)) for i in range(4)]))
    records = spans(raw)
    assert len(records) == 64
    assert len({s["attributes"]["vhr.run_id"] for s in records}) == 4
    data = summarize_statistics(raw, tmp_path / "summary.json")["data"]
    assert next(r for r in data if r["kind"] == "task")["samples"] == 60


@pytest.mark.parametrize("direction", ["horizontal", "vertical"])
def test_dask_short_tasks_are_not_lost_between_refreshes(tmp_path, direction):
    distributed = pytest.importorskip("distributed")
    config = quiet_recipe(tmp_path, direction=direction)
    events = []
    with distributed.LocalCluster(n_workers=2, threads_per_worker=1, processes=False,
                                  dashboard_address=None) as cluster:
        config["shared"].update({"core:concurrent_processing_backend": "dask",
                                  "core:dask_scheduler": ["address", cluster.scheduler_address]})
        run_workflow(config, config_dir=str(tmp_path), event_callback=events.append)
    tasks = [e for e in events if e["kind"] == "task" and e["name"] != "files"]
    assert sorted(e["name"] for e in tasks) == ["deliver", "deliver", "prepare"]
    assert all(e["attributes"]["vhr.measured"] and e["status"] == "completed" for e in tasks)


@pytest.mark.parametrize("setting", ["save_statistics_path", "load_statistics_path"])
@pytest.mark.parametrize("value", ["", "  ", True, 12, "var:history", "expr:const.history"])
def test_statistics_paths_require_a_literal_filename(tmp_path, value, setting):
    config = quiet_recipe(tmp_path)
    config["shared"]["core:" + setting] = value
    with pytest.raises(ValueError, match=setting):
        Workflow(config, config_dir=tmp_path)


def test_defaults_save_and_load_beside_yaml_and_repeated_runs_refine_eta(tmp_path, monkeypatch):
    import yaml
    from vhrharmonize import load_workflow

    config = quiet_recipe(tmp_path)
    for setting in ("save_statistics_path", "load_statistics_path"):
        del config["shared"]["core:" + setting]
    config["shared"]["core:run_from_existing"] = False
    filename = tmp_path / "workflow.yml"
    filename.write_text(yaml.safe_dump(config, sort_keys=False))
    raw = tmp_path / "statistics.jsonl"
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    first = load_workflow(filename)
    updates = []
    first.run(progress_callback=updates.append)
    assert updates[0]["total"]["eta_seconds"] is None
    original = raw.read_bytes()
    assert not (elsewhere / "statistics.jsonl").exists()
    historical = load_statistics(raw)
    for key in ("prepare", "deliver"):
        assert historical[(key, "file_source", "process_pool")][1] == 3
    updates.clear()
    load_workflow(filename).run(progress_callback=updates.append)
    assert all(r["eta_seconds"] is not None for r in updates[0]["rows"])
    assert raw.read_bytes().startswith(original)
    assert all(load_statistics(raw)[(key, "file_source", "process_pool")][1] == 6
               for key in ("prepare", "deliver"))


def test_history_seeds_logs_and_snapshots_with_saving_and_ui_disabled(tmp_path, capsys):
    raw = tmp_path / "history.jsonl"
    with StatisticsRecorder(raw) as recorder:
        for name, seconds in (("prepare", 30), ("deliver", 20)):
            recorder(event(name=name, seconds=seconds, **{"vhr.plugin": "file_source", "vhr.scene_units": 1}))
        recorder(event(kind="workflow", name="workflow"))
    before = raw.read_bytes()
    config = quiet_recipe(tmp_path)
    config["shared"].update({"core:save_statistics_path": None, "core:log_to_console": True})
    snapshots = []
    Workflow(config, config_dir=tmp_path).run(progress_callback=snapshots.append)
    assert snapshots[0]["total"]["eta_seconds"] == 70
    assert [r["eta_seconds"] for r in snapshots[0]["rows"]] == [30, 40]
    output = capsys.readouterr().out
    assert "Estimated remaining runtime: ~70s" in output
    assert "prepare: loaded: 1 | processing: 1 | unused: 1 | ETA ~30s" in output
    assert "Start 1/1/3 | ETA ~30s" in output
    assert raw.read_bytes() == before


def test_history_mean_weights_aggregates_and_prefers_live_callback_eta(tmp_path, monkeypatch):
    from vhrharmonize.workflow.progress import StepProgress, TaskProgress, WorkflowProgress

    workflow = Workflow(quiet_recipe(tmp_path), config_dir=tmp_path)
    workflow._historical_timings = {("prepare", "file_source", "process_pool"): (30, 3)}
    reporter = WorkflowProgress(workflow.plan())
    state = reporter.state
    row = state.rows["prepare"]
    assert state.mean_seconds(row) == 10
    row.run = 5
    row.done = 1
    state.tasks["done"] = TaskProgress("prepare", "old", weight=1, duration=18, done=True)
    assert state.mean_seconds(row) == 12  # (30 + 18) / (3 + 1)
    monkeypatch.setattr("vhrharmonize.workflow.progress.monotonic", lambda: 100)
    state.tasks["active"] = TaskProgress("prepare", "new", weight=2, started=95)
    assert reporter.snapshot()["active"][0]["eta_seconds"] == 19
    state.tasks["active"].stats = dict(n=1, total=3, rate=1, elapsed=1)
    assert reporter.snapshot()["active"][0]["eta_seconds"] == 2
    assert state.remaining_work(row) == 43
    assert state.remaining_work(StepProgress("pending", run=5, pending=True, history=(30, 3))) is None


def test_history_ignores_failed_unmeasured_cached_and_duplicate_records(tmp_path):
    raw = tmp_path / "history.jsonl"
    assert load_statistics(raw) == {}
    with StatisticsRecorder(raw) as recorder:
        recorder(event(seconds=12, **{"vhr.scene_units": 3, "vhr.plugin": "alignment"}))
        recorder(event(seconds=600, status="failed"))
        recorder(event(seconds=800, **{"vhr.measured": False}))
        recorder(event(kind="step_summary", seconds=0, status="reused"))
        recorder(event(name="other", seconds=8, **{"vhr.plugin": "other", "vhr.backend": "dask"}))
    raw.write_text(raw.read_text() * 2)
    assert load_statistics(raw) == {("alignment", "alignment", "process_pool"): (12, 3),
                                    ("other", "other", "dask"): (8, 1)}


@pytest.mark.parametrize("location,value", [
    (("attributes", "vhr.schema_version"), 2),
    (("attributes", "vhr.duration_seconds"), float("nan")),
    (("attributes", "vhr.scene_units"), 0),
    (("attributes", "vhr.done"), -1),
    (("attributes", "vhr.plugin"), []),
    (("context", "trace_id"), "bad"),
    (("resource", "attributes"), []),
    (("end_time",), "not a timestamp"),
])
def test_validation_is_shared_by_history_and_summary(tmp_path, location, value):
    raw = tmp_path / "history.jsonl"
    with StatisticsRecorder(raw) as recorder:
        recorder(event())
    record = spans(raw)[0]
    assert validate_statistics_record(record) is record
    container = record
    for key in location[:-1]:
        container = container[key]
    container[location[-1]] = value
    with pytest.raises(ValueError, match="Malformed statistics"):
        validate_statistics_record(record)
    raw.write_text("\n" + json.dumps(record) + "\n")
    for consume in (load_statistics, lambda path: summarize_statistics(path, tmp_path / "summary.json")):
        with pytest.raises(ValueError, match=r"history.jsonl:2:"):
            consume(raw)


def test_malformed_history_fails_before_appending(tmp_path):
    config = quiet_recipe(tmp_path)
    raw = tmp_path / "history.jsonl"
    raw.write_text('{"broken":')
    with pytest.raises(ValueError, match=r"history.jsonl:1:"):
        Workflow(config, config_dir=tmp_path).run()
    assert raw.read_text() == '{"broken":'
