"""Worker failure policies, scene removal and end-of-run diagnostics."""

from concurrent.futures import Future, ThreadPoolExecutor
from contextlib import ExitStack
import json
from pathlib import Path
import sys
from types import ModuleType

import pytest

from vhrharmonize.workflow.config import shared_settings, validate_config
from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import copy_step, import_settings, install_function


def pipeline(monkeypatch, *, direction="vertical", policy="remove_instence", failed=("bad",)):
    calls = []
    names = ["bad", "good", "other"]
    install_function(monkeypatch, "scenes", lambda: [{"id": name} for name in names],
                     var_records_return="$", var_id_return="id")

    def first(value):
        calls.append(("first", value))
        if value in failed:
            raise ValueError(f"cannot process {value}")
        return value + " processed"

    def second(value):
        calls.append(("second", value))
        return value + " finished"

    def gather(values):
        calls.append(("gather", values))
        return values

    install_function(monkeypatch, "first", first)
    install_function(monkeypatch, "second", second)
    install_function(monkeypatch, "gather", gather, scope="aggregate")
    recipe = {
        "shared": {"plugin": "shared", "core:run": True, "core:log_to_console": False,
                   "core:show_progress": False, "core:save_statistics_path": None,
                   "core:load_statistics_path": None, "core:processing_direction": direction,
                   "core:worker_error_case": policy},
        "source": {"plugin": "scenes", "core:run": True},
        "first": {"plugin": "first", "core:run": True, "core:require_outputs": True,
                  "param:value": "var:id", "var:value": "returned:$"},
        "second": {"plugin": "second", "core:run": True, "core:require_outputs": True,
                   "param:value": "var:value", "var:value": "returned:$"},
        "gather": {"plugin": "gather", "core:run": True, "core:require_outputs": True,
                   "param:values": "collect:value", "const:results": "returned:$"},
    }
    return recipe, calls


def use_backend(monkeypatch, recipe, backend):
    if backend == "serial":
        return
    if backend == "pool":
        from vhrharmonize.workflow import engine

        recipe["shared"]["core:concurrent_processing"] = 2
        monkeypatch.setattr(engine, "ProcessPoolExecutor",
                            lambda *, max_workers, mp_context: ThreadPoolExecutor(max_workers=max_workers))
        return

    class Client:
        def __init__(self, address):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def submit(self, function, payload, **kwargs):
            future = Future()
            try:
                future.set_result(function(payload))
            except Exception as exc:
                future.set_exception(exc)
            return future

        def cancel(self, futures):
            for future in futures:
                future.cancel()

    distributed = ModuleType("dask.distributed")
    distributed.Client = Client
    distributed.as_completed = lambda futures: iter(futures)
    monkeypatch.setitem(sys.modules, "dask.distributed", distributed)
    recipe["shared"].update({"core:concurrent_processing_backend": "dask",
                             "core:dask_scheduler": ["address", "tcp://scheduler:8786"]})


@pytest.mark.parametrize("direction", ["horizontal", "vertical"])
@pytest.mark.parametrize("backend", ["serial", "pool", "dask"])
@pytest.mark.parametrize("policy", ["remove_instence", "remove_and_restart"])
def test_remove_failed_images_and_rebuild_collections(monkeypatch, capsys, direction, backend, policy):
    recipe, calls = pipeline(monkeypatch, direction=direction, policy=policy)
    use_backend(monkeypatch, recipe, backend)
    workflow = Workflow(recipe)
    records = workflow.run()
    assert [r["id"] for r in records] == ["good", "other"]
    expected = [name + " processed finished" for name in ("good", "other")]
    assert ("second", "bad processed") not in calls
    assert calls[-1] == ("gather", expected)
    assert workflow.context["const"]["results"] == expected
    assert [r["context"]["var"]["value"] for r in records] == expected
    assert len(workflow.errors) == 1
    error = workflow.errors[0]
    assert (error["image_id"], error["step"], error["error_type"]) == ("bad", "first", "ValueError")
    assert "raise ValueError" in error["traceback"]
    output = capsys.readouterr().out
    assert "Worker error report (1 error)" in output
    assert "Step: first | Image ID: bad" in output
    assert "ValueError: cannot process bad" in output
    assert "Traceback (most recent call last)" in output


@pytest.mark.parametrize("direction", ["horizontal", "vertical"])
def test_default_policy_raises_original_error(monkeypatch, capsys, direction):
    recipe, calls = pipeline(monkeypatch, direction=direction)
    recipe["shared"].pop("core:worker_error_case")
    assert shared_settings(recipe)[0]["worker_error_case"] == "error"
    workflow = Workflow(recipe)
    with pytest.raises(ValueError, match="cannot process bad"):
        workflow.run()
    assert calls == [("first", "bad")]
    assert "Step: first | Image ID: bad" in capsys.readouterr().out


@pytest.mark.parametrize("direction", ["horizontal", "vertical"])
@pytest.mark.parametrize("failed", [("bad", "other"), ("bad", "good", "other")])
def test_multiple_errors_and_no_survivors(monkeypatch, capsys, direction, failed):
    recipe, calls = pipeline(monkeypatch, direction=direction, failed=failed)
    workflow = Workflow(recipe)
    records = workflow.run()
    assert [r["id"] for r in records] == (["good"] if len(failed) == 2 else [])
    assert len(workflow.errors) == len(failed)
    assert {error["image_id"] for error in workflow.errors} == set(failed)
    if not records:
        assert not any(step == "gather" for step, _ in calls)
    assert f"Worker error report ({len(failed)} errors)" in capsys.readouterr().out


@pytest.mark.parametrize("direction", ["horizontal", "vertical"])
def test_failures_at_different_steps_preserve_surviving_context(monkeypatch, direction):
    recipe, calls = pipeline(monkeypatch, direction=direction)

    def second(value):
        if value.startswith("other"):
            raise RuntimeError("second step failed")
        return value + " finished"

    install_function(monkeypatch, "second_failure", second)
    recipe["second"]["plugin"] = "second_failure"
    workflow = Workflow(recipe)
    events = []
    records = workflow.run(event_callback=events.append)
    assert [r["id"] for r in records] == ["good"]
    assert records[0]["context"]["var"]["value"] == "good processed finished"
    assert calls[-1] == ("gather", ["good processed finished"])
    assert [(e["step"], e["image_id"]) for e in workflow.errors] == [("first", "bad"), ("second", "other")]
    outcome = next(e for e in events if e["kind"] == "workflow")
    assert outcome["status"] == "completed"
    assert outcome["attributes"]["vhr.worker_error_count"] == 2


def test_later_scene_import_does_not_reintroduce_removed_ids(monkeypatch):
    recipe, _ = pipeline(monkeypatch)
    install_function(monkeypatch, "more_scenes", lambda: [{"id": name} for name in ("bad", "good", "other")],
                     var_records_return="$", var_id_return="id", var_records_mode="merge")
    recipe = {**dict(list(recipe.items())[:3]),
              "more": {"plugin": "more_scenes", "core:run": True},
              **dict(list(recipe.items())[3:])}
    workflow = Workflow(recipe)
    assert [r["id"] for r in workflow.run()] == ["good", "other"]
    assert len(workflow.errors) == 1


@pytest.mark.parametrize("direction", ["horizontal", "vertical"])
def test_cross_scene_consumers_cannot_use_partial_failed_outputs(monkeypatch, tmp_path, direction):
    recipe, _ = pipeline(monkeypatch, direction=direction)
    consumed = []

    def write(value, output):
        Path(output).write_text(value)
        if value == "bad":
            raise ValueError("failed after writing a partial output")
        return output

    def read(input_path):
        consumed.append(Path(input_path).name)
        return Path(input_path).read_text()

    install_function(monkeypatch, "write", write, output_paths={"output"})
    install_function(monkeypatch, "read", read, input_paths={"input_path"})
    recipe["shared"]["const:root"] = str(tmp_path)
    recipe["first"].update({"plugin": "write", "param:output": "expr:const.root & '/' & var.id & '.txt'"})
    recipe["second"].pop("param:value")
    recipe["second"].update({"plugin": "read",
                             "param:input_path": "expr:const.root & '/' & (var.id = 'good' ? 'bad' : var.id) & '.txt'"})
    workflow = Workflow(recipe, config_dir=tmp_path)
    assert [r["id"] for r in workflow.run()] == ["other"]
    assert consumed == ["other.txt"]
    assert [(e["step"], e["image_id"]) for e in workflow.errors] == [("first", "bad"), ("second", "good")]


def test_success_does_not_print_an_error_report(monkeypatch, capsys):
    recipe, _ = pipeline(monkeypatch, failed=())
    workflow = Workflow(recipe)
    assert len(workflow.run()) == 3
    assert workflow.errors == []
    assert "Worker error report" not in capsys.readouterr().out


@pytest.mark.parametrize("direction", ["horizontal", "vertical"])
def test_parent_output_validation_removes_image_and_progress_closes(monkeypatch, tmp_path, direction):
    recipe, _ = pipeline(monkeypatch, direction=direction, failed=())

    def write(value, output):
        if value != "bad":
            Path(output).write_text(value)
        return value

    install_function(monkeypatch, "write", write, output_paths={"output"})
    recipe["first"].update({"plugin": "write", "param:output": "expr:const.root & '/' & var.id & '.txt'"})
    recipe["shared"]["const:root"] = str(tmp_path)
    snapshots = []
    workflow = Workflow(recipe, config_dir=tmp_path)
    records = workflow.run(progress_callback=snapshots.append)
    assert [r["id"] for r in records] == ["good", "other"]
    assert "did not produce valid declared outputs" in workflow.errors[0]["error"]
    assert snapshots[-1]["status"] == snapshots[-1]["total"]["status"] == "completed"
    assert snapshots[-1]["active"] == []
    rows = {row["name"]: row for row in snapshots[-1]["rows"]}
    assert rows["first"]["status"] == "completed"
    assert rows["first"]["done"] == rows["first"]["run"] == 2
    assert rows["first"]["unused"] == 1
    assert rows["gather"]["done"] == rows["gather"]["run"] == 2


def test_invalid_worker_outputs_do_not_publish_shared_values(monkeypatch, tmp_path):
    recipe, _ = pipeline(monkeypatch, failed=())

    def write(value, output):
        if value == "bad":
            return "invalid"
        Path(output).write_text(value)
        return "valid"

    install_function(monkeypatch, "write", write, output_paths={"output"})
    recipe["shared"]["const:root"] = str(tmp_path)
    recipe["first"].update({"plugin": "write", "const:common": "returned:$",
                            "param:output": "expr:const.root & '/' & var.id & '.txt'"})
    workflow = Workflow(recipe, config_dir=tmp_path)
    records = workflow.run()
    assert [r["id"] for r in records] == ["good", "other"]
    assert workflow.context["const"]["common"] == "valid"
    assert len(workflow.errors) == 1


def test_aggregate_error_is_reported_and_raised(monkeypatch, capsys):
    recipe, _ = pipeline(monkeypatch, failed=())

    def fail(values):
        raise RuntimeError("aggregate failed")

    install_function(monkeypatch, "fail", fail, scope="aggregate")
    recipe["gather"]["plugin"] = "fail"
    workflow = Workflow(recipe)
    with pytest.raises(RuntimeError, match="aggregate failed"):
        workflow.run()
    assert workflow.errors[0]["image_id"] is None
    assert "all scenes (aggregate)" in capsys.readouterr().out


@pytest.mark.parametrize("direction", ["horizontal", "vertical"])
@pytest.mark.parametrize("backend", ["process_pool", "dask"])
def test_real_workers_preserve_traceback_and_omit_failed_metadata(tmp_path, capsys, direction, backend):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    for name in ("bad", "good"):
        (source_dir / f"{name}.txt").write_text(name)
    companion_dir = tmp_path / "companion"
    companion_dir.mkdir()
    metadata = tmp_path / "metadata.json"
    recipe = {
        "shared": {"plugin": "shared", "core:run": True, "core:log_to_console": False,
                   "core:show_progress": False, "core:save_statistics_path": None,
                   "core:load_statistics_path": None, "core:concurrent_processing": 2,
                   "core:processing_direction": direction, "core:worker_error_case": "remove_instence",
                   "core:output_metadata_path": str(metadata), "const:companion": str(companion_dir)},
        "source": import_settings(source_dir / "*.txt", tmp_path),
        "first": {**copy_step("first", "mul", suffix="_first", require_outputs=True),
                  "param:companions": "expr:var.basename = 'bad' ? [const.companion] : []"},
        "second": copy_step("second", "first", suffix="_second", folder="output_dir", require_outputs=True),
    }
    with ExitStack() as stack:
        if backend == "dask":
            distributed = pytest.importorskip("distributed")
            cluster = stack.enter_context(distributed.LocalCluster(
                n_workers=2, threads_per_worker=1, processes=False,
                dashboard_address=None, memory_limit=0,
            ))
            recipe["shared"].update({"core:concurrent_processing_backend": "dask",
                                     "core:concurrent_processing": 1,
                                     "core:dask_scheduler": ["address", cluster.scheduler_address]})
        workflow = Workflow(recipe, config_dir=tmp_path)
        snapshots = []
        records = workflow.run(progress_callback=snapshots.append)
    assert len(records) == 1
    assert records[0]["context"]["var"]["basename"] == "good"
    assert not (tmp_path / "output/bad_first_second.txt").exists()
    assert (tmp_path / "output/good_first_second.txt").read_text() == "good"
    assert len(json.loads(metadata.read_text())) == 1
    assert workflow.errors[0]["image_id"] == str(source_dir / "bad.txt")
    assert "copy_file" in workflow.errors[0]["traceback"]
    assert "IsADirectoryError" in capsys.readouterr().out
    assert snapshots[-1]["active"] == []
    assert snapshots[-1]["status"] == snapshots[-1]["total"]["status"] == "completed"


@pytest.mark.parametrize("policy,failed,exit_code", [
    ("error", (), 0),
    ("remove_and_restart", ("bad",), 0),
    ("error", ("bad",), 1),
])
def test_done_logs_exit_code_for_success_recovery_and_failure(monkeypatch, capsys, policy, failed, exit_code):
    recipe, _ = pipeline(monkeypatch, policy=policy, failed=failed)
    recipe["shared"]["core:log_to_console"] = True
    workflow = Workflow(recipe)
    if exit_code:
        with pytest.raises(ValueError, match="cannot process bad"):
            workflow.run()
    else:
        workflow.run()
    output = capsys.readouterr().out
    assert output.splitlines()[-1] == f"[core:done] Exit code: {exit_code}"
    assert output.count("[core:done]") == 1


@pytest.mark.parametrize("value", [None, True, [], {}, "continue", "remove_instance"])
def test_invalid_policy_rejected(value):
    with pytest.raises(ValueError, match="worker_error_case must be"):
        validate_config({"shared": {"plugin": "shared", "core:run": True,
                                    "core:worker_error_case": value}})
