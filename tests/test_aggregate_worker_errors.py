"""Aggregate restarts identify failed jobs and preserve ordinary reuse behavior."""

from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack, contextmanager
import json
from pathlib import Path

import pytest

from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import install_function


def copy_job(source, destination, fail_names, reuse):
    if Path(source).stem in fail_names:
        raise ValueError(f"cannot process {Path(source).stem}")
    output = Path(destination)
    if not (reuse and output.exists()):
        with output.open("a") as stream:
            stream.write(Path(source).stem + "\n")
    return str(output)


def aggregate_recipe(monkeypatch, tmp_path, *, policy="remove_and_restart", reuse=True,
                     core_reuse=True, fail_names=("bad",), require_return=True):
    attempts = []
    scenes = []
    for name in ("good", "bad", "other"):
        filename = tmp_path / f"{name}.txt"
        filename.write_text(name)
        scenes.append({"id": name, "image": str(filename)})
    install_function(monkeypatch, "scenes", lambda: scenes,
                     var_records_return="$", var_id_return="id")

    def aggregate(input_paths, output_paths, *, executor=None, reuse=True):
        attempts.append((tuple(Path(p).stem for p in input_paths), reuse))
        futures = [executor.submit(copy_job, source, output, fail_names, reuse)
                   for source, output in zip(input_paths, output_paths)]
        return [future.result() for future in futures]

    adapter = install_function(monkeypatch, "aggregate", aggregate, scope="aggregate",
                               input_paths={"input_paths"}, output_paths={"output_paths"},
                               executor_parameter="executor")
    adapter.worker_image_paths = lambda function, args, kwargs: args[0] if function is copy_job else ()
    config = {
        "shared": {"plugin": "shared", "core:run": True, "core:worker_error_case": policy,
                   "core:show_progress": False, "core:log_to_console": False,
                   "core:save_statistics_path": None, "core:load_statistics_path": None,
                   "core:run_from_existing": core_reuse, "const:root": str(tmp_path)},
        "source": {"plugin": "scenes", "core:run": True},
        "aggregate": {"plugin": "aggregate", "core:run": True,
                      "core:require_outputs": "param:output_paths",
                      "param:input_paths": "collect:image",
                      "var:product": "expr:const.root & '/' & var.id & '.out'",
                      "param:output_paths": "collect:product", "param:reuse": reuse},
    }
    if require_return:
        config["aggregate"]["const:products"] = "returned:$"
        install_function(monkeypatch, "consume", lambda products: products, scope="aggregate")
        config["consume"] = {"plugin": "consume", "core:run": True, "core:require_outputs": True,
                              "param:products": "const:products", "const:final_products": "returned:$"}
    return config, attempts, adapter


def thread_pool(monkeypatch):
    from vhrharmonize.workflow import engine

    @contextmanager
    def executor(workers, backend, scheduler=None):
        with ThreadPoolExecutor(max_workers=workers) as pool:
            yield pool

    monkeypatch.setattr(engine, "image_executor", executor)


@pytest.mark.parametrize("reuse", [True, False])
def test_restart_rebuilds_inputs_and_keeps_native_reuse_optional(monkeypatch, tmp_path, capsys, reuse):
    thread_pool(monkeypatch)
    config, attempts, _ = aggregate_recipe(monkeypatch, tmp_path, reuse=reuse)
    config["shared"]["core:output_metadata_path"] = str(tmp_path / "metadata.json")
    workflow = Workflow(config, config_dir=tmp_path)
    events, snapshots = [], []
    records = workflow.run(event_callback=events.append, progress_callback=snapshots.append)
    assert [record["id"] for record in records] == ["good", "other"]
    assert attempts == [(("good", "bad", "other"), reuse), (("good", "other"), reuse)]
    for name in ("good", "other"):
        assert (tmp_path / f"{name}.out").read_text().splitlines() == [name] * (1 if reuse else 2)
    assert workflow.context["const"]["final_products"] == [str(tmp_path / f"{name}.out") for name in ("good", "other")]
    assert len(json.loads((tmp_path / "metadata.json").read_text())) == 2
    assert [(e["step"], e["image_id"], e["error_type"]) for e in workflow.errors] == [("aggregate", "bad", "ValueError")]
    assert "copy_job" in workflow.errors[0]["traceback"]
    assert "Step: aggregate | Image ID: bad" in capsys.readouterr().out
    task_events = [event for event in events if event["kind"] == "task" and event["name"] == "aggregate"]
    assert sorted(event["status"] for event in task_events) == ["completed", "failed"]
    assert len({event["attributes"]["vhr.task_id"] for event in task_events}) == 2
    assert snapshots[-1]["active"] == []
    row = next(row for row in snapshots[-1]["rows"] if row["name"] == "aggregate")
    assert row["done"] == row["run"] == 2


@pytest.mark.parametrize("core_reuse", [True, False])
def test_restart_respects_core_reuse_without_forcing_reexecution(monkeypatch, tmp_path, core_reuse):
    thread_pool(monkeypatch)
    config, attempts, _ = aggregate_recipe(monkeypatch, tmp_path, core_reuse=core_reuse,
                                           reuse=False, require_return=False)
    workflow = Workflow(config, config_dir=tmp_path)
    assert [record["id"] for record in workflow.run()] == ["good", "other"]
    assert len(attempts) == (1 if core_reuse else 2)
    assert (tmp_path / "good.out").read_text().splitlines() == ["good"] * len(attempts)


@pytest.mark.parametrize("policy", ["error", "remove_instence"])
def test_existing_policies_still_raise_aggregate_worker_errors(monkeypatch, tmp_path, policy):
    thread_pool(monkeypatch)
    config, attempts, _ = aggregate_recipe(monkeypatch, tmp_path, policy=policy)
    workflow = Workflow(config, config_dir=tmp_path)
    with pytest.raises(ValueError, match="cannot process bad"):
        workflow.run()
    assert len(attempts) == 1
    assert workflow.errors[0]["image_id"] == "bad"


def test_multiple_failed_jobs_are_reported_and_removed_together(monkeypatch, tmp_path):
    thread_pool(monkeypatch)
    config, attempts, _ = aggregate_recipe(monkeypatch, tmp_path, fail_names=("bad", "other"))
    workflow = Workflow(config, config_dir=tmp_path)
    assert [record["id"] for record in workflow.run()] == ["good"]
    assert attempts[-1][0] == ("good",)
    assert {error["image_id"] for error in workflow.errors} == {"bad", "other"}


@pytest.mark.parametrize("wrapped_error", [False, True])
def test_repeated_restarts_keep_upstream_results_and_exception_identity(monkeypatch, tmp_path, wrapped_error):
    thread_pool(monkeypatch)
    config, attempts, adapter = aggregate_recipe(monkeypatch, tmp_path)
    prepared = []

    def prepare(input_path):
        prepared.append(Path(input_path).stem)
        return input_path

    install_function(monkeypatch, "prepare", prepare)
    config = {**dict(list(config.items())[:2]),
              "prepare": {"plugin": "prepare", "core:run": True, "core:require_outputs": True,
                          "param:input_path": "var:image", "var:image": "returned:$"},
              **dict(list(config.items())[2:])}

    def aggregate(input_paths, output_paths, *, executor=None, reuse=True):
        names = tuple(Path(filename).stem for filename in input_paths)
        attempts.append(names)
        fail_names = ("bad",) if "bad" in names else ("other",)
        futures = [executor.submit(copy_job, source, output, fail_names, reuse)
                   for source, output in zip(input_paths, output_paths)]
        try:
            return [future.result() for future in futures]
        except ValueError as exc:
            if wrapped_error:
                raise RuntimeError("aggregate worker failed") from exc
            raise

    adapter.function = lambda: aggregate
    workflow = Workflow(config, config_dir=tmp_path)
    events = []
    assert [record["id"] for record in workflow.run(event_callback=events.append)] == ["good"]
    assert prepared == ["good", "bad", "other"]
    assert attempts == [("good", "bad", "other"), ("good", "other"), ("good",)]
    assert [(error["image_id"], error["error_type"]) for error in workflow.errors] == [("bad", "ValueError"), ("other", "ValueError")]
    assert (tmp_path / "good.out").read_text().splitlines() == ["good"]
    tasks = [event for event in events if event["kind"] == "task" and event["name"] == "aggregate"]
    assert len({event["attributes"]["vhr.task_id"] for event in tasks}) == 3


def test_all_images_failing_terminates_without_an_empty_retry(monkeypatch, tmp_path):
    thread_pool(monkeypatch)
    config, attempts, _ = aggregate_recipe(monkeypatch, tmp_path, fail_names=("good", "bad", "other"))
    workflow = Workflow(config, config_dir=tmp_path)
    assert workflow.run() == []
    assert len(attempts) == 1
    assert len(workflow.errors) == 3


@pytest.mark.parametrize("ambiguous", [False, True])
def test_unknown_or_ambiguous_job_identity_stops_without_guessing(monkeypatch, tmp_path, ambiguous):
    thread_pool(monkeypatch)
    config, attempts, adapter = aggregate_recipe(monkeypatch, tmp_path)
    if ambiguous:
        adapter.worker_image_paths = lambda function, args, kwargs: [str(tmp_path / "good.txt"), str(tmp_path / "bad.txt")]
    else:
        adapter.worker_image_paths = lambda function, args, kwargs: ()
    workflow = Workflow(config, config_dir=tmp_path)
    with pytest.raises(ValueError, match="cannot process bad"):
        workflow.run()
    assert len(attempts) == 1
    assert workflow.errors[0]["image_id"] is None
    if ambiguous:
        assert set(workflow.errors[0]["image_ids"]) == {"good", "bad"}


def test_coordinator_error_is_not_blamed_on_an_unrelated_failed_job(monkeypatch, tmp_path):
    thread_pool(monkeypatch)
    config, _, adapter = aggregate_recipe(monkeypatch, tmp_path)

    def aggregate(input_paths, output_paths, *, executor=None, reuse=True):
        executor.submit(copy_job, str(tmp_path / "bad.txt"), str(tmp_path / "bad.out"), ("bad",), reuse)
        raise RuntimeError("coordinator failed")

    adapter.function = lambda: aggregate
    workflow = Workflow(config, config_dir=tmp_path)
    with pytest.raises(RuntimeError, match="coordinator failed"):
        workflow.run()
    assert len(workflow.errors) == 1
    assert workflow.errors[0]["image_id"] is None


@pytest.mark.parametrize("backend", ["process_pool", "dask"])
def test_restart_with_real_worker_backends(monkeypatch, tmp_path, backend):
    config, attempts, _ = aggregate_recipe(monkeypatch, tmp_path)
    config["shared"]["core:concurrent_processing"] = 2
    with ExitStack() as stack:
        if backend == "dask":
            distributed = pytest.importorskip("distributed")
            cluster = stack.enter_context(distributed.LocalCluster(
                n_workers=2, threads_per_worker=1, processes=False,
                dashboard_address=None, memory_limit=0,
            ))
            config["shared"].update({"core:concurrent_processing_backend": "dask",
                                     "core:concurrent_processing": 1,
                                     "core:dask_scheduler": ["address", cluster.scheduler_address]})
        workflow = Workflow(config, config_dir=tmp_path)
        snapshots = []
        assert [record["id"] for record in workflow.run(progress_callback=snapshots.append)] == ["good", "other"]
    assert len(attempts) == 2
    assert "copy_job" in workflow.errors[0]["traceback"]
    assert snapshots[-1]["active"] == []


def test_spectralmatch_job_adapter_recovers_a_real_aggregate(monkeypatch, tmp_path):
    pytest.importorskip("spectralmatch")
    from spectralmatch.utils_multiprocessing import _run_image_tasks
    from vhrharmonize.plugins.spectralmatch_base import SpectralMatchPlugin

    config, attempts, adapter = aggregate_recipe(monkeypatch, tmp_path)

    def aggregate(input_paths, output_paths, *, executor=None, reuse=True):
        attempts.append(tuple(input_paths))
        return _run_image_tasks(
            copy_job, [(source, output, ("bad",), reuse) for source, output in zip(input_paths, output_paths)],
            input_paths=input_paths, output_paths=output_paths, executor=executor,
        )

    adapter.function = lambda: aggregate
    adapter.worker_image_paths = SpectralMatchPlugin().worker_image_paths
    workflow = Workflow(config, config_dir=tmp_path)
    assert [record["id"] for record in workflow.run(progress_callback=lambda snapshot: None)] == ["good", "other"]
    assert len(attempts) == 2
    assert workflow.errors[0]["image_id"] == "bad"


def test_seamline_metadata_restarts_with_normal_record_reuse(monkeypatch, tmp_path, make_test_raster):
    pytest.importorskip("osgeo")
    import geopandas as gpd

    good = make_test_raster(tmp_path / "good.tif")
    bad = tmp_path / "bad.tif"
    bad.write_text("invalid raster")
    install_function(monkeypatch, "scenes", lambda: [{"id": name, "image": str(filename), "metadata": {}}
                     for name, filename in (("good", good), ("bad", bad))],
                     var_records_return="$", var_id_return="id")
    config = {
        "shared": {"plugin": "shared", "core:run": True, "core:worker_error_case": "remove_and_restart",
                   "core:show_progress": False, "core:log_to_console": False,
                   "core:save_statistics_path": None, "core:load_statistics_path": None},
        "source": {"plugin": "scenes", "core:run": True},
        "footprints": {"plugin": "seamline_metadata", "core:run": True, "core:require_outputs": True,
                       "param:image_paths": "collect:image", "param:metadata_records": "collect:metadata",
                       "param:footprint_source": "calculate_bounds", "param:epsg": 4326,
                       "param:output_path": str(tmp_path / "footprints.gpkg"), "param:reuse": True},
    }
    workflow = Workflow(config, config_dir=tmp_path)
    assert [record["id"] for record in workflow.run(progress_callback=lambda snapshot: None)] == ["good"]
    assert workflow.errors[0]["image_id"] == "bad"
    frame = gpd.read_file(tmp_path / "footprints.gpkg")
    assert list(frame["image_basename"]) == ["good.tif"]
