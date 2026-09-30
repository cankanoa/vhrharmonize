"""Scheduling order, synchronization and concurrent scene progress."""

from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path
from threading import Event

import pytest

from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.config import validate_config
from workflow_helpers import install_function, import_settings, copy_step, stage, transfer


def pipeline(monkeypatch, direction="horizontal"):
    events = []
    install_function(monkeypatch, "scenes", lambda: [{"name": "a"}, {"name": "b"}],
                     scene_records_return="$")

    def function(stage):
        def call(value):
            events.append((stage, value))
            return value
        return call

    recipe = {
        "defaults": {"plugin": "shared", "core:run": True,
                     "core:log_to_console": False, "core:processing_direction": direction},
        "source": {"plugin": "scenes", "core:run": True},
    }
    for index in range(1, 4):
        install_function(monkeypatch, f"stage{index}", function(index))
        recipe[f"stage{index}"] = {
            "plugin": f"stage{index}", "core:run": True,
            "param:value": "var:name" if index == 1 else "var:value",
            "var:value": "returned:$",
        }
    return recipe, events


@pytest.mark.parametrize("direction", ["horizontal", "vertical", None])
def test_direction_controls_call_order_and_preserves_values(monkeypatch, direction):
    recipe, events = pipeline(monkeypatch, direction)
    if direction is None:
        recipe["defaults"].pop("core:processing_direction")
    result = Workflow(recipe).run()
    expected = ([(stage, name) for stage in (1, 2, 3) for name in ("a", "b")]
                if direction == "horizontal" else
                [(stage, name) for name in ("a", "b") for stage in (1, 2, 3)])
    assert events == expected
    assert [r["context"]["var"]["value"] for r in result] == ["a", "b"]


def test_step_direction_overrides_shared_and_disabled_steps_do_not_block(monkeypatch):
    recipe, events = pipeline(monkeypatch)
    recipe["stage1"]["core:processing_direction"] = "vertical"
    recipe["stage2"]["core:processing_direction"] = "vertical"
    recipe = {
        key: value for key, value in [
            *list(recipe.items())[:3],
            ("disabled", {"core:run": False, "core:processing_direction": "horizontal"}),
            *list(recipe.items())[3:],
        ]
    }
    Workflow(recipe).run()
    assert events == [(1, "a"), (2, "a"), (1, "b"), (2, "b"), (3, "a"), (3, "b")]


def test_horizontal_override_synchronizes_vertical_scenes(monkeypatch):
    recipe, events = pipeline(monkeypatch, "vertical")
    recipe["stage2"]["core:processing_direction"] = "horizontal"
    Workflow(recipe).run()
    assert events == [(stage, name) for stage in (1, 2, 3) for name in ("a", "b")]


@pytest.mark.parametrize("aggregate", [False, True])
def test_collect_waits_for_all_preceding_scene_values(monkeypatch, aggregate):
    recipe, events = pipeline(monkeypatch, "vertical")

    def gather(values):
        assert values == ["a", "b"]
        events.append(("gather", tuple(values)))

    install_function(monkeypatch, "gather", gather, scope="aggregate" if aggregate else "scene")
    recipe["stage2"] = {"plugin": "gather", "core:run": True, "param:values": "collect:value"}
    Workflow(recipe).run()
    assert events == [(1, "a"), (1, "b"), *[("gather", ("a", "b"))] * (1 if aggregate else 2),
                      (3, "a"), (3, "b")]


def test_shared_returned_constants_synchronize_before_consumption(monkeypatch):
    recipe, events = pipeline(monkeypatch, "vertical")

    def shared_value(value):
        events.append((1, value))
        return 7

    install_function(monkeypatch, "shared_value", shared_value)
    recipe["stage1"] = {"plugin": "shared_value", "core:run": True,
                        "param:value": "var:name", "const:common": "returned:$"}
    recipe["stage2"]["param:value"] = "const:common"
    Workflow(recipe).run()
    assert events[:2] == [(1, "a"), (1, "b")]
    assert events[2:] == [(2, 7), (3, 7), (2, 7), (3, 7)]


def test_fast_scene_advances_while_another_scene_is_still_running(monkeypatch):
    from vhrharmonize.workflow import engine

    recipe, events = pipeline(monkeypatch, "vertical")
    recipe["defaults"]["core:concurrent_processing"] = 2
    monkeypatch.setattr(
        engine, "ProcessPoolExecutor",
        lambda *, max_workers, mp_context: ThreadPoolExecutor(max_workers=max_workers),
    )
    advanced = Event()

    def first(value):
        if value == "b":
            assert advanced.wait(5), "scene a could not advance while scene b was running"
        events.append((1, value))
        return value

    def second(value):
        events.append((2, value))
        if value == "a":
            advanced.set()
        return value

    install_function(monkeypatch, "first", first)
    install_function(monkeypatch, "second", second)
    recipe["stage1"]["plugin"] = "first"
    recipe["stage2"]["plugin"] = "second"
    Workflow(recipe).run()
    assert events.index((2, "a")) < events.index((1, "b"))


def test_dask_prioritizes_downstream_work_and_resolves_scene_returns(monkeypatch):
    import sys
    from types import ModuleType

    recipe, events = pipeline(monkeypatch, "vertical")
    recipe["defaults"].update({"core:concurrent_processing_backend": "dask",
                              "core:show_progress": False,  # Scheduling-only client stub.
                              "core:dask_scheduler_address": "tcp://scheduler:8786"})
    priorities = []

    class Client:
        def __init__(self, address):
            assert address == "tcp://scheduler:8786"

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def submit(self, function, payload, *, pure, priority):
            assert pure is False
            priorities.append(priority)
            future = Future()
            future.set_result(function(payload))
            return future

        def cancel(self, futures):
            for future in futures:
                future.cancel()

    distributed = ModuleType("dask.distributed")
    distributed.Client = Client
    distributed.as_completed = lambda futures: iter(sorted(futures, key=lambda f: f.result()))
    monkeypatch.setitem(sys.modules, "dask.distributed", distributed)
    result = Workflow(recipe).run()
    assert [r["context"]["var"]["value"] for r in result] == ["a", "b"]
    assert events.index((3, "a")) < events.index((2, "b"))
    assert max(priorities) > min(priorities)


@pytest.mark.parametrize("workers", [1, 2])
def test_vertical_copy_chain_reuses_checkpoints_cleans_up_and_stages(tmp_path, workers):
    for name in ("a", "b"):
        source = tmp_path / "source" / f"{name}.txt"
        source.parent.mkdir(exist_ok=True)
        source.write_text(name)
    recipe = {
        "shared": {"plugin": "shared", "core:run": True, "core:log_to_console": False,
                   "core:processing_direction": "vertical", "core:concurrent_processing": workers},
        "files": import_settings(tmp_path / "source/*.txt", tmp_path),
        "first": {**copy_step("first", "mul", suffix="_first"), "var:copied": "returned:$"},
        "second": copy_step("second", "first", suffix="_second", folder="output_dir"),
    }
    Workflow(recipe).run()
    assert sorted(p.read_text() for p in (tmp_path / "output").glob("*.txt")) == ["a", "b"]
    assert not list((tmp_path / "temp").glob("*.txt"))
    resumed = Workflow(recipe)
    assert resumed.counts()["first"]["processing"] == 0
    resumed.run()
    # Cached-only vertical sections do not require a live Dask scheduler.
    cached = {**recipe, "shared": {**recipe["shared"], "core:concurrent_processing_backend": "dask"}}
    Workflow(cached).run()
    staged, uploads, downloads = stage(recipe, tmp_path)
    assert staged["shared"]["core:processing_direction"] == "vertical"
    transfer(uploads)
    Workflow(staged).run()
    assert sorted(Path(remote).read_text() for local, remote in downloads.items()
                  if local.endswith(".txt")) == ["a", "b"]


def test_cross_scene_file_dependencies_are_respected_without_collect(tmp_path):
    for name in ("a", "b"):
        source = tmp_path / "source" / f"{name}.txt"
        source.parent.mkdir(exist_ok=True)
        source.write_text(name)
    recipe = {
        "shared": {"plugin": "shared", "core:run": True, "core:processing_direction": "vertical"},
        "files": import_settings(tmp_path / "source/*.txt", tmp_path),
        "first": copy_step("first", "mul", suffix="_first"),
        "second": {
            **copy_step("second", "first", suffix="_second", folder="output_dir"),
            "param:input_path": "expr:const.temp_dir & '/' & (var.basename = 'a' ? 'b' : 'a') & '_first.txt'",
        },
    }
    records = Workflow(recipe).run()
    assert [Path(r["context"]["var"]["second"]).read_text() for r in records] == ["b", "a"]


@pytest.mark.parametrize("shared", [False, True])
@pytest.mark.parametrize("value", ["diagonal", None, [], True])
def test_invalid_direction_is_rejected(shared, value):
    settings = {"core:run": True, "core:processing_direction": value}
    if shared:
        settings["plugin"] = "shared"
    with pytest.raises(ValueError, match="processing_direction must be horizontal or vertical"):
        validate_config({"settings": settings})
