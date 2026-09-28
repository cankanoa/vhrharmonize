"""Named steps, shared settings and steps without an implementation."""

import pytest
import yaml
from vhrharmonize.workflow.config import load_config, validate_config
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.api import run_plugin
from workflow_helpers import install_function, stage


def test_names_select_instances_and_plugin_selects_implementation(monkeypatch):
    calls = []
    install_function(monkeypatch, "record", lambda value: calls.append(value))
    config = {
        "defaults": {
            "plugin": "shared",
            "core:run": True,
            "core:log_to_console": False,
            "param:value": 4,
        },
        "first run": {"plugin": "record", "core:run": True},
        "second run": {"plugin": "record", "core:run": True, "param:value": 8},
    }
    counts = run_plugin("record", config)
    assert calls == [4, 8]
    assert list(counts) == ["first run", "second run"]


def test_context_steps_work_before_and_after_scene_discovery(monkeypatch):
    install_function(
        monkeypatch, "discover", lambda: [{"name": "a"}, {"name": "b"}], scene_records_return="$"
    )
    calls = []
    install_function(monkeypatch, "record", lambda value: calls.append(value))
    recipe = {
        "initial": {"core:run": True, "const:prefix": "start_"},
        "files": {"plugin": "discover", "core:run": True},
        "names": {"core:run": True, "var:label": "expr:const.prefix & var.name"},
        "setup": {"core:run": True, "core:scope": "aggregate", "const:labels": "collect:label"},
        "record": {"plugin": "record", "core:run": True, "param:value": "const:labels"},
    }
    workflow = Workflow(recipe).plan()
    assert len(workflow.records) == 2 and not calls
    workflow.run()
    assert calls == [["start_a", "start_b"], ["start_a", "start_b"]]


@pytest.mark.parametrize(
    "settings",
    [
        {"param:value": 1},
        {"const:value": "returned:$"},
        {"plugin": ["file_source", "alignment"]},
        {"plugin": "missing"},
    ],
)
def test_invalid_step_contracts_fail(settings):
    with pytest.raises(ValueError):
        validate_config({"arbitrary": settings})


def test_no_implicit_plugin_selection_and_no_list_steps():
    assert Workflow({"file_source": {"core:run": True, "const:x": 1}}).run() == []
    with pytest.raises(ValueError, match="without a plugin"):
        Workflow({"file_source": {"core:run": True, "param:input_path": "input"}})
    with pytest.raises(ValueError, match="one mapping"):
        validate_config({"one": [{"plugin": "file_source"}]})


@pytest.mark.parametrize(
    "text",
    [
        "one: {}\none: {}\n",
        "one:\n  plugin: file_source\n  plugin: alignment\n",
    ],
)
def test_duplicate_step_names_and_plugin_selectors_fail(tmp_path, text):
    filename = tmp_path / "recipe.yml"
    filename.write_text(text)
    with pytest.raises(ValueError, match="Duplicate YAML key"):
        load_config(filename)


def test_shared_is_identified_by_plugin_not_step_name(monkeypatch):
    calls = []
    install_function(monkeypatch, "record", lambda value=1: calls.append(value))
    Workflow(
        {
            "defaults": {"plugin": "shared", "core:run": True, "param:value": 9},
            "shared": {"plugin": "record", "core:run": True},
        }
    ).run()
    assert calls == [9]


def test_disabled_context_and_shared_steps_do_nothing(monkeypatch):
    calls = []
    install_function(monkeypatch, "record", lambda value=1: calls.append(value))
    workflow = Workflow(
        {
            "defaults": {"plugin": "shared", "param:value": 9},
            "setup": {"const:value": "var:missing"},
            "go": {"plugin": "record", "core:run": True},
        }
    )
    workflow.run()
    assert calls == [1] and workflow.context["const"] == {}


def test_named_steps_and_context_setup_survive_hpc(monkeypatch, tmp_path):
    calls = []
    install_function(monkeypatch, "discover", lambda: [{"n": 4}], scene_records_return="$")
    install_function(monkeypatch, "record", lambda n: calls.append(n))
    recipe = {
        "settings": {"plugin": "shared", "core:run": True, "core:log_to_console": False},
        "import anything": {"plugin": "discover", "core:run": True},
        "setup": {"core:run": True, "const:gain": 3},
        "shared": {"plugin": "record", "core:run": True, "param:n": "expr:var.n * const.gain"},
    }
    staged, _, _ = stage(recipe, tmp_path)
    assert staged["shared"]["plugin"] == "record" and staged["settings"]["plugin"] == "shared"
    assert staged["import anything"]["core:run"] is False
    Workflow(staged).run()
    assert calls == [12]


def test_core_controls_do_not_supply_function_parameters(monkeypatch):
    calls = []
    install_function(
        monkeypatch,
        "native",
        lambda concurrent_processing_backend=None, log_to_console=False: calls.append(
            (concurrent_processing_backend, log_to_console)
        ),
    )
    recipe = {
        "defaults": {"plugin": "shared", "core:run": True, "core:log_to_console": True},
        "call": {"plugin": "native", "core:run": True},
    }
    Workflow(recipe).run()
    assert calls == [(None, False)]
    recipe["defaults"]["param:concurrent_processing_backend"] = "process_pool"
    Workflow(recipe).run()
    assert calls[-1] == ("process_pool", False)


def test_multiple_shared_blocks_preserve_assignment_order_on_hpc(monkeypatch, tmp_path):
    calls = []
    install_function(monkeypatch, "native", lambda amount: calls.append(amount))
    recipe = {
        "defaults": {"plugin": "shared", "core:run": True, "const:amount": 2},
        "more_defaults": {
            "plugin": "shared",
            "core:run": True,
            "const:amount": "expr:const.amount + 3",
            "param:amount": "const:amount",
        },
        "call": {"plugin": "native", "core:run": True},
    }
    Workflow(recipe).run()
    staged, _, _ = stage(recipe, tmp_path)
    Workflow(staged).run()
    assert calls == [5, 5]
