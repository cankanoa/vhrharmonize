"""Scene discovery and replacement use the same ordinary function contract."""

from pathlib import Path
import pytest
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.api import run_plugin, run_workflow
from workflow_helpers import install_function, stage


def test_any_function_sets_plain_scenes_and_keeps_constants(monkeypatch, tmp_path):
    calls = []
    install_function(
        monkeypatch,
        "source",
        lambda: {
            "items": [{"name": "a", "nested": {"gain": 2}}, {"name": "b", "nested": {"gain": 3}}],
            "destination": str(tmp_path / "results"),
        },
        scene_records_return="items",
    )
    install_function(
        monkeypatch, "consume", lambda name, gain, common: calls.append((name, gain, common))
    )
    recipe = {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False, "const:common": 7},
        "source": {"plugin": 'source', 
            "core:run": True,
            "const:output_dir": "returned:destination",
            "var:label": "returned:name",
        },
        "consume": {"plugin": 'consume', 
            "core:run": True,
            "param:name": "var:label",
            "param:gain": "var:nested.gain",
            "param:common": "const:common",
        },
    }
    workflow = Workflow(recipe, config_dir=tmp_path)
    assert len(workflow.records) == 2
    assert workflow.counts()["consume"]["processing"] == 2
    workflow.run()
    assert calls == [("a", 2, 7), ("b", 3, 7)]
    assert workflow.context["const"]["output_dir"] == str(tmp_path / "results")
    assert workflow.records[0]["context"]["var"] == {
        "name": "a",
        "nested": {"gain": 2},
        "label": "a",
    }
    assert workflow.counts()["consume"]["processing"] == 2


def test_later_replacement_rebuilds_dependencies_and_scene_count(monkeypatch, tmp_path):
    seen = []
    install_function(
        monkeypatch,
        "seed",
        lambda: [{"n": 1}, {"n": 2}],
        scene_records_return="$",
    )
    install_function(monkeypatch, "double", lambda n: n * 2)
    install_function(
        monkeypatch, "reset", lambda values: [{"total": sum(values)}], scene_records_return="$"
    )
    install_function(monkeypatch, "consume", lambda total, label: seen.append((total, label)))
    recipe = {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False, "const:label": "unchanged"},
        "seed": {"plugin": 'seed', "core:run": True},
        "double": {"plugin": 'double', "core:run": True, "param:n": "var:n", "var:doubled": "returned:$"},
        "reset": {"plugin": 'reset', "core:run": True, "param:values": "collect:doubled"},
        "consume": {"plugin": 'consume', "core:run": True, "param:total": "var:total", "param:label": "const:label"},
    }
    workflow = Workflow(recipe, config_dir=tmp_path)
    workflow.run()
    assert seen == [(6, "unchanged")]
    assert len(workflow.records) == 1
    assert "n" not in workflow.records[0]["context"]["var"]
    assert workflow.counts()["double"]["processing"] == 2
    assert workflow.counts()["consume"]["processing"] == 1


def test_scene_replacement_checkpoint_restores_without_calling_function(monkeypatch, tmp_path):
    calls, seen = [], []

    def source(output_path):
        calls.append(output_path)
        Path(output_path).write_text("saved")
        return {"data": {"items": [{"n": 8}]}, "work": str(tmp_path / "work")}

    install_function(
        monkeypatch,
        "source",
        source,
        output_paths={"output_path"},
        scene_records_return="data.items",
        output_temporary_cleanup_paths=(),
    )
    install_function(monkeypatch, "consume", lambda n: seen.append(n))
    recipe = {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False},
        "source": {"plugin": 'source', 
            "core:run": True,
            "param:output_path": str(tmp_path / "manifest.txt"),
            "const:temp_dir": "returned:work",
        },
        "consume": {"plugin": 'consume', "core:run": True, "param:n": "var:n"},
    }
    for _ in range(2):
        Workflow(recipe, config_dir=tmp_path).run()
    assert len(calls) == 1
    assert seen == [8, 8]
    checkpoint = tmp_path / "manifest.txt.context.json"
    checkpoint.write_text('{"values": {}}')  # Missing scene snapshot cannot justify reuse.
    Workflow(recipe, config_dir=tmp_path).run()
    assert len(calls) == 2


@pytest.mark.parametrize("returned", [None, {}, [3], ["bad"]])
def test_invalid_scene_returns_fail_in_python_api(monkeypatch, tmp_path, returned):
    install_function(monkeypatch, "source", lambda: returned, scene_records_return="$")
    with pytest.raises(ValueError, match="list of plain dictionaries"):
        run_workflow({"source": {"plugin": 'source', "core:run": True}}, config_dir=tmp_path)


def test_empty_replacement_still_runs_aggregate_functions(monkeypatch, tmp_path):
    calls = []
    install_function(monkeypatch, "source", lambda: [], scene_records_return="$")
    install_function(monkeypatch, "scene", lambda: pytest.fail("No scenes should run"))
    install_function(
        monkeypatch, "aggregate", lambda values: calls.append(values), scope="aggregate"
    )
    workflow = Workflow(
        {
            "source": {"plugin": 'source', "core:run": True},
            "scene": {"plugin": 'scene', "core:run": True},
            "aggregate": {"plugin": 'aggregate', "core:run": True, "param:values": "collect:$"},
        },
        config_dir=tmp_path,
    )
    assert workflow.run() == []
    assert calls == [[]]


def test_dry_run_initializes_scenes_without_a_planning_opt_in(monkeypatch, tmp_path):
    calls = []
    install_function(
        monkeypatch,
        "source",
        lambda: calls.append("source") or [{"n": 1}],
        scene_records_return="$",
    )
    install_function(
        monkeypatch, "next", lambda n: pytest.fail("Processing does not run in dry-run")
    )
    recipe = {"source": {"plugin": 'source', "core:run": True}, "next": {"plugin": 'next', "core:run": True, "param:n": "var:n"}}
    counts = run_workflow(recipe, config_dir=tmp_path, dry_run=True)
    assert calls == ["source"]
    assert counts["next"]["processing"] == 1
    staged, _, _ = stage(recipe, tmp_path)
    assert staged["source"]["core:run"] is False


def test_hpc_snapshot_works_with_a_custom_discovery_plugin(monkeypatch, tmp_path):
    seen = []
    install_function(
        monkeypatch,
        "source",
        lambda: [{"n": 4}],
        scene_records_return="$",
    )
    install_function(monkeypatch, "consume", lambda n: seen.append(n))
    recipe = {"source": {"plugin": 'source', "core:run": True}, "consume": {"plugin": 'consume', "core:run": True, "param:n": "var:n"}}
    staged, uploads, downloads = stage(recipe, tmp_path)
    assert not uploads and not downloads
    assert staged["source"]["core:run"] is False
    Workflow(staged, config_dir=tmp_path).run()
    assert seen == [4]


def test_disabled_source_neither_runs_nor_resets_scenes(monkeypatch, tmp_path):
    install_function(
        monkeypatch,
        "seed",
        lambda: [{"n": 1}],
        scene_records_return="$",
    )
    install_function(
        monkeypatch, "disabled", lambda: pytest.fail("Disabled"), scene_records_return="$"
    )
    workflow = Workflow(
        {"seed": {"plugin": 'seed', "core:run": True}, "disabled": {"plugin": 'disabled', "core:run": False}}, config_dir=tmp_path
    )
    assert workflow.run()[0]["context"]["var"]["n"] == 1


def test_scene_reset_keeps_original_files_protected(monkeypatch, tmp_path):
    original = tmp_path / "original.txt"
    original.write_text("source")
    install_function(
        monkeypatch,
        "seed",
        lambda: [{"files": [str(original)]}],
        scene_records_return="$",
        source_file_protection_paths_return="files",
    )
    install_function(monkeypatch, "reset", lambda: [{}], scene_records_return="$")
    install_function(
        monkeypatch,
        "write",
        lambda output_path: Path(output_path).write_text("overwrite"),
        output_path_resolution_paths={"output_path"},
    )
    workflow = Workflow(
        {
            "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False},
            "seed": {"plugin": 'seed', "core:run": True},
            "reset": {"plugin": 'reset', "core:run": True},
            "write": {"plugin": 'write', "core:run": True, "param:output_path": str(original)},
        },
        config_dir=tmp_path,
    )
    with pytest.raises(ValueError, match="protected input"):
        workflow.run()
    assert original.read_text() == "source"


def test_scene_setters_require_aggregate_scope(monkeypatch, tmp_path):
    install_function(
        monkeypatch,
        "source",
        lambda: [{}],
        scene_records_return="$",
    )
    with pytest.raises(ValueError, match="aggregate scope"):
        Workflow({"source": {"plugin": 'source', "core:run": True, "core:scope": "scene"}}, config_dir=tmp_path)


def test_plugin_only_run_restores_scene_setter_or_reports_missing_outputs(monkeypatch, tmp_path):
    seen = []

    def source(output_path):
        Path(output_path).write_text("cached scenes")
        return [{"n": 5}]

    install_function(
        monkeypatch,
        "source",
        source,
        scene_records_return="$",
        output_paths={"output_path"},
        output_temporary_cleanup_paths=(),
    )
    install_function(monkeypatch, "consume", lambda n: seen.append(n))
    recipe = {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False},
        "source": {"plugin": 'source', "core:run": True, "param:output_path": str(tmp_path / "manifest.txt")},
        "consume": {"plugin": 'consume', "core:run": True, "param:n": "var:n"},
    }
    with pytest.raises(ValueError, match="Unselected step source"):
        run_plugin("consume", recipe, config_dir=tmp_path)
    run_workflow(recipe, config_dir=tmp_path)
    counts = run_plugin("consume", recipe, config_dir=tmp_path)
    assert seen == [5, 5]
    assert counts["source"]["loaded"] == 1


def test_planning_scene_replacement_preserves_independent_processing(monkeypatch, tmp_path):
    calls = []
    install_function(monkeypatch, "process", lambda: calls.append("process"))
    install_function(
        monkeypatch,
        "source",
        lambda: calls.append("source") or [{"name": "a"}],
        scene_records_return="$",
    )
    recipe = {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False},
        "process": {"plugin": 'process', "core:run": True},
        "source": {"plugin": 'source', "core:run": True},
    }
    workflow = Workflow(recipe, config_dir=tmp_path).plan()
    assert not calls
    workflow.run()
    assert calls == ["process", "source"]
    assert len(workflow.records) == 1
