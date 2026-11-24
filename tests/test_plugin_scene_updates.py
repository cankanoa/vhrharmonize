"""Processing plugins can split or merge scenes without discovery-side effects."""
from pathlib import Path

import pytest

from vhrharmonize.plugins.base import PluginResult
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.values import remap_paths
from workflow_helpers import import_settings, install_function


def recipe(tmp_path):
    for name in ("a", "b"):
        (tmp_path / f"{name}.txt").write_text(name)
    return {
        "shared": {"plugin": "shared", "core:run": True,
                   "core:show_progress": False, "core:log_to_console": False,
                   "core:save_statistics_path": None, "core:load_statistics_path": None,
                   "core:output_metadata_path": None},
        "files": {**import_settings(tmp_path / "*.txt", tmp_path),
                  "var:common": "same", "var:sensor": "var:basename"},
    }


def test_merge_preserves_provenance_and_runs_only_during_execution(tmp_path, monkeypatch):
    calls = []

    def merge():
        calls.append("merge")
        return PluginResult("normal-return", [{"source_indices": [0, 1], "scene_id": "joined",
                                               "values": {"basename": "joined", "answer": 42}}])

    plugin = install_function(monkeypatch, "merge_scenes_test", merge, scope="aggregate")
    plugin.updates_scenes = lambda **kwargs: True
    config = recipe(tmp_path)
    config["merge"] = {"plugin": "merge_scenes_test", "core:run": True,
                       "const:normal": "returned:$"}
    config["next"] = {"core:run": True, "var:after": "expr:var.answer + 1"}
    workflow = Workflow(config).plan()
    assert calls == []
    assert workflow.barrier_index is not None
    workflow.run()
    assert calls == ["merge"]
    assert len(workflow.records) == 1
    record = workflow.records[0]
    v = record["context"]["var"]
    assert record["id"] == "joined"
    assert v["common"] == "same"
    assert v["after"] == 43
    assert "sensor" not in v
    assert {r["sensor"] for r in v["parent_records"]} == {"a", "b"}
    assert len(v["parent_scene_ids"]) == 2
    assert set(record["source_paths"]) == {str(tmp_path / "a.txt"), str(tmp_path / "b.txt")}
    assert workflow.context["const"]["normal"] == "normal-return"


@pytest.mark.parametrize("updates, error", [
    ([{"source_indices": [4], "scene_id": "x", "values": {}}], "source_indices"),
    ([{"source_indices": [0], "scene_id": "x", "values": {}}] * 2, "unique"),
    ([{"source_indices": [True], "scene_id": "x", "values": {}}], "source_indices"),
    ([{"source_indices": [0], "scene_id": "x", "values": []}], "dictionary"),
])
def test_invalid_scene_updates_fail_before_replacing_records(tmp_path, monkeypatch, updates, error):
    plugin = install_function(monkeypatch, "bad_scenes", lambda: PluginResult(None, updates), scope="aggregate")
    plugin.updates_scenes = lambda **kwargs: True
    config = recipe(tmp_path)
    config["bad"] = {"plugin": "bad_scenes", "core:run": True}
    workflow = Workflow(config)
    original = [r["id"] for r in workflow.records]
    with pytest.raises(ValueError, match=error):
        workflow.run()
    assert [r["id"] for r in workflow.records] == original


@pytest.mark.parametrize("updates, remaining", [(None, 2), ([], 0)])
def test_no_change_and_empty_scene_update(tmp_path, monkeypatch, updates, remaining):
    plugin = install_function(monkeypatch, "scene_result", lambda: PluginResult(None, updates), scope="aggregate")
    plugin.updates_scenes = lambda **kwargs: True
    config = recipe(tmp_path)
    config["result"] = {"plugin": "scene_result", "core:run": True}
    workflow = Workflow(config)
    workflow.run()
    assert len(workflow.records) == remaining


def test_result_paths_remapped_in_values_and_scene_updates():
    value = PluginResult(["/scratch/output.tif"], [{"source_indices": [0], "scene_id": "output",
                          "values": {"file_path": "/scratch/output.tif"}}])
    mapped = remap_paths(value, {"/scratch/output.tif": "/final/output.tif"})
    assert mapped.value == ["/final/output.tif"]
    assert mapped.scene_updates[0]["values"]["file_path"] == "/final/output.tif"
