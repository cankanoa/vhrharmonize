"""Plugin scene removal happens during execution, before later collections run."""
from pathlib import Path

import pytest

from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import install_function


def recipe(monkeypatch, tmp_path, *, scope="aggregate", cached=False, remove_all=False):
    for name in ("a", "b"):
        (tmp_path / f"{name}.txt").write_text(name)
    calls, seen = [], []

    def select(input_images, output_path):
        calls.append(input_images)
        with open(output_path, "w") as output:
            output.write('"selected"')
        return "done"

    plugin = install_function(monkeypatch, "select_scenes", select, scope=scope,
                              input_paths={"input_images"}, output_paths={"output_path"})
    plugin.removes_scenes = True

    def remove(*, params, returned, scene_ids):
        assert returned in (None, "done")
        assert Path(params["output_path"]).exists()
        return [s for s in scene_ids if remove_all or s.endswith("a.txt")]

    plugin.scene_ids_to_remove = remove
    install_function(monkeypatch, "see_scenes", lambda input_images: seen.extend(input_images),
                     scope="aggregate", input_paths={"input_images"})
    output = str(tmp_path / "selection.json")
    if cached:
        Path(output).write_text('"selected"')
    config = {
        "shared": {"plugin": "shared", "core:run": True, "core:log_to_console": False,
                   "core:show_progress": False, "core:save_statistics_path": None,
                   "core:load_statistics_path": None},
        "input": {"plugin": "import_files", "core:run": True,
                  "param:search_glob": str(tmp_path / "*.txt"),
                  "var:image": "returned:file_path"},
        "select": {"plugin": "select_scenes", "core:run": True,
                   "param:input_images": "collect:image" if scope == "aggregate" else "expr:[var.image]",
                   "param:output_path": output if scope == "aggregate" else "expr:var.image & '.json'"},
        "see": {"plugin": "see_scenes", "core:run": True, "core:require_outputs": True,
                "param:input_images": "collect:image"},
    }
    return Workflow(config), calls, seen, plugin


@pytest.mark.parametrize("scope", ["aggregate", "var"])
@pytest.mark.parametrize("cached", [False, True])
def test_removal_at_execution_boundary(monkeypatch, tmp_path, scope, cached):
    if scope == "var" and cached:
        for name in ("a", "b"):
            (tmp_path / f"{name}.txt.json").write_text('"selected"')
    workflow, calls, seen, _ = recipe(monkeypatch, tmp_path, scope=scope, cached=cached)
    workflow.counts()
    assert len(workflow.records) == 2
    assert calls == []
    result = workflow.run()
    assert [r["id"] for r in result] == [str(tmp_path / "b.txt")]
    assert seen == [str(tmp_path / "b.txt")]
    assert bool(calls) is not cached
    assert workflow.errors == []
    assert (tmp_path / "a.txt").exists()


def test_remove_all_skips_downstream_aggregate(monkeypatch, tmp_path):
    workflow, _, seen, _ = recipe(monkeypatch, tmp_path, remove_all=True)
    assert workflow.run() == []
    assert seen == []


@pytest.mark.parametrize("invalid", [["unknown"], [1], ["duplicate", "duplicate"], None])
def test_invalid_scene_removal_fails(monkeypatch, tmp_path, invalid):
    workflow, _, seen, plugin = recipe(monkeypatch, tmp_path)
    plugin.scene_ids_to_remove = lambda **kwargs: invalid
    with pytest.raises(ValueError, match="unique scene IDs"):
        workflow.run()
    assert seen == []


def test_disabled_selection_does_not_remove_scenes(monkeypatch, tmp_path):
    workflow, calls, seen, plugin = recipe(monkeypatch, tmp_path)
    config = workflow.config
    # A disabled step must not invoke either its function or its hook.
    config['select']['core:run'] = False
    plugin.scene_ids_to_remove = lambda **kwargs: pytest.fail('disabled removal hook')
    result = Workflow(config).run()
    assert len(result) == 2
    assert len(seen) == 2
    assert calls == []


def test_removed_scenes_statistics_can_be_loaded(monkeypatch, tmp_path):
    from vhrharmonize.statistics import load_statistics

    workflow, _, _, _ = recipe(monkeypatch, tmp_path, scope="var")
    history = tmp_path / "history.jsonl"
    workflow.controls["save_statistics_path"] = str(history)
    workflow.run()
    assert history.exists()
    load_statistics(history)


def test_removal_plugin_without_output_files_still_executes(monkeypatch, tmp_path):
    from vhrharmonize.plugins.base import OUTPUT_PATH_FEATURES

    workflow, _, seen, plugin = recipe(monkeypatch, tmp_path)
    plugin.function = lambda: lambda input_images: input_images[:1]
    plugin.scene_ids_to_remove = lambda *, params, returned, scene_ids: returned
    for feature in OUTPUT_PATH_FEATURES:
        setattr(plugin, feature, frozenset())
    del workflow.config['select']['param:output_path']
    workflow = Workflow(workflow.config)
    assert len(workflow.records) == 2
    workflow.counts()
    assert len(workflow.records) == 2
    assert [r['id'] for r in workflow.run()] == [str(tmp_path / 'b.txt')]
    assert seen == [str(tmp_path / 'b.txt')]


def test_failed_step_does_not_call_removal_hook(monkeypatch, tmp_path):
    workflow, _, seen, plugin = recipe(monkeypatch, tmp_path)

    def fail(input_images, output_path):
        raise RuntimeError('selection failed')

    plugin.function = lambda: fail
    plugin.scene_ids_to_remove = lambda **kwargs: pytest.fail('hook ran after failure')
    with pytest.raises(RuntimeError, match='selection failed'):
        workflow.run()
    assert len(workflow.records) == 2
    assert seen == []
