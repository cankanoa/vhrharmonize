"""Aggregate assignments distribute values to scenes without losing ordering."""

from pathlib import Path
import json
import shutil

import pytest

from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.values import evaluate_settings
from workflow_helpers import import_settings, install_function, stage, transfer


def test_ordered_collect_reads_before_and_after_assignment():
    records = [{"var": {"image": "a"}}, {"var": {"image": "b"}}]
    params, updates, context, _ = evaluate_settings(
        {"param:before": "collect:image", "var:image": "expr:'new-' & var.image",
         "param:after": "collect:image"},
        {"var": {}, "const": {}}, records=records, scene_ids=["A", "B"],
    )
    assert params == {"before": ["a", "b"], "after": ["new-a", "new-b"]}
    assert updates["var.image"] == ["new-a", "new-b"]
    assert records[0]["var"]["image"] == "a"


@pytest.mark.parametrize("value", [[], [1], [1, 2, 3], {"A": 1}, {"A": 1, "C": 2}, 7])
def test_invalid_aggregate_assignment_is_rejected(value):
    with pytest.raises(ValueError, match="var:image must map exactly 2 scenes"):
        evaluate_settings({"var:image": "returned:$"}, {"const": {}, "var": {}},
                          records=[{"var": {}}, {"var": {}}], scene_ids=["A", "B"], returned=value)


def test_returned_collect_cannot_supply_current_function_input():
    with pytest.raises(ValueError, match="cannot depend on this function's returned"):
        evaluate_settings({"var:image": "returned:$", "param:input_images": "collect:image"},
                          {"const": {}, "var": {}}, records=[{"var": {}}], scene_ids=["A"])


@pytest.mark.parametrize("template", [
    "var:image", "var:$", "expr:var.image", "expr:$.var.image",
    "expr:$lookup($, 'var').image", "expr:$lookup($, $join(['v', 'ar'])).image",
    "expr:$", {"nested": ["var:image"]},
])
@pytest.mark.parametrize("count", [0, 1, 2])
def test_aggregate_parameters_require_explicit_collection(template, count):
    with pytest.raises(ValueError, match="use collect:"):
        evaluate_settings({"param:images": template}, {"const": {}, "var": {}},
                          records=[{"var": {"image": "a"}} for _ in range(count)])


@pytest.mark.parametrize("value", [7, [1, 2, 3], {"nested": ["a", "b"]}])
def test_literal_assignment_is_copied_to_each_scene(value):
    params, _, _, _ = evaluate_settings(
        {"var:value": value, "param:values": "collect:value"}, {"const": {}, "var": {}},
        records=[{"var": {}}, {"var": {}}],
    )
    assert params["values"] == [value, value]


def test_nested_batch_returns_map_before_scene_expressions():
    records = [{"var": {"label": "a"}}, {"var": {"label": "b"}}]
    _, updates, _, _ = evaluate_settings(
        {"var:result": {"score": "returned:scores", "label": "var:label"},
         "var:adjusted": "expr:var.result.score + 1"},
        {"const": {}, "var": {}}, records=records, scene_ids=["A", "B"],
        returned={"scores": {"B": 20, "A": 10}},
    )
    assert updates["var.result"] == [{"score": 10, "label": "a"}, {"score": 20, "label": "b"}]
    assert updates["var.adjusted"] == [11, 21]


def test_collect_assignment_explicitly_stores_whole_list_in_each_scene():
    _, updates, _, _ = evaluate_settings(
        {"var:neighbors": "collect:name"}, {"const": {}, "var": {}},
        records=[{"var": {"name": "a"}}, {"var": {"name": "b"}}],
    )
    assert updates["var.neighbors"] == [["a", "b"], ["a", "b"]]


def test_constants_and_whole_context_work_before_scene_initialization():
    params, updates, _, _ = evaluate_settings(
        {"param:context": "expr:$", "const:snapshot": "expr:$"},
        {"const": {"gain": 2}}, aggregate=True,
    )
    assert params["context"] == updates["const.snapshot"] == {"const": {"gain": 2}}


@pytest.mark.parametrize("returned", [False, True])
@pytest.mark.parametrize("mapping", [False, True])
def test_scene_batch_scene_roundtrip_and_hpc(tmp_path, monkeypatch, returned, mapping):
    sources = []
    for name in ("a", "b"):
        source = tmp_path / "source" / f"{name}.txt"
        source.parent.mkdir(exist_ok=True)
        source.write_text(name)
        sources.append(str(source))
    outputs = [str(tmp_path / "temp" / f"{name}.txt") for name in ("a", "b")]
    calls = []

    def batch(input_images, output_images):
        calls.append(input_images)
        for source, destination in zip(input_images, output_images):
            shutil.copy2(source, destination)
        return dict(zip(reversed(sources), reversed(output_images))) if mapping else output_images

    install_function(monkeypatch, "batch", batch, scope="aggregate",
                     input_paths={"input_images"}, output_paths={"output_images"})
    assignment = {"param:output_images": outputs, "var:image": "returned:$"} if returned else {
        "var:image": "expr:const.temp_dir & '/' & var.basename & '.txt'",
        "param:output_images": "collect:image",
    }
    recipe = {
        "shared": {"plugin": "shared", "core:run": True, "core:log_to_console": False,
                   "core:delete_temp_steps_proactively": False},
        "files": {**import_settings(tmp_path / "source/*.txt", tmp_path), "var:image": "returned:file_path"},
        "batch": {"plugin": "batch", "core:run": True, "param:input_images": "collect:image", **assignment},
        "finish": {"plugin": "file_source", "core:run": True, "param:input_path": "var:image",
                   "var:image": "expr:const.output_dir & '/' & var.basename & '.txt'",
                   "param:output_path": "var:image"},
    }
    workflow = Workflow(recipe)
    workflow.run()
    assert calls == [sources]
    assert [Path(record["context"]["var"]["image"]).read_text() for record in workflow.records] == ["a", "b"]
    # Cached aggregate values must restore to the appropriate scene.
    for path in (tmp_path / "output").glob("*.txt"):
        path.unlink()
    calls.clear()
    Workflow(recipe).run()
    assert not calls
    # Known paths can be staged even when their assignment came from a return.
    for path in (tmp_path / "output").glob("*.txt"):
        path.unlink()
    staged, uploads, downloads = stage(recipe, tmp_path)
    transfer(uploads)
    remote = Workflow(staged)
    remote.run()
    assert not calls
    assert [Path(downloads[str(tmp_path / 'output' / f'{name}.txt')]).read_text() for name in ("a", "b")] == ["a", "b"]


def test_final_metadata_uses_aggregate_updated_scene_values(tmp_path):
    source = tmp_path / "input.txt"
    source.write_text("input")
    destination = tmp_path / "metadata.json"
    recipe = {
        "shared": {"plugin": "shared", "core:run": True, "core:log_to_console": False,
                   "core:output_metadata_path": str(destination)},
        "files": import_settings(source, tmp_path),
        "update": {"core:run": True, "core:scope": "aggregate", "var:score": 9},
    }
    Workflow(recipe).run()
    assert json.loads(destination.read_text())[0]["var"]["score"] == 9


def test_nested_assignments_preserve_each_scenes_other_fields():
    records = [{"var": {"settings": {"gain": 1, "label": "a"}}},
               {"var": {"settings": {"gain": 2, "label": "b"}}}]
    params, _, context, _ = evaluate_settings(
        {"var:settings.gain": "expr:var.settings.gain + 2", "param:settings": "collect:settings"},
        {"var": {}, "const": {}}, records=records, scene_ids=["a", "b"],
    )
    assert params["settings"] == [{"gain": 3, "label": "a"}, {"gain": 4, "label": "b"}]


def test_relative_batch_outputs_update_nested_scene_paths(tmp_path, monkeypatch):
    sources = []
    for name in ("a", "b"):
        source = tmp_path / "source" / f"{name}.txt"
        source.parent.mkdir(exist_ok=True)
        source.write_text(name)
        sources.append(str(source))

    def batch(input_images, output_images):
        for source, destination in zip(input_images, output_images):
            shutil.copy2(source, destination)

    install_function(monkeypatch, "batch", batch, scope="aggregate",
                     input_paths={"input_images"}, output_paths={"output_images"})
    seen = []
    install_function(monkeypatch, "inspect", lambda image: seen.append(image),
                     input_paths={"image"})
    recipe = {
        "files": {**import_settings(tmp_path / "source/*.txt", tmp_path),
                  "var:image": {"path": "returned:file_path", "label": "original"}},
        "batch": {"plugin": "batch", "core:run": True,
                  "param:input_images": "collect:image.path",
                  "var:image.path": "expr:var.basename & '.txt'", "param:output_images": "collect:image.path"},
        "inspect": {"plugin": "inspect", "core:run": True, "param:image": "var:image.path"},
    }
    workflow = Workflow(recipe)
    workflow.run()
    assert seen == [str(tmp_path / "output" / f"{name}.txt") for name in ("a", "b")]
    assert [Path(path).read_text() for path in seen] == ["a", "b"]
    assert all(record["context"]["var"]["image"]["label"] == "original"
               for record in workflow.records)


@pytest.mark.parametrize("invalid", [[1], {"wrong": 1, "keys": 2}, 1])
def test_returned_mapping_is_validated_after_invocation(monkeypatch, tmp_path, invalid):
    install_function(monkeypatch, "scenes", lambda: [{}, {}], scene_records_return="$")
    install_function(monkeypatch, "result", lambda: invalid, scope="aggregate")
    recipe = {"source": {"plugin": "scenes", "core:run": True},
              "result": {"plugin": "result", "core:run": True, "var:data": "returned:$"}}
    workflow = Workflow(recipe)
    with pytest.raises(ValueError, match="must map exactly 2 scenes"):
        workflow.run()


def test_scene_collect_can_write_same_named_shared_constant(monkeypatch):
    install_function(monkeypatch, "scenes", lambda: [{"gain": 2}, {"gain": 3}], scene_records_return="$")
    seen = []
    install_function(monkeypatch, "inspect", lambda gain, all_gains: seen.append((gain, all_gains)))
    recipe = {
        "source": {"plugin": "scenes", "core:run": True},
        "inspect": {"plugin": "inspect", "core:run": True, "const:gain": "collect:gain",
                    "param:gain": "var:gain", "param:all_gains": "collect:gain"},
    }
    workflow = Workflow(recipe)
    workflow.run()
    assert workflow.context["const"]["gain"] == [2, 3]
    assert seen == [(2, [2, 3]), (3, [2, 3])]
    assert [r["context"]["var"]["gain"] for r in workflow.records] == [2, 3]


@pytest.mark.parametrize("assignment", ["var:gain", "returned:$", "expr:$", "expr:$lookup($, 'var').gain"])
def test_scene_constants_reject_conflicting_scalar_values(monkeypatch, assignment):
    install_function(monkeypatch, "scenes", lambda: [{"gain": 2}, {"gain": 3}], scene_records_return="$")
    install_function(monkeypatch, "inspect", lambda gain: gain)
    recipe = {
        "source": {"plugin": "scenes", "core:run": True},
        "inspect": {"plugin": "inspect", "core:run": True, "param:gain": "var:gain",
                    "const:shared_gain": assignment},
    }
    with pytest.raises(ValueError, match="conflicting scene values"):
        Workflow(recipe).run()
