import json
from pathlib import Path

import pytest

from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.context_io import select_fields
from vhrharmonize.workflow.config import validate_config
from workflow_helpers import install_function


def shared(**values):
    return {"plugin": "shared", "core:run": True, "core:log_to_console": False,
            "core:show_progress": False, "core:save_statistics_path": None,
            "core:load_statistics_path": None, **values}


def importer(source):
    return {"plugin": "import_files", "core:run": True, "param:search_glob": str(source),
            "param:scene_id": "literal:expr:$split(var.file_path, '/')[-1]",
            "var:raw": "returned:file_path"}


def test_save_is_explicit_selective_and_load_precedes_reuse(monkeypatch, tmp_path):
    calls, seen = [], []
    source = tmp_path / "scene.txt"
    source.write_text("source")
    output, snapshot = tmp_path / "result.txt", tmp_path / "metadata.json"

    def produce(input_path, output_path):
        calls.append(1)
        Path(output_path).write_text(Path(input_path).read_text())
        return {"gain": 2, "unused": "not selected"}

    install_function(monkeypatch, "producer", produce, input_paths=("input_path",), output_paths=("output_path",))
    install_function(monkeypatch, "consumer", lambda gain: seen.append(gain))
    config = {
        "settings": shared(), "import": importer(source),
        "produce": {"plugin": "producer", "core:run": True,
                    "param:input_path": "var:raw", "param:output_path": str(output),
                    "var:gain": "returned:gain", "var:unused": "returned:unused",
                    "core:load_context": {str(snapshot): "var.gain"},
                    "core:save_context": {str(snapshot): "var.gain"}},
        "consume": {"plugin": "consumer", "core:run": True, "core:require_outputs": True,
                    "param:gain": "var:gain"},
    }
    Workflow(config, config_dir=tmp_path).run()
    assert json.loads(snapshot.read_text()) == {"const": {}, "scenes": {"scene.txt": {"gain": 2}}}
    assert not Path(str(output) + ".context.json").exists()
    Workflow(config, config_dir=tmp_path).run()
    assert calls == [1]
    assert seen == [2, 2]


@pytest.mark.parametrize("logging", [False, True])
def test_missing_load_warns_only_when_enabled_and_does_not_run_producer(monkeypatch, tmp_path, capsys, logging):
    calls = []
    install_function(monkeypatch, "unused_producer", lambda: calls.append("producer"))
    install_function(monkeypatch, "use", lambda gain: calls.append(gain))
    config = {"settings": shared(**{"core:log_to_console": logging, "const:gain": 3}),
              "load": {"core:run": True, "core:load_context": {str(tmp_path / "missing.json"): "const.gain"}},
              "unused": {"plugin": "unused_producer", "core:run": True},
              "use": {"plugin": "use", "core:run": True, "core:require_outputs": True, "param:gain": "const:gain"}}
    Workflow(config, config_dir=tmp_path).run()
    output = capsys.readouterr()
    assert ("Context file not found" in output.out + output.err) == logging
    assert calls == [3]


@pytest.mark.parametrize("control,expected", [
    ("load_context", {"new": 2, "nested": {"b": 2}}),
    ("load_upsert_context", {"old": 1, "new": 2, "nested": {"a": 1, "b": 2}}),
])
def test_load_replaces_or_merges_selected_values(tmp_path, control, expected):
    filename = tmp_path / "context.json"
    filename.write_text(json.dumps({"const": {"metadata": {"new": 2, "nested": {"b": 2}}}}))
    config = {"settings": shared(**{"const:metadata": {"old": 1, "nested": {"a": 1}}, "const:keep": 8}),
              "load": {"core:run": True, "core:" + control: {str(filename): "const.metadata"}}}
    workflow = Workflow(config, config_dir=tmp_path)
    assert workflow.initial_context["const"] == {"metadata": expected, "keep": 8}


@pytest.mark.parametrize("control,expected", [
    ("save_context", {"const": {"metadata": {"new": 2}}, "scenes": {}}),
    ("save_upsert_context", {"const": {"keep": 7, "metadata": {"old": 1, "new": 2}}, "scenes": {}}),
])
def test_save_replaces_or_upserts_files(tmp_path, control, expected):
    filename = tmp_path / "context.json"
    filename.write_text(json.dumps({"const": {"keep": 7, "metadata": {"old": 1}}, "scenes": {}}))
    config = {"settings": shared(**{"const:metadata": {"new": 2}}),
              "save": {"core:run": True, "core:require_outputs": True,
                       "core:" + control: {str(filename): "const.metadata"}}}
    Workflow(config, config_dir=tmp_path).run()
    assert json.loads(filename.read_text()) == expected


def test_selectors_expand_without_assignment_history():
    settings = {"param:image": "var:image", "var:label": "ok"}
    definitions = {"var.image": "expr:const.root & '/' & var.basename", "const.root": "/data",
                   "var.basename": "scene", "var.label": "ok", "const.unused": 5}
    assert select_fields(["dependencies", "defined", "var.image"], settings, definitions) == [
        "const.root", "var.basename", "var.image", "var.label"]
    with pytest.raises(ValueError, match="Invalid context selector"):
        validate_config({"save": {"core:run": True, "core:save_context": {"file.json": "dependancies"}}})


def test_dependency_scope_selectors_expand_actual_fields_and_their_dependencies():
    definitions = {"var.gain": "const:gain", "const.gain": 4, "const.unused": 7}
    context = {"var": {"gain": 4, "metadata": {"angle": 2}}, "const": {"gain": 4, "unused": 7}}
    assert select_fields("dependencies", {"param:variables": "var:$"}, definitions, context) == [
        "const.gain", "var.gain", "var.metadata"]


def test_import_snapshots_store_constants_once_and_match_scenes(tmp_path):
    for name in ("a.txt", "b.txt"):
        (tmp_path / name).write_text(name)
    filename = tmp_path / "context.json"
    config = {"settings": shared(**{"const:bands": [1, 2]}),
              "import": {**importer(tmp_path / "*.txt"), "core:save_context": {str(filename): ["defined", "const.bands"]}}}
    Workflow(config, config_dir=tmp_path)
    data = json.loads(filename.read_text())
    assert data["const"] == {"bands": [1, 2]}
    assert set(data["scenes"]) == {"a.txt", "b.txt"}
    assert all(set(scene) == {"raw"} for scene in data["scenes"].values())
    loaded = Workflow({"settings": shared(), "load": {"core:run": True,
        "core:load_context": {str(filename): ["var.raw", "const.bands"]}}}, config_dir=tmp_path)
    assert [record["id"] for record in loaded.records] == ["a.txt", "b.txt"]
    assert loaded.records[1]["context"]["var"]["raw"] == str(tmp_path / "b.txt")


def test_satisfies_connects_cached_output_and_skips_missing_upstream_metadata(monkeypatch, tmp_path):
    source = tmp_path / "scene_cloudmasked.txt"
    source.write_text("cached")
    calls = []
    install_function(monkeypatch, "earlier", lambda amount, output: calls.append("earlier"), output_paths=("output",))
    install_function(monkeypatch, "mask", lambda image, output: calls.append("mask"), input_paths=("image",), output_paths=("output",))
    install_function(monkeypatch, "use_image", lambda image: calls.append(Path(image).read_text()), input_paths=("image",))
    config = {"settings": shared(), "import": {**importer(source), "core:satisfies": {"mask": "output"}},
              "earlier": {"plugin": "earlier", "core:run": True, "param:amount": "var:missing_metadata",
                          "var:current": str(tmp_path / "earlier.txt"), "param:output": "var:current"},
              "mask": {"plugin": "mask", "core:run": True, "param:image": "var:current",
                       "var:current": str(tmp_path / "mask.txt"), "param:output": "var:current"},
              "use": {"plugin": "use_image", "core:run": True, "core:require_outputs": True,
                      "param:image": "var:current"}}
    workflow = Workflow(config, config_dir=tmp_path)
    assert workflow.counts()["earlier"]["unused"] == 1
    workflow.run()
    assert calls == ["cached"]


def test_partial_scene_load_keeps_other_scenes_metadata_dependencies(monkeypatch, tmp_path):
    calls, seen = [], []
    for name in ("a", "b"):
        (tmp_path / f"{name}.txt").write_text(name)
        (tmp_path / f"{name}.out").write_text("cached")

    def measure(image, output):
        calls.append(Path(image).stem)
        return 9

    install_function(monkeypatch, "measure", measure, input_paths=("image",), output_paths=("output",))
    install_function(monkeypatch, "inspect", lambda score: seen.append(score))
    snapshot = tmp_path / "partial.json"
    snapshot.write_text(json.dumps({"scenes": {"a.txt": {"score": 4}}}))
    config = {"settings": shared(), "import": importer(tmp_path / "*.txt"),
              "measure": {"plugin": "measure", "core:run": True, "param:image": "var:raw",
                          "param:output": "expr:$replace(var.raw, '.txt', '.out')", "var:score": "returned:$"},
              "load": {"core:run": True, "core:load_context": {str(snapshot): "var.score"}},
              "inspect": {"plugin": "inspect", "core:run": True, "core:require_outputs": True,
                          "param:score": "var:score"}}
    Workflow(config, config_dir=tmp_path).run()
    assert calls == ["b"]
    assert seen == [4, 9]


def test_existing_invalid_context_or_missing_selected_field_is_an_error(tmp_path):
    snapshot = tmp_path / "invalid.json"
    config = {"settings": shared(), "load": {"core:run": True,
        "core:load_context": {str(snapshot): "const.metadata"}}}
    for invalid in ("{broken", '{"values":{}}', '{"const":{}}'):
        snapshot.write_text(invalid)
        with pytest.raises(ValueError):
            Workflow(config, config_dir=tmp_path)


def test_context_only_loader_publishes_new_constants_to_following_steps(tmp_path, monkeypatch):
    snapshot = tmp_path / "context.json"
    snapshot.write_text(json.dumps({"scenes": {"a": {"score": 2}}}))
    seen = []
    install_function(monkeypatch, "inspect", lambda score: seen.append(score))
    config = {"settings": shared(), "load": {"core:run": True,
              "core:load_context": {str(snapshot): "var.score"}, "const:gain": 4},
              "inspect": {"plugin": "inspect", "core:run": True, "core:require_outputs": True,
                          "param:score": "expr:var.score * const.gain"}}
    Workflow(config, config_dir=tmp_path).run()
    assert seen == [8]


def test_loaded_metadata_alone_satisfies_a_value_only_producer(tmp_path, monkeypatch):
    calls, seen = [], []
    snapshot = tmp_path / "calibration.json"
    install_function(monkeypatch, "calibrate", lambda: calls.append(1) or 3, scope="aggregate")
    install_function(monkeypatch, "inspect", lambda gain: seen.append(gain), scope="aggregate")
    config = {"settings": shared(),
              "calibrate": {"plugin": "calibrate", "core:run": True, "const:gain": "returned:$",
                            "core:load_context": {str(snapshot): "const.gain"},
                            "core:save_context": {str(snapshot): "const.gain"}},
              "inspect": {"plugin": "inspect", "core:run": True, "core:require_outputs": True,
                          "param:gain": "const:gain"}}
    Workflow(config, config_dir=tmp_path).run()
    Workflow(config, config_dir=tmp_path).run()
    assert calls == [1]
    assert seen == [3, 3]


@pytest.mark.parametrize("selection", ["const.metadata", "all"])
def test_requested_save_tracks_selected_metadata_dependencies(tmp_path, monkeypatch, selection):
    calls = []
    install_function(monkeypatch, "measure", lambda: calls.append(1) or {"gain": 5}, scope="aggregate")
    snapshot = tmp_path / "selected.json"
    config = {"settings": shared(),
              "measure": {"plugin": "measure", "core:run": True, "const:metadata": "returned:$"},
              "save": {"core:run": True, "core:require_outputs": True,
                       "core:save_context": {str(snapshot): selection}}}
    Workflow(config, config_dir=tmp_path).run()
    assert calls == [1]
    assert json.loads(snapshot.read_text())["const"] == {"metadata": {"gain": 5}}


def test_multiple_import_contexts_merge_new_scenes_without_searching_raw_files(tmp_path):
    config = {"settings": shared()}
    for name in ("a", "b"):
        snapshot = tmp_path / f"{name}.json"
        snapshot.write_text(json.dumps({"scenes": {name: {"raw": f"/missing/{name}.txt"}}}))
        config[name] = {**importer(tmp_path / "no-raw-files"),
                        "core:load_context": {str(snapshot): "var.raw"}}
    workflow = Workflow(config, config_dir=tmp_path)
    assert [record["id"] for record in workflow.records] == ["a", "b"]
    assert [record["context"]["var"]["raw"] for record in workflow.records] == ["/missing/a.txt", "/missing/b.txt"]
    assert workflow.discovery_sources == set()


@pytest.mark.parametrize("control", ["load_context", "load_upsert_context"])
@pytest.mark.parametrize("selection", ["all", ["all", "defined", "dependencies", "const.metadata.new"]])
def test_load_all_uses_saved_fields_per_scene_and_preserves_load_semantics(tmp_path, control, selection):
    seed = tmp_path / "seed.json"
    seed.write_text(json.dumps({"scenes": {
        "a": {"metadata": {"old": 1}, "local_only": 8}, "b": {"other": 0}
    }}))
    snapshot = tmp_path / "snapshot.json"
    snapshot.write_text(json.dumps({"const": {"metadata": {"new": 2}}, "scenes": {
        "a": {"metadata": {"new": 3}, "items": [1, 2]}, "b": {"other": 4}
    }}))
    workflow = Workflow({
        "settings": shared(**{"const:metadata": {"old": 1}, "const:outside": 9}),
        "seed": {"core:run": True, "core:load_context": {str(seed): "all"}},
        "load": {"core:run": True, "core:" + control: {str(snapshot): selection}},
    }, config_dir=tmp_path)
    upsert = control == "load_upsert_context"
    assert workflow.initial_context["const"] == {
        "metadata": {**({"old": 1} if upsert else {}), "new": 2}, "outside": 9
    }
    assert {record["id"]: record["context"]["var"] for record in workflow.records} == {
        "a": {"metadata": {**({"old": 1} if upsert else {}), "new": 3}, "local_only": 8, "items": [1, 2]},
        "b": {"other": 4},
    }


@pytest.mark.parametrize("control", ["save_context", "save_upsert_context"])
def test_save_all_includes_imported_fields_and_combines_scene_writes(tmp_path, control):
    snapshot = {"const": {"metadata": {"new": 1}}, "scenes": {
        "a": {"score": 3, "metadata": {"new": 4}}, "b": {"different_field": [5]}, "empty": {}
    }}
    seed = tmp_path / "seed.json"
    seed.write_text(json.dumps(snapshot))
    destination = tmp_path / "saved.json"
    destination.write_text(json.dumps({"const": {"metadata": {"old": 0}}, "scenes": {"previous": {"score": 0}}}))
    Workflow({
        "settings": shared(),
        "seed": {"core:run": True, "core:load_context": {str(seed): "all"}},
        "save": {"core:run": True, "core:require_outputs": True,
                 "core:" + control: {str(destination): ["all", "defined"]}},
    }, config_dir=tmp_path).run()
    if control == "save_upsert_context":
        snapshot["const"]["metadata"]["old"] = 0
        snapshot["scenes"]["previous"] = {"score": 0}
    assert json.loads(destination.read_text()) == snapshot


def test_all_round_trip_skips_raw_discovery(tmp_path):
    source = tmp_path / "a.txt"
    source.write_text("raw")
    snapshot = tmp_path / "import.json"
    config = {"settings": shared(**{"const:bands": [1, 2]}),
              "import": {**importer(source), "core:save_context": {str(snapshot): "all"},
                         "core:load_context": {str(snapshot): "all"}}}
    original = Workflow(config, config_dir=tmp_path)
    saved = json.loads(snapshot.read_text())
    assert saved == {"const": original.initial_context["const"], "scenes": {
        record["id"]: record["context"]["var"] for record in original.records
    }}
    source.unlink()
    restored = Workflow(config, config_dir=tmp_path)
    assert restored.discovery_sources == set()
    assert restored.records[0]["context"] == original.records[0]["context"]


def test_all_preserves_empty_scene_collection_without_rediscovery(tmp_path):
    snapshot = tmp_path / "empty.json"
    snapshot.write_text(json.dumps({"const": {}, "scenes": {}}))
    source = tmp_path / "unexpected.txt"
    source.write_text("must not be discovered")
    workflow = Workflow({"settings": shared(), "import": {
        **importer(source), "core:load_context": {str(snapshot): "all"}
    }}, config_dir=tmp_path)
    assert workflow.records == []
    assert workflow.discovery_sources == set()
