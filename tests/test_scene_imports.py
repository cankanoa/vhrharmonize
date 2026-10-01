"""Later imports add scenes and metadata while preserving completed workflow state."""

import json
from pathlib import Path

import pytest

from vhrharmonize.plugins.base import FunctionPlugin
from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import context_controls
from workflow_helpers import install_function, stage, transfer


def importing(pattern, tmp_path, **settings):
    return {
        "plugin": "import_files",
        "core:run": True,
        "param:search_glob": str(pattern),
        "param:temp_dir": str(tmp_path / "work"),
        "param:output_dir": str(tmp_path / "out"),
        **settings,
    }


def quiet():
    return {"plugin": "shared", "core:run": True, "core:log_to_console": False}


def test_later_import_merges_nested_metadata_and_preserves_constants(tmp_path):
    for name in ("a", "b"):
        (tmp_path / f"{name}.txt").write_text(name)
        (tmp_path / f"{name}.txt.extra.json").write_text(
            json.dumps(
                {
                    "gain": 99,
                    "nested": {"old": False, "new": True},
                    "bands": [9],
                    "empty": 4,
                }
            )
        )
    (tmp_path / "a.txt.json").write_text(
        json.dumps(
            {
                "gain": 2,
                "nested": {"old": True},
                "bands": [1, 2],
                "empty": None,
            }
        )
    )
    recipe = {
        "shared": quiet(),
        "initial": importing(
            tmp_path / "a.txt",
            tmp_path,
            **{
                "param:create_metadata_json": {
                    "metadata": {"to_json": "literal:expr:var.file_path & '.json'"},
                },
                "var:suffix": "_already_processed",
                "const:calibration": {"gain": 5},
            },
        ),
        "additional": importing(
            tmp_path / "*.txt",
            tmp_path,
            **{
                "param:temp_dir": str(tmp_path / "ignored_work"),
                "param:output_dir": str(tmp_path / "ignored_out"),
                "param:create_metadata_json": {
                    "metadata": {"to_json": "literal:expr:var.file_path & '.extra.json'"},
                },
                "var:suffix": "",
                "var:quality": "expr:var.suffix & '_ready'",
            },
        ),
    }
    workflow = Workflow(recipe)
    records = workflow.run()
    assert [Path(r["id"]).name for r in records] == ["a.txt", "b.txt"]
    first, second = [r["context"]["var"] for r in records]
    assert first["metadata"] == {
        "gain": 2,
        "nested": {"old": True, "new": True},
        "bands": [1, 2],
        "empty": None,
    }
    assert (
        first["suffix"] == "_already_processed" and first["quality"] == "_already_processed_ready"
    )
    assert second["metadata"]["gain"] == 99 and second["suffix"] == ""
    assert (
        first["source_paths"]
        == records[0]["source_paths"]
        == [
            str(tmp_path / "a.txt"),
            str(tmp_path / "a.txt.json"),
            str(tmp_path / "a.txt.extra.json"),
        ]
    )
    assert set(first["source_paths"]) <= workflow.protected_paths
    assert workflow.context["const"] == {
        "temp_dir": str(tmp_path / "work"),
        "calibration": {"gain": 5},
    }
    assert first["output_dir"] == str(tmp_path / "out")
    assert second["output_dir"] == str(tmp_path / "ignored_out")
    assert all(r["context"]["const"] == workflow.context["const"] for r in records)


@pytest.mark.parametrize("enabled", [True, False])
def test_empty_or_disabled_later_import_keeps_existing_scenes(tmp_path, enabled):
    source = tmp_path / "a.txt"
    source.write_text("source")
    recipe = {
        "shared": quiet(),
        "initial": importing(source, tmp_path, **{"var:value": 7}),
        "later": importing(
            tmp_path / "missing/*",
            tmp_path,
            **{
                "core:run": enabled,
                "var:value": 9,
            },
        ),
    }
    records = Workflow(recipe).run()
    assert len(records) == 1 and records[0]["context"]["var"]["value"] == 7


@pytest.mark.parametrize("direction", ["horizontal", "vertical"])
def test_later_import_preserves_runtime_and_cached_results_and_replans(tmp_path, monkeypatch, direction):
    calls, seen = [], []
    for name in ("a", "b", "c"):
        (tmp_path / f"{name}.txt").write_text(name)

    def process(input_path, output_path):
        calls.append(input_path)
        Path(output_path).write_text("processed")
        return {"score": 10 if Path(input_path).stem == "a" else 20}

    install_function(
        monkeypatch, "process", process, input_paths={"input_path"}, output_paths={"output_path"}
    )
    install_function(
        monkeypatch, "consume", lambda name, score, total: seen.append((name, score, total))
    )
    recipe = {
        "shared": quiet(),
        "initial": importing(tmp_path / "[ab].txt", tmp_path),
        "process": {
            "plugin": "process",
            "core:run": True, "core:require_outputs": True,
            "param:input_path": "var:file_path",
            "var:product": "expr:const.temp_dir & '/' & $split(var.file_path, '/')[-1]",
            "param:output_path": "var:product",
            "var:score": "returned:score",
        },
        "sum": {
            "plugin": "sum_scores",
            "core:run": True, "core:require_outputs": True,
            "param:iterable": "collect:score",
            "const:total": "returned:$",
        },
        "later": importing(tmp_path / "c.txt", tmp_path, **{"var:score": 0}),
        "consume": {
            "plugin": "consume",
            "core:run": True, "core:require_outputs": True,
            "param:name": "expr:$split(var.file_path, '/')[-1]",
            "param:score": "var:score",
            "param:total": "const:total",
        },
    }
    # Use a plain Python header rather than sum's positional-only argument.
    install_function(monkeypatch, "sum_scores", lambda iterable: sum(iterable), scope="aggregate")
    recipe["shared"]["core:processing_direction"] = direction
    recipe["process"].update(context_controls(tmp_path / "scores.json", "var.score"))
    first = Workflow(recipe)
    assert first.counts()["consume"]["pending"] == 1
    assert not calls
    first.run()
    assert first.counts()["process"]["processing"] == 2
    assert first.counts()["consume"]["processing"] == 3
    second = Workflow(recipe)
    second.run()
    assert len(calls) == 2
    assert second.counts()["process"]["loaded"] == 2
    assert seen == [("a.txt", 10, 30), ("b.txt", 20, 30), ("c.txt", 0, 30)] * 2
    assert all("product" in r["context"]["var"] for r in second.records[:2])


def test_consecutive_imports_stage_all_scenes_and_context_for_hpc(tmp_path, monkeypatch):
    for name in ("a", "b"):
        (tmp_path / f"{name}.txt").write_text(name)
    seen = []
    install_function(
        monkeypatch,
        "consume",
        lambda input_path, label: seen.append((Path(input_path).read_text(), label)),
        input_paths={"input_path"},
    )
    recipe = {
        "shared": quiet(),
        "initial": importing(tmp_path / "a.txt", tmp_path, **{"var:label": "original"}),
        "additional": importing(tmp_path / "*.txt", tmp_path, **{"var:label": "new"}),
        "consume": {
            "plugin": "consume",
            "core:run": True, "core:require_outputs": True,
            "param:input_path": "var:file_path",
            "param:label": "var:label",
        },
    }
    staged, uploads, _ = stage(recipe, tmp_path)
    assert staged["initial"]["core:run"] is staged["additional"]["core:run"] is True
    assert {str(tmp_path / "a.txt"), str(tmp_path / "b.txt")} <= uploads.keys()
    transfer(uploads)
    Workflow(staged).run()
    assert seen == [("a", "original"), ("b", "new")]


def test_var_directory_roots_survive_reimport_and_new_scenes_get_own_roots(tmp_path):
    for name in ("a", "b"):
        (tmp_path / name).mkdir()
        (tmp_path / name / "image.txt").touch()
    recipe = {
        "shared": quiet(),
        "initial": importing(
            tmp_path / "a/image.txt",
            tmp_path,
            **{
                "param:temp_dir_scope": "var",
                "param:output_dir_scope": "var",
                "param:temp_dir": "first",
                "param:output_dir": "products",
            },
        ),
        "additional": importing(
            tmp_path / "*/image.txt",
            tmp_path,
            **{
                "param:temp_dir_scope": "var",
                "param:output_dir_scope": "var",
                "param:temp_dir": "later",
                "param:output_dir": "new_products",
            },
        ),
    }
    first, second = Workflow(recipe).run()
    assert first["context"]["var"]["temp_dir"] == str(tmp_path / "a/first")
    assert first["context"]["var"]["output_dir"] == str(tmp_path / "a/products")
    assert second["context"]["var"]["temp_dir"] == str(tmp_path / "b/later")
    assert second["context"]["var"]["output_dir"] == str(tmp_path / "b/new_products")


def test_merge_declaration_requires_stable_id_and_valid_mode():
    plugin = FunctionPlugin()
    plugin.var_records_return = "scenes"
    plugin.var_records_mode = "merge"
    with pytest.raises(ValueError, match="var_id_return"):
        plugin.file_features()
    plugin.var_id_return = "name"
    plugin.file_features()
    plugin.var_records_mode = "invalid"
    with pytest.raises(ValueError, match="var_records_mode"):
        plugin.file_features()


def test_later_import_discovers_files_created_by_previous_processing(tmp_path, monkeypatch):
    source = tmp_path / "original.txt"
    source.write_text("source")
    generated = tmp_path / "work/generated.txt"
    seen = []
    install_function(
        monkeypatch,
        "consume",
        lambda input_path: seen.append(input_path),
        input_paths={"input_path"},
    )
    recipe = {
        "shared": quiet(),
        "initial": importing(source, tmp_path),
        "produce": {
            "plugin": "file_source",
            "core:run": True, "core:require_outputs": True,
            "param:input_path": "var:file_path",
            "param:output_path": str(generated),
            "var:processed_path": "returned:$",
        },
        "later": importing(tmp_path / "work/*.txt", tmp_path),
        "consume": {"plugin": "consume", "core:run": True, "core:require_outputs": True, "param:input_path": "var:file_path"},
    }
    workflow = Workflow(recipe)
    assert workflow.counts()["consume"]["pending"] == 1
    assert not generated.exists()
    workflow.run()
    assert seen == [str(source), str(generated)]
    assert workflow.records[0]["context"]["var"]["processed_path"] == str(generated)
    assert workflow.counts()["produce"]["processing"] == 1
