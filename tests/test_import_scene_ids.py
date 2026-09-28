"""Import IDs are function parameters, independent of file paths and plugin labels."""

from pathlib import Path

import pytest
import yaml

from vhrharmonize.cli.main import main
from vhrharmonize.plugins.import_files import import_file, import_files
from vhrharmonize.workflow.api import run_workflow
from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import install_function, stage, transfer


def settings(pattern, tmp_path, **extra):
    return {
        "plugin": "import_files",
        "core:run": True,
        "param:search_glob": str(pattern),
        "param:temp_dir": str(tmp_path / "temp"),
        "param:output_dir": str(tmp_path / "out"),
        **extra,
    }


def test_id_parameter_reads_parsed_metadata_before_filtering(tmp_path):
    source = tmp_path / "scene.txt"
    source.touch()
    source.with_suffix(".json").write_text('{"identity": {"id": "acquisition-17"}}')
    arguments = {
        "scene_id": "var:metadata.identity.id",
        "create_metadata_json": {"metadata": {"to_json": "scene.json"}},
    }
    scene = import_file(str(source), **arguments)
    assert scene["scene_id"] == "acquisition-17"
    assert import_files(
        str(source), **arguments, temp_dir="temp", where="expr:var.scene_id = 'acquisition-17'"
    )["scenes"] == [{**scene, "output_dir": str(tmp_path / "output")}]
    assert import_file(str(source), scene_id="my-explicit-id")["scene_id"] == "my-explicit-id"


@pytest.mark.parametrize("invalid", [None, "", "   ", 17, True, {}, []])
def test_invalid_id_fails_in_the_function_for_literals_and_expressions(tmp_path, invalid):
    source = tmp_path / "scene.txt"
    source.touch()
    for function in (import_file, import_files):
        with pytest.raises(ValueError, match="scene_id must resolve to a nonempty string"):
            function(str(source), scene_id=invalid)
    for expression in ("expr:null", "expr:17", "expr:[]", "expr:{}", "expr:''"):
        with pytest.raises(ValueError, match="scene_id must resolve to a nonempty string"):
            import_file(str(source), scene_id=expression)


def test_different_paths_merge_on_configured_id_without_losing_old_values(tmp_path):
    for name in ("mul", "pan"):
        folder = tmp_path / name
        folder.mkdir()
        (folder / "scene.txt").write_text(name)
    identifier = r"literal:expr:$replace($split(var.file_path, '/')[-1], /\.[^.]*$/, '')"
    recipe = {
        "mul": settings(
            tmp_path / "mul/scene.txt",
            tmp_path,
            **{
                "param:scene_id": identifier,
                "var:mul": "returned:file_path",
                "var:quality": 8,
            },
        ),
        "pan": settings(
            tmp_path / "pan/scene.txt",
            tmp_path,
            **{
                "param:scene_id": identifier,
                "var:pan": "returned:file_path",
                "var:quality": 2,
            },
        ),
    }
    (record,) = Workflow(recipe).run()
    assert record["id"] == record["context"]["var"]["scene_id"] == "scene"
    fields = record["context"]["var"]
    assert fields["mul"] == fields["file_path"] == str(tmp_path / "mul/scene.txt")
    assert fields["pan"] == str(tmp_path / "pan/scene.txt")
    assert fields["quality"] == 8
    assert fields["source_paths"] == record["source_paths"] == [fields["mul"], fields["pan"]]


def test_same_file_with_different_ids_creates_independent_scenes(tmp_path):
    source = tmp_path / "scene.txt"
    source.touch()
    recipe = {
        "first": settings(source, tmp_path, **{"param:scene_id": "first", "var:label": "A"}),
        "second": settings(source, tmp_path, **{"param:scene_id": "second", "var:label": "B"}),
    }
    records = Workflow(recipe).run()
    assert [r["id"] for r in records] == ["first", "second"]
    assert [r["context"]["var"]["label"] for r in records] == ["A", "B"]


def test_duplicate_ids_within_an_import_raise_the_same_python_and_cli_error(tmp_path):
    for name in ("a", "b"):
        (tmp_path / f"{name}.txt").touch()
    pattern = str(tmp_path / "*.txt")
    with pytest.raises(ValueError, match="unique within an import") as direct:
        import_files(pattern, scene_id="duplicate", temp_dir=str(tmp_path / "temp"))
    config = {"import": settings(pattern, tmp_path, **{"param:scene_id": "duplicate"})}
    filename = tmp_path / "recipe.yml"
    filename.write_text(yaml.safe_dump(config, sort_keys=False))
    with pytest.raises(ValueError) as python:
        run_workflow(filename, dry_run=True)
    with pytest.raises(ValueError) as cli:
        main(["workflow", "--config", str(filename), "--dry-run"])
    assert str(direct.value) == str(python.value) == str(cli.value)


@pytest.mark.parametrize("identifier", ["acquisition/17", "~/logical-id"])
def test_hpc_preserves_custom_ids_and_merged_input_aliases(tmp_path, monkeypatch, identifier):
    first = tmp_path / "mul.txt"
    second = tmp_path / "pan.txt"
    first.write_text("mul")
    second.write_text("pan")
    seen = []
    install_function(
        monkeypatch,
        "inspect_scene",
        lambda scene_id, mul, pan: seen.append(
            (scene_id, Path(mul).read_text(), Path(pan).read_text())
        ),
        input_paths={"mul", "pan"},
    )
    recipe = {
        "mul": settings(
            first, tmp_path, **{"param:scene_id": identifier, "var:mul": "returned:file_path"}
        ),
        "pan": settings(
            second, tmp_path, **{"param:scene_id": identifier, "var:pan": "returned:file_path"}
        ),
        "inspect": {
            "plugin": "inspect_scene",
            "core:run": True,
            "param:scene_id": "var:scene_id",
            "param:mul": "var:mul",
            "param:pan": "var:pan",
        },
    }
    staged, uploads, _ = stage(recipe, tmp_path)
    transfer(uploads)
    (record,) = Workflow(staged).run()
    assert record["id"] == record["context"]["var"]["scene_id"] == identifier
    assert seen == [(identifier, "mul", "pan")]
