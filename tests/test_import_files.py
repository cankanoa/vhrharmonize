"""Public importer contract: explicit metadata rules and minimal scene dictionaries."""

import json

import pytest

from vhrharmonize.plugins.import_files import import_file, import_files
from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import import_settings, copy_step


def test_minimal_return_without_metadata_and_explicit_directory_arguments(tmp_path):
    source = tmp_path / "scene.v1.TIF"
    source.write_text("image")
    expected = {"scene_id": str(source), "file_path": str(source), "source_paths": [str(source)]}
    assert import_file(str(source)) == expected
    assert import_files(str(source), temp_dir="work", output_dir="products") == {
        "const": {
            "temp_dir": str(tmp_path / "work"),
            "output_dir": str(tmp_path / "products"),
        },
        "scenes": [expected],
    }
    assert list(tmp_path.iterdir()) == [source]


def test_rules_decode_metadata_and_resolve_relative_absolute_and_home_patterns(
    tmp_path, monkeypatch
):
    source = tmp_path / "scene.v1.TIF"
    source.write_text("image")
    document = source.with_suffix(".json")
    document.write_text(json.dumps({"companions": "*.{RPB,rpb}", "cloud_cover": 0.1}))
    sidecars = [tmp_path / "scene.v1.RPB", tmp_path / "scene.v2.rpb"]
    for sidecar in sidecars:
        sidecar.write_text("rpc")
    (tmp_path / "directory.RPB").mkdir()  # Only files are returned.
    monkeypatch.setenv("HOME", str(tmp_path))
    rules = {
        "metadata": {"to_json": "expr:$replace(var.file_path, '.TIF', '.json')"},
        "companions": {"path": "var:metadata.companions"},
        "absolute": {"path": str(sidecars[0])},
        "home": {"path": "~/scene.v2.rpb"},
        "multiple": {"path": ["*.RPB", "*.rpb", "*.RPB"]},
        "absent": {"path": "*.missing"},
    }
    scene = import_file(str(source), create_metadata_json=rules)
    assert set(scene) == {"scene_id", "file_path", "source_paths", *rules}
    assert scene["metadata"] == {"companions": "*.{RPB,rpb}", "cloud_cover": 0.1}
    assert scene["companions"] == scene["multiple"] == list(map(str, sidecars))
    assert scene["absolute"] == [str(sidecars[0])]
    assert scene["home"] == [str(sidecars[1])]
    assert scene["absent"] == []
    assert scene["source_paths"] == list(map(str, [source, document, *sidecars]))
    returned = import_files(
        str(source),
        create_metadata_json=rules,
        temp_dir="work",
        where="expr:var.metadata.cloud_cover < 0.5",
    )
    assert returned["scenes"] == [scene]
    assert (
        import_files(
            str(source),
            create_metadata_json=rules,
            temp_dir="work",
            where="expr:var.metadata.cloud_cover < 0.05",
        )["scenes"]
        == []
    )


@pytest.mark.parametrize(
    "rules",
    [
        [],
        {"metadata": "scene.IMD"},
        {"metadata": {}},
        {"metadata": {"path": "*.IMD", "to_json": "scene.IMD"}},
        {"metadata": {"unknown": "scene.IMD"}},
        {"": {"path": "*"}},
        {3: {"path": "*"}},
        *[
            {name: {"path": "*"}}
            for name in ("scene_id", "file_path", "source_paths", "temp_dir", "output_dir")
        ],
    ],
)
def test_invalid_rules_fail_even_when_no_files_match(tmp_path, rules):
    with pytest.raises(ValueError, match="create_metadata_json"):
        import_files(str(tmp_path / "*.TIF"), create_metadata_json=rules, temp_dir="work")


def test_metadata_missing_invalid_or_multiple_paths_raise_in_python_api(tmp_path):
    source = tmp_path / "scene.TIF"
    source.touch()
    with pytest.raises(FileNotFoundError):
        import_file(str(source), create_metadata_json={"metadata": {"to_json": "missing.IMD"}})
    invalid = tmp_path / "invalid.json"
    invalid.write_text("not json")
    with pytest.raises(ValueError):
        import_file(str(source), create_metadata_json={"metadata": {"to_json": str(invalid)}})
    with pytest.raises(ValueError, match="one path"):
        import_file(str(source), create_metadata_json={"metadata": {"to_json": [str(invalid)]}})


@pytest.mark.parametrize(
    "name", ["metadata_import", "companions", "constants", "metadata_output", "append_to_name"]
)
def test_removed_arguments_are_rejected_by_python_and_workflow(tmp_path, name):
    source = tmp_path / "scene.TIF"
    source.touch()
    with pytest.raises(TypeError, match="unexpected keyword"):
        import_files(str(source), **{name: "unused"})
    recipe = {"import": {**import_settings(source, tmp_path), f"param:{name}": "unused"}}
    with pytest.raises(ValueError, match=name):
        Workflow(recipe)


@pytest.mark.parametrize("kind", ["path", "to_json"])
def test_both_rule_types_protect_sources_from_output_overwrite(tmp_path, kind):
    source = tmp_path / "scene.TIF"
    source.write_text("image")
    companion = tmp_path / "metadata.json"
    companion.write_text('{"gain": 2}')
    recipe = {
        "shared": {"plugin": "shared", "core:run": True, "core:run_from_existing": False},
        "import": {
            **import_settings(source, tmp_path),
            "param:create_metadata_json": {"metadata": {kind: str(companion)}},
        },
        "copy": copy_step("product", "mul", str(companion)),
    }
    with pytest.raises(ValueError, match="protected input"):
        Workflow(recipe).run()
    assert companion.read_text() == '{"gain": 2}'


def test_scene_roots_resolve_per_file_without_derived_return_fields(tmp_path):
    for name in ("first", "second"):
        (tmp_path / name).mkdir()
        (tmp_path / name / "scene.TIF").touch()
    result = import_files(
        str(tmp_path / "*/scene.TIF"),
        directory_scope="var",
        temp_dir="expr:'./work/' & $split(var.file_path, '/')[-1]",
        output_dir=".",
    )
    assert result["const"] == {}
    for name, scene in zip(("first", "second"), result["scenes"]):
        assert set(scene) == {"scene_id", "file_path", "source_paths", "temp_dir", "output_dir"}
        assert scene["output_dir"] == str(tmp_path / name)
        assert scene["temp_dir"] == str(tmp_path / name / "work/scene.TIF")
