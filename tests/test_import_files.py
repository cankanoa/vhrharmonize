"""Public importer contract: explicit metadata rules and minimal scene dictionaries."""

import json

import pytest

from vhrharmonize.plugins.import_files import import_file, import_files
from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import import_settings, copy_step


def test_minimal_return_has_no_directory_defaults(tmp_path):
    source = tmp_path / "scene.v1.TIF"
    source.write_text("image")
    expected = {"scene_id": str(source), "file_path": str(source), "source_paths": [str(source)]}
    assert import_file(str(source)) == expected
    assert import_files(str(source)) == {"scenes": [expected]}
    assert import_files(str(tmp_path / "missing*.tif")) == {"scenes": []}
    assert list(tmp_path.iterdir()) == [source]


@pytest.mark.parametrize("temp_scope", ["const", "var"])
@pytest.mark.parametrize("output_scope", ["const", "var"])
def test_directory_scopes_are_ordinary_yaml_assignments(tmp_path, temp_scope, output_scope):
    for name in ("a", "b"):
        source = tmp_path / name / "image.tif"
        source.parent.mkdir()
        source.touch()
    settings = {"plugin": "import_files", "core:run": True,
                "param:search_glob": str(tmp_path / "*/image.tif")}
    for field, folder, scope in (("work", "work", temp_scope), ("products", "products", output_scope)):
        settings[f"{scope}:{field}"] = (
            "path:./" + folder if scope == "const"
            else "path:expr:$replace(var.file_path, /[^\\/]+$/, '') & '" + folder + "'"
        )
    workflow = Workflow({"discover": settings}, config_dir=tmp_path)
    for field, folder, scope in (("work", "work", temp_scope), ("products", "products", output_scope)):
        if scope == "const":
            assert workflow.initial_context["const"][field] == str(tmp_path / folder)
        else:
            assert [r["context"]["var"][field] for r in workflow.records] == [
                str(tmp_path / name / folder) for name in ("a", "b")
            ]


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
        where="expr:var.metadata.cloud_cover < 0.5",
    )
    assert returned["scenes"] == [scene]
    assert (
        import_files(
            str(source),
            create_metadata_json=rules,
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
            for name in ("scene_id", "file_path", "source_paths")
        ],
    ],
)
def test_invalid_rules_fail_even_when_no_files_match(tmp_path, rules):
    with pytest.raises(ValueError, match="create_metadata_json"):
        import_files(str(tmp_path / "*.TIF"), create_metadata_json=rules)


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
    "name", ["metadata_import", "companions", "constants", "metadata_output", "append_to_name", "directory_scope",
             "temp_dir", "output_dir", "temp_dir_scope", "output_dir_scope"]
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
        "copy": copy_step("product", "mul", str(companion), require_outputs=True),
    }
    with pytest.raises(ValueError, match="protected input"):
        Workflow(recipe).run()
    assert companion.read_text() == '{"gain": 2}'


def test_directory_names_are_available_for_ordinary_metadata(tmp_path):
    source = tmp_path / "scene.TIF"
    source.touch()
    companion = tmp_path / "scene.RPB"
    companion.touch()
    result = import_files(str(source), create_metadata_json={
        "temp_dir": {"path": "*.RPB"}, "output_dir": {"path": "*.RPB"},
    })
    assert result["scenes"][0]["temp_dir"] == [str(companion)]
    assert result["scenes"][0]["output_dir"] == [str(companion)]
