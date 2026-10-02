"""Directory roots are YAML values; cleanup requires an explicit core selection."""

from pathlib import Path
import shutil

import pytest

from vhrharmonize.workflow.config import validate_config
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.staging import stage_workflow


def recipe(tmp_path):
    for name in ("a", "b"):
        source = tmp_path / "inputs" / name / "source.txt"
        source.parent.mkdir(parents=True)
        source.write_text(name)
    return {
        "shared": {
            "plugin": "shared", "core:run": True, "core:log_to_console": False,
            "core:show_progress": False, "core:save_statistics_path": None,
            "core:load_statistics_path": None,
            "const:root": "path:" + str(tmp_path),
            "const:products": "path:./products", "const:cache": "path:./cache",
        },
        "discover": {
            "plugin": "import_files", "core:run": True,
            "param:search_glob": str(tmp_path / "inputs/*/source.txt"),
            "var:name": "expr:$split(var.file_path, '/')[-2]",
            "var:work": "path:expr:const.cache & '/' & var.name",
            "var:relative_output_dir": r"path:expr:$replace(var.file_path, /[^\/]+$/, '') & '../processed/' & var.name",
        },
        "prepare": {
            "plugin": "file_source", "core:run": True,
            "param:input_path": "var:file_path",
            "var:intermediate": "expr:var.work & '/intermediate.txt'",
            "param:output_path": "var:intermediate",
        },
        "finish": {
            "plugin": "file_source", "core:run": True,
            "core:require_outputs": "param:output_path",
            "param:input_path": "var:intermediate",
            "param:output_path": "expr:const.products & '/' & var.name & '.txt'",
        },
    }


@pytest.mark.parametrize("roots", [[], "const:cache", ["const:cache"], ["var:work"]])
@pytest.mark.parametrize("required", [False, True])
def test_cleanup_uses_explicit_roots_and_preserves_requested_outputs(tmp_path, roots, required):
    config = recipe(tmp_path)
    config["shared"]["core:cleanup_dirs"] = roots
    config["prepare"]["core:require_outputs"] = "param:output_path" if required else False
    unrelated = tmp_path / "cache/unrelated.txt"
    unrelated.parent.mkdir()
    unrelated.write_text("keep")
    workflow = Workflow(config, config_dir=tmp_path)
    assert workflow.initial_context["const"]["products"] == str(tmp_path / "products")
    assert workflow.directory_locations["output_dir"] == []
    workflow.run()
    for name in ("a", "b"):
        assert (tmp_path / "products" / f"{name}.txt").read_text() == name
        assert (tmp_path / "inputs" / name / "source.txt").read_text() == name
        assert (tmp_path / "cache" / name / "intermediate.txt").exists() == (required or not roots)
    assert unrelated.read_text() == "keep"


def test_hpc_rebases_yaml_roots_and_context_without_importer_directory_arguments(tmp_path):
    config = recipe(tmp_path)
    config["shared"]["core:cleanup_dirs"] = ["const:cache"]
    context_file = "path:expr:const.cache & '/context.json'"
    config["discover"].update({"core:save_context": {context_file: "all"},
                              "core:load_context": {context_file: "all"}})
    Workflow(config, config_dir=tmp_path)  # Explicit discovery context save.
    remote = tmp_path / "remote"
    staged, uploads, downloads = stage_workflow(
        config, config_dir=tmp_path, remote_work_dir=str(remote),
        path_mappings={"const:root": str(remote / "inputs"),
                       "const:cache": str(remote / "work"),
                       "const:products": str(remote / "output"),
                       "var:relative_output_dir": "expr:'" + str(remote / "scenes") + "/' & var.name"},
        context_staging_dir=tmp_path / "staged-context",
    )
    assert not {"param:temp_dir", "param:output_dir", "param:temp_dir_scope", "param:output_dir_scope"} & staged["discover"].keys()
    for local, destination in uploads.items():
        Path(destination).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(local, destination)
    shutil.rmtree(tmp_path / "inputs")
    resumed = Workflow(staged, config_dir=tmp_path)
    resumed.run()
    for record in resumed.records:
        variables = record["context"]["var"]
        name = variables["name"]
        assert variables["relative_output_dir"] == str(remote / "scenes" / name)
        assert (remote / "output" / f"{name}.txt").read_text() == name
        assert not (remote / "work" / name / "intermediate.txt").exists()
        assert downloads[str(tmp_path / "products" / f"{name}.txt")] == str(remote / "output" / f"{name}.txt")


@pytest.mark.parametrize("value", [None, True, "temp_dir", "path:./temp", "param:temp", "const:", ["var:bad-name"], [12], {}])
def test_cleanup_selectors_are_validated(value):
    with pytest.raises(ValueError, match="core:cleanup_dirs"):
        validate_config({"shared": {"plugin": "shared", "core:run": True,
                                    "core:cleanup_dirs": value}})
