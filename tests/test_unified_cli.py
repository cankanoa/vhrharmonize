from copy import deepcopy
from pathlib import Path
import shutil
import tempfile

import pytest
import yaml

from vhrharmonize.cli.main import main
from vhrharmonize.workflow.api import load_workflow, run_plugin, run_workflow
from vhrharmonize.workflow.config import validate_config
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow import registry
from workflow_helpers import import_settings, copy_step


@pytest.fixture
def recipe(tmp_path):
    source = tmp_path / "inputs" / "sample.bin"
    source.parent.mkdir()
    source.write_bytes(b"Any file can be imported, not only an image.")
    return {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False},
        "import_files": {"plugin": 'import_files', 
            **import_settings(source, tmp_path),
            "const:output_dir": "./products",
            "const:temp_dir": "./work",
        },
        "file_source": copy_step("copied", "mul", "expr:const.output_dir & '/' & var.filename"),
    }


def test_default_system_temp_is_shared_per_run_and_unique_between_runs(
    recipe, tmp_path, monkeypatch
):
    system = tmp_path / "system"
    system.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(system))
    (tmp_path / "inputs" / "second.bin").write_bytes(b"second")
    discovery = recipe["import_files"]
    discovery["param:search_glob"] = str(tmp_path / "inputs" / "*.bin")
    discovery.pop("const:output_dir")
    discovery.pop("const:temp_dir")
    first = load_workflow(recipe, config_dir=str(tmp_path))
    second = load_workflow(recipe, config_dir=str(tmp_path))
    records = [r["context"]["const"] for r in first.records]
    root = records[0]["temp_dir"]
    assert Path(root).is_dir()
    assert Path(root).parent == system
    assert {r["temp_dir"] for r in records} == {root}
    assert records[0]["output_dir"] == str(tmp_path / "inputs/output")
    assert second.records[0]["context"]["const"]["temp_dir"] != root
    assert "temp_dir" not in first.shared
    assert "temp" not in first.records[0]["context"]["var"]


def test_imports_any_file_and_resolves_companions_relative_to_it(recipe, tmp_path):
    (tmp_path / "inputs/sample.json").write_text('{"gain":2}')
    recipe["import_files"].update(
        {
            "param:create_metadata_json": {
                "source": {"to_json": "literal:expr:$replace(var.file_path, '.bin', '.json')"},
                "docs": {"path": "literal:expr:$replace(var.file_path, '.bin', '.json')"},
            },
            "var:gain": "var:source.gain",
        }
    )
    workflow = load_workflow(recipe, config_dir=str(tmp_path))
    metadata = {**workflow.records[0]["context"]["const"], **workflow.records[0]["context"]["var"]}
    assert metadata["docs"] == [str(tmp_path / "inputs/sample.json")]
    assert metadata["gain"] == 2 and metadata["temp_dir"] == str(tmp_path / "work")
    assert not (tmp_path / "work").exists()
    workflow.run()
    assert (tmp_path / "products/sample.bin").read_bytes().startswith(b"Any file")


def test_steps_require_explicit_run_and_ids_are_removed(recipe, tmp_path):
    recipe["file_source"].pop("core:run")
    run_workflow(recipe, config_dir=str(tmp_path))
    assert not (tmp_path / "products").exists()
    recipe["file_source"]["id"] = "another_name"
    with pytest.raises(ValueError, match="prefix"):
        validate_config(recipe)
    recipe["file_source"].pop("id")
    recipe["import_files"].pop("core:run")
    assert Workflow(validate_config(recipe)).records == []


def test_shared_settings_require_prefixes(recipe):
    recipe["shared"]["temp_dir"] = "./work"
    with pytest.raises(ValueError, match="prefix"):
        validate_config(recipe)


def test_repeated_plugin_uses_names_and_position_without_ids(recipe, tmp_path):
    recipe.update({'file_source_1': {**(recipe.pop("file_source")), "plugin": 'file_source'}, 'file_source_2': {**(copy_step("final", "copied", "expr:const.output_dir & '/final.bin'")), "plugin": 'file_source'}})
    counts = run_plugin("file_source", recipe, config_dir=str(tmp_path))
    assert list(counts) == ["file_source_1", "file_source_2"]
    assert (tmp_path / "products/final.bin").is_file()


def test_new_registration_automatically_gets_cli_and_only_runs_selected_plugin(
    recipe, tmp_path, monkeypatch, capsys
):
    from vhrharmonize.plugins.base import FunctionPlugin

    calls = []

    class Consumer(FunctionPlugin):
        input_dependency_paths = {"source"}
        input_existence_check_paths = {"source"}
        output_path_resolution_paths = {"destination"}
        output_parent_creation_paths = {"destination"}
        output_dependency_paths = {"destination"}
        output_reuse_paths = {"destination"}
        output_validation_paths = {"destination"}

        def run(self, *, params, shared):
            calls.append(({"source": params["source"]}, params["settings"]))
            shutil.copy2(params["source"], params["destination"])
            return {"written": params["destination"]}

    class Entry:
        name = "consumer"

        def load(self):
            return Consumer

    monkeypatch.setattr(registry, "_entries", lambda: [Entry()])
    with pytest.raises(SystemExit) as help_result:
        main(["--help"])
    assert help_result.value.code == 0
    assert "consumer" in capsys.readouterr().out
    recipe["consumer"] = {"plugin": 'consumer', 
        "core:run": True,
        "param:source": "var:copied",
        "param:settings": "const:$",
        "var:result": "expr:const.output_dir & '/final.bin'",
        "param:destination": "var:result",
    }
    recipe["alignment"] = {"plugin": 'alignment', 
        "core:run": True,
        "param:moving_image_path": "var:result",
        "param:fixed_image_path": "var:mul",
        "param:output_image_path": "expr:const.output_dir & '/later.bin'",
    }
    config = tmp_path / "recipe.yml"
    config.write_text(yaml.safe_dump(recipe, sort_keys=False))
    with pytest.raises(ValueError, match="Unselected step file_source") as direct:
        run_plugin("consumer", str(config))
    with pytest.raises(ValueError) as cli:
        main(["consumer", "--config", str(config)])
    assert str(cli.value) == str(direct.value)
    assert calls == []
    cached = tmp_path / "products" / "sample.bin"
    cached.parent.mkdir()
    cached.write_bytes(b"cached upstream output")
    assert main(["consumer", "--config", str(config), "--dry-run"]) == 0
    assert calls == []
    capsys.readouterr()
    assert main(["consumer", "--config", str(config)]) == 0
    assert len(calls) == 1
    assert calls[0][0]["source"] == str(cached)
    assert calls[0][1]["temp_dir"] == str(tmp_path / "work")
    assert (tmp_path / "products" / "final.bin").read_bytes() == cached.read_bytes()
    assert not (tmp_path / "products" / "later.bin").exists()
    recipe["consumer"]["core:run"] = False
    calls.clear()
    run_plugin("consumer", recipe, config_dir=str(tmp_path))
    assert calls == []


def test_hpc_staging_rewrites_metadata_roots(recipe, tmp_path):
    from vhrharmonize.workflow.staging import stage_workflow

    recipe["import_files"].pop("const:temp_dir")
    original = deepcopy(recipe)
    remote = tmp_path / "remote"
    staged, uploads, downloads = stage_workflow(
        recipe,
        config_dir=str(tmp_path),
        remote_output_dir=str(remote / "out"),
        remote_temp_dir=str(remote / "temp"),
        remote_reference_dir=str(remote / "refs"),
    )
    assert recipe == original
    assert "temp_dir" not in staged["shared"]
    assert staged["import_files"]["core:run"] is False
    assert staged["restore_scenes"]["core:run"] is True
    assert "id" not in staged["import_files"]
    metadata = Workflow(staged).initial_context["const"]
    assert metadata["temp_dir"] == str(remote / "temp")
    assert metadata["output_dir"] == str(remote / "out")
    for local, destination in uploads.items():
        Path(destination).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(local, destination)
    run_workflow(staged)
    assert (remote / "out" / "sample.bin").is_file()
    assert downloads[str(tmp_path / "products" / "sample.bin")] == str(
        remote / "out" / "sample.bin"
    )


def test_all_hpc_commands_have_config_generated_from_python_signature(capsys):
    for command in ["prepare", "upload", "start", "status", "stop", "close", "download"]:
        with pytest.raises(SystemExit) as help_result:
            main(["hpc-" + command, "--help"])
        assert help_result.value.code == 0
        help_text = capsys.readouterr().out
        assert "--config " in help_text
        assert "--config-path" not in help_text


def test_missing_plugin_in_recipe_is_a_python_error(recipe, tmp_path):
    with pytest.raises(ValueError, match="not present"):
        run_plugin("alignment", recipe, config_dir=str(tmp_path))


def test_staged_metadata_roots_expand_on_the_execution_host(recipe, tmp_path, monkeypatch):
    import os
    from vhrharmonize.workflow.staging import stage_workflow

    expanduser = os.path.expanduser
    monkeypatch.setattr(
        os.path,
        "expanduser",
        lambda value: (
            str(tmp_path / "remote" / value[len("~/remote/") :])
            if str(value).startswith("~/remote/")
            else expanduser(value)
        ),
    )
    staged, uploads, _ = stage_workflow(
        recipe,
        config_dir=str(tmp_path),
        remote_output_dir="~/remote/out",
        remote_temp_dir="~/remote/temp",
        remote_reference_dir="~/remote/refs",
    )
    for local, destination in uploads.items():
        destination = Path(os.path.expanduser(destination))
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(local, destination)
    workflow = load_workflow(staged)
    metadata = {**workflow.records[0]["context"]["const"], **workflow.records[0]["context"]["var"]}
    assert metadata["temp_dir"] == str(tmp_path / "remote" / "temp")
    assert metadata["output_dir"] == str(tmp_path / "remote" / "out")
    workflow.run()
    assert (tmp_path / "remote" / "out" / "sample.bin").is_file()


def test_cached_metadata_preserves_current_system_temp_root(recipe, tmp_path, monkeypatch):
    import json

    system = tmp_path / "system"
    system.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(system))
    recipe["import_files"].pop("const:temp_dir")
    recipe["import_files"]["var:custom"] = {"sensor_setting": 3}
    recipe["file_source"]["var:saved"] = "returned:$"
    first = load_workflow(recipe, config_dir=str(tmp_path))
    first.run()
    old_temp = first.records[0]["context"]["const"]["temp_dir"]
    checkpoint = json.loads((tmp_path / "products" / "sample.bin.context.json").read_text())
    assert checkpoint["values"]["var.saved"].endswith("sample.bin")
    assert "temp_dir" not in checkpoint["values"]
    assert "output_dir" not in checkpoint["values"]
    second = load_workflow(recipe, config_dir=str(tmp_path))
    new_temp = second.records[0]["context"]["const"]["temp_dir"]
    assert new_temp != old_temp
    assert second.counts()["file_source"]["loaded"] == 1
    second.run()
    metadata = {**second.records[0]["context"]["const"], **second.records[0]["context"]["var"]}
    assert metadata["temp_dir"] == new_temp
    assert metadata["saved"] == str(tmp_path / "products" / "sample.bin")
