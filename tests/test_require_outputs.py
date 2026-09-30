"""Recipe-owned targets, explicit path parsing and reference-based HPC placement."""

from pathlib import Path
import shutil

import pytest
import yaml

from vhrharmonize.workflow.config import validate_config
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.staging import stage_workflow
from vhrharmonize.workflow.values import PathResolver, references, resolve
from workflow_helpers import context_controls
from workflow_helpers import copy_step, import_settings, install_function


def recipe(tmp_path):
    source = tmp_path / "source.txt"
    source.write_text("source")
    return {
        "shared": {"plugin": "shared", "core:run": True, "core:log_to_console": False,
                   "core:show_progress": False, "core:save_statistics_path": None,
                   "core:load_statistics_path": None, "const:root": str(tmp_path)},
        "files": import_settings(source, tmp_path),
        "prepare": copy_step("prepared", "mul", str(tmp_path / "temp/first.txt")),
        "deliver": {**copy_step("final", "prepared", str(tmp_path / "output/final.txt")),
                    "core:require_outputs": "param:output_path"},
    }


def test_only_yaml_targets_request_work_and_cached_outputs_skip_dependencies(tmp_path):
    config = recipe(tmp_path)
    final = tmp_path / "output/final.txt"
    final.parent.mkdir()
    final.write_text("cached")
    assert Workflow(config).counts()["prepare"]["unused"] == 1
    config["prepare"]["core:require_outputs"] = "param:output_path"
    workflow = Workflow(config)
    workflow.run()
    assert (tmp_path / "temp/first.txt").read_text() == "source"
    assert final.read_text() == "cached"
    config["prepare"]["core:require_outputs"] = False
    config["deliver"]["core:require_outputs"] = False
    workflow = Workflow(config)
    assert all(not node.needed for node in workflow.plan().nodes)


def test_require_outputs_temporary_output_survives_cleanup(tmp_path):
    config = recipe(tmp_path)
    config["prepare"]["core:require_outputs"] = "param:output_path"
    Workflow(config).run()
    assert (tmp_path / "temp/first.txt").exists()
    assert (tmp_path / "output/final.txt").read_text() == "source"


@pytest.mark.parametrize("selection", ["param:missing", "param:input_path"])
def test_require_outputs_selector_rejects_unknown_and_input_parameters(tmp_path, selection):
    config = recipe(tmp_path)
    config["deliver"]["core:require_outputs"] = selection
    with pytest.raises(ValueError, match="declared output parameters"):
        Workflow(config)


@pytest.mark.parametrize("value", [None, 42, {"file": "x"}, [], [["x"]]])
def test_require_outputs_selector_rejects_non_path_values(tmp_path, value):
    config = recipe(tmp_path)
    config["deliver"]["param:output_path"] = value
    with pytest.raises(ValueError):
        Workflow(config)


@pytest.mark.parametrize("value", ["output_path", "var:result", [], [False], 1, None])
def test_invalid_require_outputs_configuration(value):
    with pytest.raises(ValueError, match="core:require_outputs"):
        validate_config({"step": {"core:require_outputs": value}})


def test_path_wrapper_tracks_references_and_resolves_lists(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    resolver = PathResolver(tmp_path)
    context = {"var": {"images": ["./a.tif", "../b.tif", "~/c.tif"]}}
    result = resolve("path:var:images", context, path_resolver=resolver)
    assert result == [str(tmp_path / "a.tif"), str(tmp_path.parent / "b.tif"), str(tmp_path / "home/c.tif")]
    assert references("path:var:images") == {"var.images"}
    assert resolve("var:images", context) == context["var"]["images"]


def test_path_sys_stays_stable_between_planning_and_execution(tmp_path):
    config = recipe(tmp_path)
    config["prepare"] = {
        "const:scratch": "path:sys", **config["prepare"],
        "core:require_outputs": "param:output_path",
        "var:prepared": "path:expr:const.scratch & '/result.txt'",
    }
    config["deliver"]["core:require_outputs"] = False
    workflow = Workflow(config, config_dir=tmp_path)
    output = Path(workflow.nodes[0].params["output_path"])
    try:
        workflow.run()
        assert output.read_text() == "source"
        assert workflow.nodes[0].params["output_path"] == str(output)
        assert len(workflow.path_resolver.system_directories) == 1
    finally:
        for directory in workflow.path_resolver.system_directories.values():
            assert Path(directory).name.startswith("vhr-") and Path(directory).parent != Path("/")
            shutil.rmtree(directory)


def test_path_relative_uses_yaml_directory_for_inputs_and_outputs(tmp_path, monkeypatch):
    config = recipe(tmp_path)
    config["deliver"]["param:input_path"] = "path:./source.txt"
    config["deliver"]["param:output_path"] = "path:../delivered.txt"
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    workflow = Workflow(config, config_dir=tmp_path)
    workflow.run()
    assert (tmp_path.parent / "delivered.txt").read_text() == "source"


def test_require_outputs_lists_request_multiple_outputs(tmp_path, monkeypatch):
    def write(outputs, report):
        for name in [*outputs, report]:
            Path(name).write_text("done")

    install_function(monkeypatch, "multiple", write, output_paths=("outputs", "report"))
    config = {"write": {"plugin": "multiple", "core:run": True,
                        "core:require_outputs": ["param:outputs", "param:report"],
                        "param:outputs": [str(tmp_path / "a"), str(tmp_path / "b")],
                        "param:report": str(tmp_path / "report")}}
    workflow = Workflow(config, config_dir=tmp_path)
    workflow.run()
    assert len(workflow.nodes[0].demanded_paths) == 3


def test_hpc_references_route_required_files_and_download_targets(tmp_path):
    config = recipe(tmp_path)
    config["prepare"]["core:require_outputs"] = "param:output_path"
    staged, uploads, downloads = stage_workflow(
        config, config_dir=tmp_path, remote_work_dir=str(tmp_path / "remote"),
        path_mappings={"const:root": str(tmp_path / "remote/raw"),
                       "const:temp_dir": str(tmp_path / "remote/intermediates"),
                       "const:output_dir": str(tmp_path / "remote/results")},
    )
    assert list(uploads) == [str(tmp_path / "source.txt")]
    assert uploads[str(tmp_path / "source.txt")].startswith(str(tmp_path / "remote/raw"))
    assert downloads[str(tmp_path / "temp/first.txt")].startswith(str(tmp_path / "remote/intermediates"))
    assert downloads[str(tmp_path / "output/final.txt")].startswith(str(tmp_path / "remote/results"))
    for local, remote in uploads.items():
        Path(remote).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(local, remote)
    Workflow(staged, config_dir=tmp_path).run()
    assert all(Path(remote).read_text() == "source" for remote in downloads.values())


def test_hpc_rejects_missing_reference(tmp_path):
    with pytest.raises(ValueError, match="Undefined HPC path_mappings"):
        stage_workflow(recipe(tmp_path), config_dir=tmp_path, remote_work_dir="/remote",
                       path_mappings={"const:typo": "/remote/files"})


def test_hpc_keeps_companions_beside_a_mapped_input(tmp_path):
    config = recipe(tmp_path)
    companion = tmp_path / "source.RPB"
    companion.write_text("RPC")
    config["prepare"]["core:requires"] = "path:./source.RPB"
    _, uploads, _ = stage_workflow(
        config, config_dir=tmp_path, remote_work_dir="/remote",
        path_mappings={"const:root": "/remote/images"},
    )
    assert Path(uploads[str(companion)]).parent == Path(uploads[str(tmp_path / "source.txt")]).parent


def test_hpc_maps_per_scene_lists_at_aggregate_steps(tmp_path, monkeypatch):
    config = recipe(tmp_path)
    config["files"]["var:images"] = ["returned:file_path"]
    config["deliver"]["core:require_outputs"] = False
    seen = []
    install_function(monkeypatch, "batch", lambda images: seen.append(images), scope="aggregate")
    config["batch"] = {"plugin": "batch", "core:run": True, "core:require_outputs": True,
                       "param:images": "collect:images"}
    from workflow_helpers import transfer
    remote = str(tmp_path / "remote")
    staged, uploads, _ = stage_workflow(config, config_dir=tmp_path, remote_work_dir=remote,
                                 path_mappings={"const:root": remote})
    transfer(uploads)
    # Variables, including nested per-scene collections, retain the remote values.
    workflow = Workflow(staged, config_dir=tmp_path)
    assert workflow.records[0]["context"]["var"]["images"][0].startswith(remote + "/")


def test_hpc_downloads_explicit_context_with_selected_output(tmp_path, monkeypatch):
    def produce(first, second):
        Path(first).write_text("first")
        Path(second).write_text("second")
        return 9

    install_function(monkeypatch, "paired", produce, output_paths=("first", "second"))
    first, second = tmp_path / "a.txt", tmp_path / "b.txt"
    config = {
        "shared": {"plugin": "shared", "core:run": True, "core:log_to_console": False,
                   "core:show_progress": False, "core:save_statistics_path": None,
                   "core:load_statistics_path": None, "const:root": str(tmp_path)},
        "produce": {"plugin": "paired", "core:run": True, "core:require_outputs": "param:second",
                    "param:first": str(first), "param:second": str(second),
                    "const:gain": "returned:$"},
    }
    config["produce"].update(context_controls(tmp_path / "gain.json", "const.gain"))
    staged, _, downloads = stage_workflow(config, config_dir=tmp_path,
                                         remote_work_dir=str(tmp_path / "remote"),
                                         path_mappings={"const:root": str(tmp_path / "remote")})
    assert str(first) not in downloads
    assert str(tmp_path / "gain.json") in downloads
    Workflow(staged, config_dir=tmp_path).run()
    for local, remote in downloads.items():
        shutil.copy2(remote, local)
    assert Workflow(config, config_dir=tmp_path).counts()["produce"]["loaded"] == 1


def test_hpc_yaml_accepts_reference_mappings_without_legacy_roots(tmp_path):
    from vhrharmonize.slurm import prepare_slurm_plan

    workflow_file = tmp_path / "workflow.yml"
    workflow_file.write_text(yaml.safe_dump(recipe(tmp_path), sort_keys=False))
    sbatch = tmp_path / "job.sbatch"
    sbatch.write_text('#!/bin/sh\nvhr workflow --config "$1"\n')
    hpc = tmp_path / "hpc.yml"
    hpc.write_text(yaml.safe_dump({
        "workflow_config": str(workflow_file), "slurm_start_file": str(sbatch),
        "ssh_host": "example.invalid", "ssh_user": "test", "run_id": "test",
        "remote_work_dir": "/remote/{run_id}", "remote_log_dir": "/remote/logs",
        "path_mappings": {"const:root": "/remote/{run_id}/sources"},
    }))
    plan = prepare_slurm_plan(str(hpc))
    assert plan["path_mappings"] == {"const:root": "/remote/test/sources"}
    assert plan["uploaded_input_paths"][str(tmp_path / "source.txt")].startswith("/remote/test/sources/")
    assert "remote_output_dir" not in plan and "remote_temp_dir" not in plan
