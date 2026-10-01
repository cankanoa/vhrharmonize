from copy import deepcopy
import json
from pathlib import Path
import shutil

import pytest
import yaml

from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.staging import stage_workflow
from vhrharmonize.slurm import prepare_slurm_plan
from test_explicit_context import shared, importer


def recipe(tmp_path):
    source = tmp_path / "raw" / "scene.txt"
    source.parent.mkdir()
    source.write_text("data")
    config = {"settings": shared(**{"const:root": "path:" + str(tmp_path)}),
              "import": importer(source),
              "copy": {"plugin": "file_source", "core:run": True, "core:require_outputs": "param:output_path",
                       "param:input_path": "var:raw", "var:result": "expr:const.root & '/products/scene.txt'",
                       "param:output_path": "var:result"}}
    return config, source


def transfer(uploads):
    for local, remote in uploads.items():
        Path(remote).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(local, remote)


def test_direct_root_mapping_preserves_expressions_and_relative_layout(tmp_path):
    config, source = recipe(tmp_path)
    remote = tmp_path / "remote"
    staged, uploads, downloads = stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote),
                                               path_mappings={"const:root": str(remote)})
    assert list(staged) == list(config)
    assert staged["copy"] == config["copy"]
    assert staged["settings"]["const:root"] == "path:" + str(remote)
    assert uploads[str(source)] == str(remote / "raw/scene.txt")
    assert downloads[str(tmp_path / "products/scene.txt")] == str(remote / "products/scene.txt")
    transfer(uploads)
    Workflow(staged, config_dir=tmp_path).run()
    assert (remote / "products/scene.txt").read_text() == "data"


def test_missing_mapping_and_collisions_fail_but_identity_mapping_is_allowed(tmp_path):
    config, source = recipe(tmp_path)
    with pytest.raises(ValueError, match="no HPC root mapping"):
        stage_workflow(config, config_dir=tmp_path, remote_work_dir="/remote")
    staged, uploads, _ = stage_workflow(config, config_dir=tmp_path, remote_work_dir="/remote",
                                       path_mappings={"const:root": str(tmp_path)})
    assert uploads[str(source)] == str(source)
    (tmp_path / "other").mkdir()
    other = tmp_path / "other/scene.txt"
    other.write_text("different")
    config["settings"].update({"const:raw_root": str(source.parent), "const:other_root": str(other.parent)})
    config["import"]["param:search_glob"] = [str(source), str(other)]
    config["import"].pop("param:scene_id")
    with pytest.raises(ValueError, match="collision"):
        stage_workflow(config, config_dir=tmp_path, remote_work_dir="/remote",
                       path_mappings={"const:root": "/remote", "const:raw_root": "/remote/files", "const:other_root": "/remote/files"})


def test_explicit_import_context_moves_without_raw_discovery_snapshot(tmp_path):
    config, source = recipe(tmp_path)
    snapshot = tmp_path / "context/import.json"
    config["import"]["core:save_context"] = {str(snapshot): "defined"}
    config["import"]["core:load_context"] = {str(snapshot): "defined"}
    Workflow(config, config_dir=tmp_path)  # explicitly save the import as it runs
    original = snapshot.read_bytes()
    remote = tmp_path / "remote"
    staged, uploads, _ = stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote),
                                       path_mappings={"const:root": str(remote)}, context_staging_dir=tmp_path / "staged-context")
    assert snapshot.read_bytes() == original
    assert "restore_scenes" not in staged
    assert not any("staged_" in key for block in staged.values() for key in block)
    assert str(remote / "context/import.json") in uploads.values()
    transfer(uploads)
    Workflow(staged, config_dir=tmp_path).run()
    assert (remote / "products/scene.txt").read_text() == "data"


@pytest.mark.parametrize("debug_logs", [False, True])
def test_prepare_cutoff_preserves_both_yaml_documents_and_does_not_process_later(tmp_path, capsys, debug_logs):
    config, source = recipe(tmp_path)
    config["settings"].update({"core:log_to_console": True, "core:show_progress": True})
    context = tmp_path / "context/import.json"
    config["import"]["core:save_context"] = {str(context): "defined"}
    config["import"]["core:load_context"] = {str(context): "defined"}
    workflow = tmp_path / "workflow.yml"
    text = "# original recipe\n" + yaml.safe_dump(config, sort_keys=False)
    text = text.replace("  plugin: file_source", "  plugin: file_source # keep this comment")
    workflow.write_text(text)
    job = tmp_path / "job.sbatch"
    job.write_text('#!/bin/sh\nvhr workflow --config "$1"\n')
    hpc = tmp_path / "hpc.yml"
    hpc.write_text("# user HPC comments\n" + yaml.safe_dump({
        "workflow_config": str(workflow), "slurm_start_file": str(job), "ssh_host": "example.invalid", "ssh_user": "user",
        "run_id": "test", "remote_work_dir": str(tmp_path / "remote"), "remote_log_dir": str(tmp_path / "remote/logs"),
        "debug_logs": debug_logs,
        "run_to_step_before_prepare": "import", "path_mappings": {"const:root": str(tmp_path / "remote")},
    }, sort_keys=False))
    plan = prepare_slurm_plan(str(hpc))
    output = capsys.readouterr().out
    if debug_logs:
        assert output.startswith("[hpc:prepare] Start\n")
        local = output.index("[hpc:prepare] Running local steps | through=import")
        staging = output.index("[hpc:prepare] Staging workflow and planning transfers")
        assert local < output.index("[core:workflow] Start") < staging
        assert staging < output.rindex("[core:workflow] Start")
    else:
        assert "[hpc:prepare]" not in output
    assert "VHRHarmonize Workflow Progress" not in output
    assert "[core:cleanup] Start" in output
    assert context.exists()
    assert not (tmp_path / "products/scene.txt").exists()
    assert workflow.read_text() == text
    prepared = Path(plan["staged_workflow_file"]).with_suffix(".prepare.yml").read_text()
    assert "# original recipe" in prepared and "# keep this comment" in prepared
    assert yaml.safe_load(prepared)["copy"]["core:run"] is False
    staged = Path(plan["staged_workflow_file"]).read_text()
    assert "# original recipe" in staged and "# keep this comment" in staged
    assert yaml.safe_load(staged)["copy"]["core:run"] is True
    assert Path(plan["staged_hpc_file"]).read_text().startswith("# user HPC comments\n")


def test_root_mapping_changes_first_assignment_and_keeps_later_mutations(tmp_path):
    config, source = recipe(tmp_path)
    before_copy = config.pop("copy")
    config["nested"] = {"core:run": True, "const:root": "expr:const.root & '/nested'"}
    config["copy"] = before_copy
    remote = tmp_path / "remote"
    staged, uploads, downloads = stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote),
                                               path_mappings={"const:root": str(remote)})
    assert staged["nested"] == config["nested"]
    transfer(uploads)
    Workflow(staged, config_dir=tmp_path).run()
    assert (remote / "nested/products/scene.txt").read_text() == "data"
    assert downloads[str(tmp_path / "nested/products/scene.txt")] == str(remote / "nested/products/scene.txt")


def test_unused_context_save_needs_no_mapping_and_adds_no_transfer(tmp_path):
    config, _ = recipe(tmp_path)
    config["unused"] = {"core:run": True, "const:gain": 8,
                        "core:save_context": {"/unmapped/unused.json": "const.gain"}}
    _, _, downloads = stage_workflow(config, config_dir=tmp_path, remote_work_dir="/remote",
                                     path_mappings={"const:root": "/remote"})
    assert "/unmapped/unused.json" not in downloads


def test_all_import_context_and_cached_intermediate_avoid_source_upload(tmp_path):
    config, source = recipe(tmp_path)
    snapshot = tmp_path / "context/import.json"
    config["import"].update({"core:save_context": {str(snapshot): "all"},
                             "core:load_context": {str(snapshot): "all"}})
    Workflow(config, config_dir=tmp_path).run()
    intermediate = tmp_path / "products/scene.txt"
    config["copy"]["core:require_outputs"] = False
    config["finish"] = {"plugin": "file_source", "core:run": True,
                        "core:require_outputs": "param:output_path",
                        "param:input_path": "var:result",
                        "param:output_path": "expr:const.root & '/finished.txt'"}
    source.unlink()
    remote = tmp_path / "remote"
    staged, uploads, _ = stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote),
                                       path_mappings={"const:root": str(remote)},
                                       context_staging_dir=tmp_path / "staged-context")
    assert str(source) not in uploads
    assert str(intermediate) in uploads
    assert set(uploads.values()) == {str(remote / "products/scene.txt"), str(remote / "context/import.json")}
    transfer(uploads)
    Workflow(staged, config_dir=tmp_path).run()
    assert (remote / "finished.txt").read_text() == "data"
