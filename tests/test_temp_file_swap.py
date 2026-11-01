"""Selected file arguments use isolated local copies without changing the plan."""

from pathlib import Path
from uuid import UUID

import pytest

from vhrharmonize.workflow.config import validate_config
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.temp_files import TempFileSwap, TempFileSwapError
from workflow_helpers import copy_step, import_settings, install_function, stage, transfer


def swap(tmp_path, names=("inputs", "outputs")):
    return TempFileSwap(names, str(tmp_path / "swap"), frozenset({"inputs"}), frozenset({"outputs"}))


def test_each_file_gets_its_own_uuid_and_original_basename(tmp_path):
    inputs, outputs = [], []
    for name in ("a", "b"):
        source = tmp_path / name / "same.txt"
        source.parent.mkdir()
        source.write_text(name)
        inputs.append(str(source))
        outputs.append(str(tmp_path / "result" / name / "same.txt"))
    params = {"inputs": inputs, "outputs": outputs, "untouched": inputs[0]}
    with swap(tmp_path).apply(params) as (local, originals):
        paths = [Path(p) for p in local["inputs"] + local["outputs"]]
        assert len({p.parent for p in paths}) == 4
        assert all(p.name == "same.txt" and UUID(p.parent.name) for p in paths)
        assert all(p.parent.parent == tmp_path / "swap" for p in paths)
        assert local["untouched"] == inputs[0]
        assert originals == dict(zip(map(str, paths), inputs + outputs))
        for source, output in zip(local["inputs"], local["outputs"]):
            Path(output).write_text(Path(source).read_text() + " processed")
    assert [Path(p).read_text() for p in outputs] == ["a processed", "b processed"]
    assert [Path(p).read_text() for p in inputs] == ["a", "b"]
    assert not list((tmp_path / "swap").iterdir())
    assert params["inputs"] == inputs and params["outputs"] == outputs


def test_plugin_failure_preserves_originals_and_cleans_copies(tmp_path):
    source, destination = tmp_path / "source.txt", tmp_path / "output.txt"
    source.write_text("original input")
    destination.write_text("existing output")
    with pytest.raises(RuntimeError, match="plugin failed"):
        with swap(tmp_path).apply({"inputs": str(source), "outputs": str(destination)}) as (local, _):
            assert Path(local["outputs"]).read_text() == "existing output"
            Path(local["outputs"]).write_text("incomplete")
            raise RuntimeError("plugin failed")
    assert source.read_text() == "original input"
    assert destination.read_text() == "existing output"
    assert not list((tmp_path / "swap").iterdir())


def test_failed_copy_back_retains_completed_output_but_removes_input_copy(tmp_path, monkeypatch):
    from vhrharmonize.workflow import temp_files

    source, destination = tmp_path / "source.txt", tmp_path / "output.txt"
    source.write_text("original")
    copy = temp_files.shutil.copyfile
    def copy_file(before, after, **kwargs):
        if str(after) == str(destination):
            raise OSError("NAS unavailable")
        return copy(before, after, **kwargs)
    monkeypatch.setattr(temp_files.shutil, "copyfile", copy_file)
    with pytest.raises(OSError, match="retained at"):
        with swap(tmp_path).apply({"inputs": str(source), "outputs": str(destination)}) as (local, _):
            Path(local["outputs"]).write_text("complete output")
    assert not Path(local["inputs"]).exists()
    assert Path(local["outputs"]).read_text() == "complete output"
    assert not destination.exists()


@pytest.mark.parametrize("existing_output", [False, True])
def test_copy_back_succeeds_when_destination_rejects_metadata_changes(
    tmp_path, monkeypatch, existing_output,
):
    from vhrharmonize.workflow import temp_files

    source, destination = tmp_path / "source.txt", tmp_path / "output.txt"
    source.write_text("original input")
    if existing_output:
        destination.write_text("previous output")
    copystat = temp_files.shutil.copystat

    def restricted_copystat(before, after, **kwargs):
        if Path(after) == destination:
            raise PermissionError(1, "Operation not permitted")
        return copystat(before, after, **kwargs)

    monkeypatch.setattr(temp_files.shutil, "copystat", restricted_copystat)
    with swap(tmp_path).apply({"inputs": str(source), "outputs": str(destination)}) as (local, _):
        Path(local["outputs"]).write_text("complete output")
    assert destination.read_text() == "complete output"
    assert source.read_text() == "original input"
    assert not list((tmp_path / "swap").iterdir())


def test_missing_output_is_not_published_and_optional_files_are_unchanged(tmp_path):
    with pytest.raises(TempFileSwapError, match="Temporary output was not created"):
        with swap(tmp_path).apply({"inputs": None, "outputs": str(tmp_path / "missing.txt")}) as (local, _):
            assert local["inputs"] is None
    assert not list((tmp_path / "swap").iterdir())


def test_declared_companions_stay_beside_the_correct_input_only(tmp_path):
    inputs, requirements = [], []
    for name in ("a", "b"):
        folder = tmp_path / name
        folder.mkdir()
        source = folder / "image.TIF"
        source.write_text(name)
        companion = folder / "image.RPB"
        companion.write_text(name + " RPC")
        (folder / "image.IMD").write_text("not declared")
        inputs.append(str(source))
        requirements.append(str(companion))
    spec = swap(tmp_path)
    spec.requirements = tuple(requirements)
    with spec.apply({"inputs": inputs, "outputs": None}) as (local, _):
        for source in map(Path, local["inputs"]):
            assert source.with_suffix(".RPB").read_text() == source.read_text() + " RPC"
            assert not source.with_suffix(".IMD").exists()
            assert source.parent.parent == tmp_path / "swap"
    assert not list((tmp_path / "swap").iterdir())


def test_companion_copy_failure_errors_and_cleans_up_before_plugin_call(tmp_path, monkeypatch):
    from vhrharmonize.workflow import temp_files

    source = tmp_path / "image.TIF"
    source.write_text("input")
    companion = source.with_suffix(".RPB")
    companion.write_text("RPC")
    spec = swap(tmp_path)
    spec.requirements = (str(companion),)
    copy = temp_files.shutil.copy2
    def fail_companion(before, after):
        if Path(before) == companion:
            raise OSError("companion transfer failed")
        return copy(before, after)
    monkeypatch.setattr(temp_files.shutil, "copy2", fail_companion)
    with pytest.raises(TempFileSwapError, match="companion transfer failed"):
        with spec.apply({"inputs": str(source), "outputs": None}):
            pytest.fail("Plugin must not run without the required companion")
    assert not list((tmp_path / "swap").iterdir())
    assert companion.read_text() == "RPC"


@pytest.mark.parametrize("value", ["input_path", ["param:input_path"], ["input_path", "input_path"], [1]])
def test_selection_requires_bare_unique_parameter_names(value):
    with pytest.raises(ValueError, match="swap_for_temp_file"):
        validate_config({"copy": {"plugin": "file_source", "core:swap_for_temp_file": value}})


def recipe(tmp_path, *, workers=1, direction="vertical"):
    for name in ("a", "b"):
        filename = tmp_path / "source" / (name + ".txt")
        filename.parent.mkdir(exist_ok=True)
        filename.write_text(name)
    return {
        "shared": {"plugin": "shared", "core:run": True, "core:swap_temp_dir": "const:temp_dir",
                   "core:concurrent_processing": workers, "core:processing_direction": direction,
                   "core:show_progress": True, "core:log_to_console": False,
                   "core:save_statistics_path": None, "core:load_statistics_path": None},
        "source": import_settings(tmp_path / "source/*.txt", tmp_path),
        "copy": {**copy_step("copied", "mul", suffix="_copy", folder="output_dir", require_outputs=True),
                 "core:swap_for_temp_file": ["input_path", "output_path"],
                 "var:reported": "returned:$"},
    }


@pytest.mark.parametrize("workers,direction", [(1, "vertical"), (2, "vertical"), (2, "horizontal")])
def test_workflow_keeps_original_paths_for_returns_reuse_and_worker_progress(tmp_path, monkeypatch, workers, direction):
    config = recipe(tmp_path, workers=workers, direction=direction)
    workflow = Workflow(config, config_dir=tmp_path)
    originals = [dict(n.params) for n in workflow.nodes]
    records = workflow.run()
    for record in records:
        values = record["context"]["var"]
        destination = str(tmp_path / "output" / (values["basename"] + "_copy.txt"))
        assert values["reported"] == values["copied"] == destination
        assert Path(destination).read_text() == values["basename"]
    assert [n.params for n in workflow.nodes] == originals
    assert not list((tmp_path / "temp").iterdir())
    assert workflow.get_progress()["status"] == "completed"
    from vhrharmonize.workflow import temp_files
    monkeypatch.setattr(temp_files, "uuid4", lambda: pytest.fail("Reused files must not be copied"))
    resumed = Workflow(config, config_dir=tmp_path)
    resumed.run()
    assert resumed.counts()["copy"]["processing"] == 0


def test_hpc_plans_original_files_and_resolves_swap_folder_on_remote(tmp_path):
    config = recipe(tmp_path)
    staged, uploads, _ = stage(config, tmp_path)
    assert not (tmp_path / "temp").exists()
    assert staged["shared"]["core:swap_temp_dir"] == "const:temp_dir"
    transfer(uploads)
    remote = Workflow(staged, config_dir=tmp_path)
    records = remote.run()
    assert all("/remote/output/" in r["context"]["var"]["reported"] for r in records)
    assert not list((tmp_path / "remote/temp").iterdir())


def test_unknown_selection_and_missing_folder_fail_without_copying(tmp_path):
    config = recipe(tmp_path)
    config["copy"]["core:swap_for_temp_file"] = ["gain"]
    with pytest.raises(ValueError, match="declared file parameters"):
        Workflow(config, config_dir=tmp_path)
    config["copy"]["core:swap_for_temp_file"] = ["input_path"]
    config["shared"].pop("core:swap_temp_dir")
    with pytest.raises(ValueError, match="swap_temp_dir must be set"):
        Workflow(config, config_dir=tmp_path).run()
    assert not (tmp_path / "temp").exists()


def test_workflow_supplies_declared_companions_to_temporary_input(tmp_path, monkeypatch):
    observed = []
    def read_companion(input_path, output_path):
        content = Path(input_path).with_suffix(".RPB").read_text()
        observed.append(input_path)
        Path(output_path).write_text(content)
        return output_path
    install_function(monkeypatch, "read_companion", read_companion,
                     input_paths={"input_path"}, output_paths={"output_path"})
    config = recipe(tmp_path)
    config["copy"].update({"plugin": "read_companion",
        "core:requires": "expr:$replace(var.mul, /\\.txt$/, '.RPB')"})
    for name in ("a", "b"):
        (tmp_path / "source" / (name + ".RPB")).write_text(name + " RPC")
    # Keep only the intended primary files in discovery.
    workflow = Workflow(config, config_dir=tmp_path)
    workflow.run()
    assert len(observed) == 2
    assert all(Path(p).parent.parent == tmp_path / "temp" for p in observed)
    assert (tmp_path / "output/a_copy.txt").read_text() == "a RPC"
    assert not list((tmp_path / "temp").iterdir())


@pytest.mark.parametrize("workers", [1, 2])
@pytest.mark.parametrize("policy", ["remove_and_restart", "remove_instence"])
def test_transfer_errors_stop_workflow_even_with_image_removal_policy(tmp_path, workers, policy):
    config = recipe(tmp_path, workers=workers)
    bad_root = tmp_path / "not_a_folder"
    bad_root.write_text("file")
    config["shared"].update({"core:swap_temp_dir": str(bad_root), "core:worker_error_case": policy})
    workflow = Workflow(config, config_dir=tmp_path)
    with pytest.raises(TempFileSwapError, match="Could not prepare temporary file"):
        workflow.run()
    assert not workflow._removed_scene_ids
    assert workflow.errors[0]["error_type"] == "TempFileSwapError"
