"""Upload reporting uses real rsync records without contacting an HPC host."""

from io import StringIO
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
from prompt_toolkit.utils import get_cwidth
import yaml

from vhrharmonize import slurm
from vhrharmonize.cli.functions import function_cli
from vhrharmonize.workflow import upload_progress as upload
from test_selective_hpc_staging import recipe


def inventory(tmp_path):
    folder = tmp_path / "source"
    folder.mkdir()
    (folder / "image one.tif").write_bytes(b"a" * 4096)
    (folder / "metadata.txt").write_bytes(b"b" * 1024)
    (folder / "empty.txt").touch()
    (folder / "café.txt").write_bytes(b"c" * 512)
    items = [("uploaded_input_paths", str(folder), "/remote/products")]
    groups = [{"step": "alignment", "variable": "var:output_dir", "files": [str(folder)]}]
    return folder, items, groups


def test_partial_bytes_rates_snapshots_and_current_files(tmp_path, monkeypatch):
    _, items, groups = inventory(tmp_path)
    clock = [10.0]
    monkeypatch.setattr(upload, "monotonic", lambda: clock[0])
    snapshots = []
    progress = upload.UploadProgress(items, groups=groups, callbacks=[snapshots.append])
    progress.begin("/remote", items)
    progress.consume("VHR:>f+++++++:4096:products/image one.tif")
    clock[0] = 12
    progress.consume("         2,048  50%  1.00kB/s 0:00:02")
    partial = progress.snapshot()
    assert partial["total"]["files"] == 4
    assert partial["total"]["done"] == 0
    assert partial["total"]["transferred_bytes"] == 2048
    assert partial["total"]["rate"] == 1024
    assert partial["total"]["eta_seconds"] == 3.5
    assert partial["rows"][0]["step"] == "alignment"
    assert partial["rows"][0]["variable"] == "var:output_dir"
    assert "image one" not in json.dumps(partial)
    progress.consume("         4,096 100%  1.00kB/s 0:00:00 (xfr#1)")
    # Repeated progress records cannot double count.
    progress.consume("         4,096 100%  1.00kB/s 0:00:00 (xfr#1)")
    clock[0] = 14
    progress.finish()
    final = progress.snapshot()
    assert final["total"]["done"] == 4
    assert final["total"]["current"] == 3
    assert final["total"]["current_bytes"] == 1536
    assert final["total"]["transferred_bytes"] == 4096
    assert final["total"]["eta_seconds"] == 0
    clock[0] = 40
    assert progress.snapshot() == final  # Final timings are stable.
    snapshots[-1]["rows"][0]["done"] = -1
    assert progress.snapshot() == final  # Callbacks receive detached snapshots.


@pytest.mark.parametrize("error,status", [(RuntimeError("rsync failed"), "failed"), (KeyboardInterrupt(), "cancelled")])
def test_failed_batches_do_not_claim_unfinished_files(tmp_path, error, status):
    _, items, groups = inventory(tmp_path)
    progress = upload.UploadProgress(items, groups=groups)
    progress.begin("/remote", items)
    progress.consume("VHR:>f+++++++:4096:products/image one.tif")
    progress.consume("2048 50% 1kB/s 0:00:02")
    progress.finish(error=error)
    result = progress.snapshot()
    assert result["status"] == result["rows"][0]["status"] == status
    assert result["total"]["done"] == result["total"]["current"] == 0
    assert result["total"]["transferred_bytes"] == 2048


def test_directory_and_individual_entry_count_once_and_collisions_fail(tmp_path):
    folder, items, groups = inventory(tmp_path)
    overlap = ("uploaded_input_paths", str(folder / "image one.tif"), "/remote/products/image one.tif")
    progress = upload.UploadProgress([*items, overlap], groups=groups)
    assert progress.snapshot()["total"]["files"] == 4
    assert progress.snapshot()["total"]["total_bytes"] == 5632
    with pytest.raises(ValueError, match="Multiple local files"):
        upload.UploadProgress([*items, (overlap[0], str(folder / "metadata.txt"), overlap[2])])


@pytest.mark.skipif(shutil.which("rsync") is None, reason="local rsync is not installed")
def test_real_local_rsync_first_transfer_and_repeat(tmp_path, capsys):
    folder, _, groups = inventory(tmp_path)
    destination = tmp_path / "destination"
    destination.mkdir()
    items = [("uploaded_input_paths", str(folder), str(destination))]
    for repeat in (False, True):
        progress = upload.UploadProgress(items, groups=groups)
        progress.begin(str(destination), items)
        slurm._run_local_command(["rsync", "-aL", "--itemize-changes", "--progress",
                                  "--out-format=" + upload.RSYNC_FORMAT, str(folder) + "/", str(destination) + "/"],
                                 output_callback=progress.consume)
        progress.finish()
        result = progress.snapshot()["total"]
        assert result["done"] == result["files"] == 4
        assert result["current"] == (4 if repeat else 0)
        assert result["transferred_bytes"] == (0 if repeat else 5632)
        assert result["current_bytes"] == (5632 if repeat else 0)
    assert (destination / "café.txt").read_bytes() == (folder / "café.txt").read_bytes()
    assert capsys.readouterr().out == ""


def test_streaming_callback_cancellation_reaps_child(monkeypatch):
    real_popen = subprocess.Popen
    children = []

    def popen(*args, **kwargs):
        child = real_popen(*args, **kwargs)
        children.append(child)
        return child

    def cancel(line):
        assert line == "ready"
        raise KeyboardInterrupt

    monkeypatch.setattr(slurm.subprocess, "Popen", popen)
    with pytest.raises(KeyboardInterrupt):
        slurm._run_local_command([sys.executable, "-c", "import time; print('ready', flush=True); time.sleep(60)"],
                                 output_callback=cancel)
    assert children[0].poll() is not None


@pytest.mark.parametrize("width", [60, 80, 120, 160])
def test_upload_rendering_preserves_columns_colors_and_rightmost_bar(tmp_path, width):
    _, items, groups = inventory(tmp_path)
    progress = upload.UploadProgress(items, groups=groups)
    progress.begin("/remote", items)
    progress.consume("VHR:>f+++++++:4096:products/image one.tif")
    progress.consume("2048 50% 1kB/s 0:00:02")
    stream = StringIO()
    display = upload.TerminalUploadDisplay(stream=stream, width=width)
    with display:
        display.update(progress.snapshot())
        display.update(progress.snapshot())
    text = stream.getvalue()
    assert text.count("VHRHarmonize HPC Upload") == 1
    assert "image one" not in text and "metadata.txt" not in text
    assert "total" in text and "From" in text and "Files" in text
    assert all(get_cwidth(line) == width for line in text.splitlines())
    row = next(line for line in text.splitlines() if "total" in line)
    assert row.endswith("▬ │")
    if width >= 120:
        assert "var:output_dir" in text and "Speed" in text and "ETA" in text
    colors = {style for style, value in display.render() if "▬" in value}
    assert f"fg:{upload.COLORS['done']}" in colors
    assert f"fg:{upload.COLORS['run']}" in colors
    assert display.read_input is False


@pytest.mark.parametrize("show_progress", [False, True])
def test_upload_cli_grouped_report_without_file_dump(tmp_path, monkeypatch, capsys, show_progress):
    folder, items, groups = inventory(tmp_path)
    config = tmp_path / "staged.hpc.yml"
    config.write_text(yaml.safe_dump(dict(ssh_user="user", ssh_host="example.invalid",
        remote_work_dir="/remote", uploaded_input_paths={items[0][1]: items[0][2]},
        uploaded_reference_paths={}, upload_groups=groups, debug_logs=True, show_progress=show_progress)))
    commands = []

    def local(command, **kwargs):
        commands.append(command)
        assert "--progress" in command
        assert "--info=progress2" not in command
        assert "--out-format=" + upload.RSYNC_FORMAT in command
        callback = kwargs["output_callback"]
        callback("VHR:>f+++++++:4096:products/image one.tif")
        callback("4096 100% 4kB/s 0:00:00")
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(slurm, "_run_local_command", local)
    assert function_cli(slurm.upload_slurm_files, ["--config", str(config)]) == 0
    output = capsys.readouterr().out
    assert "image one.tif" not in output and "uploaded_input_paths" not in output
    assert ("VHRHarmonize HPC Upload" in output) is show_progress
    assert "4/4" in output
    assert len(commands) == 1  # Still one batched rsync, not one process per file.
    saved = yaml.safe_load(config.read_text())
    assert saved["status"] == "uploaded"
    assert saved["upload_results"][str(folder)]["status"] == "synced"


def test_staging_labels_imports_products_and_mapped_roots(tmp_path):
    config, source = recipe(tmp_path)
    groups = []
    from vhrharmonize.workflow.engine import Workflow
    from vhrharmonize.workflow.staging import stage_workflow

    Workflow(config, config_dir=tmp_path).run()
    _, uploads, _ = stage_workflow(config, config_dir=tmp_path, remote_work_dir="/remote",
        path_mappings={"const:root": "/remote"}, upload_groups=groups)
    labels = {filename: (group["step"], group["variable"]) for group in groups for filename in group["files"]}
    assert set(labels) == set(uploads)
    assert labels[str(source)] == ("import", "const:root")
    assert labels[str(tmp_path / "products/scene.txt")] == ("copy", "const:root")


def test_prepare_records_control_group_and_progress_default(tmp_path):
    config, _ = recipe(tmp_path)
    workflow = tmp_path / "workflow.yml"
    workflow.write_text(yaml.safe_dump(config, sort_keys=False))
    job = tmp_path / "job.sbatch"
    job.write_text("#!/bin/sh\n")
    hpc = tmp_path / "hpc.yml"
    hpc.write_text(yaml.safe_dump(dict(workflow_config=str(workflow), slurm_start_file=str(job),
        ssh_user="user", ssh_host="example.invalid", remote_work_dir="/remote", remote_log_dir="/remote/logs",
        path_mappings={"const:root": "/remote"})))
    result = slurm.prepare_slurm_plan(str(hpc))
    assert result["show_progress"] is True
    controls = next(group for group in result["upload_groups"] if group["step"] == "controls")
    assert set(controls["files"]) == set(result["uploaded_reference_paths"])
    assert "show_progress" in slurm.SLURM_PREPARE_CONFIG_KEYS


def test_empty_transfer_finishes_and_legacy_maps_still_get_totals(tmp_path):
    progress = upload.UploadProgress([])
    progress.finish()
    assert progress.snapshot()["status"] == "done"
    assert progress.snapshot()["total"]["files"] == 0
    _, items, _ = inventory(tmp_path)
    legacy = upload.UploadProgress(items, path_mappings={"var:output_dir": "/remote/products"})
    assert legacy.snapshot()["rows"][0]["step"] == "inputs"
    assert legacy.snapshot()["rows"][0]["variable"] == "var:output_dir"


def test_live_upload_leaves_terminal_input_available(tmp_path, monkeypatch):
    from vhrharmonize.workflow import progress_terminal as terminal
    from test_progress_terminal import ScreenOutput

    class TtyStream(StringIO):
        def isatty(self):
            return True

    def forbidden_input(**kwargs):
        raise AssertionError("Upload dashboard must leave real terminal input to SSH")

    output = ScreenOutput()
    monkeypatch.setattr(terminal, "create_input", forbidden_input)
    monkeypatch.setattr(terminal, "create_output", lambda **kwargs: output)
    monkeypatch.setattr(terminal.sys, "stdin", TtyStream())
    monkeypatch.setenv("TERM", "xterm-256color")
    _, items, groups = inventory(tmp_path)
    progress = upload.UploadProgress(items, groups=groups)
    display = upload.TerminalUploadDisplay(stream=TtyStream())
    display.update(progress.snapshot())
    with display:
        assert display.ready.is_set()
        assert display.thread.is_alive()
    assert display.error is None
    assert not display.thread.is_alive()
    assert not output.entered and not output.mouse_enabled
    assert output.cursor_visible


@pytest.mark.skipif(shutil.which("rsync") is None, reason="local rsync is not installed")
def test_upload_directory_mapped_to_workspace_root_counts_actual_transfers(tmp_path, monkeypatch, capsys):
    folder, _, groups = inventory(tmp_path)
    remote = tmp_path / "remote"
    snapshots = []

    def local_rsync(data, *, stage_root, remote_root, progress_callback):
        return slurm._run_local_command(["rsync", "-aL", "--progress", "--out-format=" + upload.RSYNC_FORMAT,
            stage_root + "/", remote_root + "/"], output_callback=progress_callback)

    monkeypatch.setattr(slurm, "_rsync_upload_tree", local_rsync)
    config = tmp_path / "staged.hpc.yml"
    config.write_text(yaml.safe_dump(dict(ssh_user="user", ssh_host="example.invalid",
        remote_work_dir=str(remote), uploaded_input_paths={str(folder): str(remote)},
        uploaded_reference_paths={}, upload_groups=groups, show_progress=False)))
    slurm.upload_slurm_files(str(config), progress_callback=snapshots.append)
    assert snapshots[-1]["status"] == "done"
    assert snapshots[-1]["total"]["transferred_bytes"] == 5632
    assert snapshots[-1]["total"]["current"] == 0
    assert (remote / "image one.tif").read_bytes() == (folder / "image one.tif").read_bytes()
    assert not (remote / "remote").exists()
    assert "image one.tif" not in capsys.readouterr().out


@pytest.mark.skipif(shutil.which("rsync") is None, reason="local rsync is not installed")
def test_scene_folders_outside_workspace_share_one_rsync_batch(tmp_path, monkeypatch, capsys):
    source = tmp_path / "source"
    source.mkdir()
    input_root = tmp_path / "inputs"
    uploads, images, companions = {}, [], []
    for scene in ("one", "two"):
        local = source / (scene + ".tif")
        local.write_bytes(b"x" * 1024)
        uploads[str(local)] = str(input_root / scene / local.name)
        images.append(str(local))
    # A smaller companion group should join the enclosing input batch too.
    sidecar = source / "one.RPB"
    sidecar.write_bytes(b"rpc")
    uploads[str(sidecar)] = str(input_root / "one/one.RPB")
    companions.append(str(sidecar))
    groups = [{"step": "import", "variable": "var:images", "files": images},
              {"step": "import", "variable": "var:companions", "files": companions}]
    calls, snapshots = [], []

    def local_rsync(data, *, stage_root, remote_root, progress_callback):
        calls.append(remote_root)
        Path(remote_root).mkdir(parents=True, exist_ok=True)
        return slurm._run_local_command(["rsync", "-aL", "--progress", "--out-format=" + upload.RSYNC_FORMAT,
            stage_root + "/", remote_root + "/"], output_callback=progress_callback)

    monkeypatch.setattr(slurm, "_rsync_upload_tree", local_rsync)
    result = slurm._upload_required_files(dict(remote_work_dir=str(tmp_path / "run"), uploaded_input_paths=uploads,
        uploaded_reference_paths={}, upload_groups=groups, show_progress=False), progress_callback=snapshots.append)
    assert calls == [str(input_root)]
    assert len(result) == 3
    assert snapshots[-1]["total"]["done"] == 3
    assert snapshots[-1]["total"]["transferred_bytes"] == 2051
    assert {row["variable"] for row in snapshots[-1]["rows"]} == {"var:images", "var:companions"}
    assert (input_root / "two/two.tif").read_bytes() == b"x" * 1024
    assert "one.tif" not in capsys.readouterr().out
