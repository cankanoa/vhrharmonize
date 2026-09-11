from __future__ import annotations

from pathlib import Path

import pytest

from vhrharmonize.cli.cli_helpers import _load_yaml_config
from vhrharmonize.io.workflow_utils import (
    build_output_path_from_input,
    plan_step_outputs,
    remove_output_files,
    resolve_output_dir,
    resolve_relative_to_input,
    resolve_temp_dir,
)


def test_load_yaml_config_normalizes_keys(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yml"
    config_path.write_text("input-dir: abc\nlog-to-console: true\n", encoding="utf-8")
    loaded = _load_yaml_config(str(config_path))
    assert loaded == {"input_dir": "abc", "log_to_console": True}


def test_workflow_utils_paths(tmp_path: Path) -> None:
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    assert resolve_relative_to_input("a/b", str(input_dir)) == str(input_dir / "a/b")
    assert resolve_relative_to_input("~/a/b", str(input_dir)) == str(Path.home() / "a/b")
    assert resolve_output_dir(None, temp_dir=str(tmp_path / "temp"), step_name="step").endswith("step")
    assert Path(resolve_temp_dir(str(tmp_path / "custom_temp"), input_folder=str(input_dir))).exists()
    assert build_output_path_from_input(str(input_dir / "image.tif"), str(tmp_path), suffix="_x") == str(tmp_path / "image_x.tif")


def test_plan_step_outputs_skips_existing(tmp_path: Path) -> None:
    input_path = tmp_path / "image.tif"
    input_path.write_text("x", encoding="utf-8")
    existing_output = tmp_path / "out" / "image_done.tif"
    existing_output.parent.mkdir(parents=True, exist_ok=True)
    existing_output.write_text("y", encoding="utf-8")
    plan = plan_step_outputs([str(input_path)], output_dir=str(existing_output.parent), suffix="_done", skip_existing=True)
    assert plan.output_paths == [str(existing_output)]
    assert plan.pending_output_paths == []


@pytest.mark.parametrize("extension", [".tif", ".dat", ".img", ".vrt", ".json"])
def test_remove_output_files_preserves_sidecars(tmp_path: Path, extension: str) -> None:
    output = tmp_path / f"output{extension}"
    sidecars = [Path(f"{output}{suffix}") for suffix in (".aux.xml", ".msk", ".ovr", ".hdr", ".params.txt")]
    sidecars.append(output.with_suffix(".hdr"))
    for path in [output, *sidecars]:
        path.write_bytes(b"corrupt")
    unrelated = tmp_path / "other.tif"
    unrelated.write_bytes(b"keep")

    remove_output_files([str(output)], input_paths=[str(unrelated)])

    assert not output.exists()
    assert all(path.read_bytes() == b"corrupt" for path in sidecars)
    assert unrelated.read_bytes() == b"keep"


@pytest.mark.parametrize("alias", ["same_path", "hardlink", "symlink"])
def test_remove_output_files_protects_inputs_before_any_deletion(tmp_path: Path, alias: str) -> None:
    source = tmp_path / "source.tif"
    source.write_bytes(b"source")
    output = source
    if alias != "same_path":
        output = tmp_path / "alias.tif"
        if alias == "hardlink":
            output.hardlink_to(source)
        else:
            output.symlink_to(source)
    other_output = tmp_path / "other.tif"
    other_output.write_bytes(b"keep until all checks pass")

    with pytest.raises(ValueError, match="also an input"):
        remove_output_files([str(other_output), str(output)], input_paths=[str(source)])

    assert other_output.exists()
    assert source.read_bytes() == b"source"


def test_remove_output_files_does_not_delete_directories(tmp_path: Path) -> None:
    directory = tmp_path / "tiles.tif"
    directory.mkdir()
    tile = directory / "tile.tif"
    tile.write_bytes(b"keep")
    with pytest.raises(IsADirectoryError):
        remove_output_files([str(directory)])
    assert tile.read_bytes() == b"keep"
