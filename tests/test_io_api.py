from __future__ import annotations

from pathlib import Path

from vhrharmonize.cli.cli_helpers import load_yaml_config
from vhrharmonize.io.workflow_utils import (
    build_output_path_from_input,
    plan_step_outputs,
    resolve_output_dir,
    resolve_relative_to_input,
    resolve_temp_dir,
)


def test_load_yaml_config_normalizes_keys(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yml"
    config_path.write_text("input-dir: abc\nlog-to-console: true\n", encoding="utf-8")
    loaded = load_yaml_config(str(config_path))
    assert loaded == {"input_dir": "abc", "log_to_console": True}


def test_workflow_utils_paths(tmp_path: Path) -> None:
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    assert resolve_relative_to_input("a/b", str(input_dir)) == str(input_dir / "a/b")
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
