import importlib
import json
from pathlib import Path
from unittest.mock import Mock

import pytest

from vhrharmonize.cli import worldview


@pytest.mark.parametrize("named", [False, True])
@pytest.mark.parametrize("output_kind", ["missing", "file", "corrupt_file", "folder"])
def test_spectralmatch_always_runs_without_output_checks(tmp_path, make_test_raster, monkeypatch, named, output_kind):
    args = worldview._build_parser().parse_args([])
    args.calculate_overviews_spectralmatch = True
    args.overview_scales = [2, 4]
    args.match_steps = ["global_regression", "align"]
    args.match_shared_resume_from_steps = "validate"
    output = tmp_path / ("timelapse" if output_kind == "folder" else "mosaic.tif")
    if output_kind == "folder":
        output.mkdir()
        (output / "partial.tif").write_text("incomplete")
    elif output_kind == "file":
        make_test_raster(output)
    elif output_kind == "corrupt_file":
        output.write_bytes(b"broken TIFF")
    before = output.read_bytes() if output.is_file() else None
    pipeline = Mock()
    monkeypatch.setattr(worldview, "spectralmatch", pipeline)
    for name in ("_prepare_step_outputs", "_gdal_raster_is_valid", "calculate_raster_overviews"):
        monkeypatch.setattr(worldview, name, Mock(side_effect=AssertionError("Output checking belongs to SpectralMatch")))

    if named:
        result = worldview._run_named_spectralmatch_group(
            output.name, "auto:*.tif", available_paths=["scene.tif"], args=args,
            temp_root=str(tmp_path / "temp"), output_root=str(tmp_path),
        )
    else:
        result = worldview._run_default_spectralmatch(
            ["scene.tif"], args=args, output_path=str(output), temp_root=str(tmp_path / "temp"),
        )

    assert result == str(output)
    pipeline.assert_called_once()
    assert pipeline.call_args.kwargs["global_regression_build_overviews"] is True
    assert pipeline.call_args.kwargs["shared_resume_from_steps"] == "validate"
    if before is not None:
        assert output.read_bytes() == before


def test_spectralmatch_wrapper_trusts_successful_pipeline_return(tmp_path, monkeypatch):
    module = importlib.import_module("vhrharmonize.preprocess.spectralmatch")
    pipeline = Mock(return_value={"output": ["one.tif", "two.tif"]})
    monkeypatch.setattr(module, "spectralmatch_pipeline", pipeline)
    exists = Mock(return_value=False)
    monkeypatch.setattr(Path, "exists", exists)

    assert module.spectralmatch(["input.tif"], str(tmp_path / "timelapse")) == str(tmp_path / "timelapse")

    pipeline.assert_called_once()
    exists.assert_not_called()


def test_standalone_spectralmatch_merges_rasters(tmp_path, make_test_raster):
    module = importlib.import_module("vhrharmonize.preprocess.spectralmatch")
    source = make_test_raster(tmp_path / "input.tif", width=32, height=32)
    output = tmp_path / "mosaic.tif"

    assert module.spectralmatch(
        [str(source)], str(output), steps=["merge"], merge_rasters_build_overviews=False,
        shared_image_threads=1, shared_io_threads=1, shared_tile_threads=1,
    ) == str(output)

    with worldview.gdal.Open(str(output)) as dataset:
        assert dataset.RasterXSize == dataset.RasterYSize == 32
        assert (dataset.ReadAsArray() == 1).all()


@pytest.mark.parametrize("steps, expected", [
    (None, "merge_rasters_build_overviews"),
    (["joint_coregistration"], "joint_coregistration_build_overviews"),
    (["global_regression", "local_block_adjustment", "align"], "local_block_adjustment_build_overviews"),
    (["local_block_adjustment", "global_regression", "mask"], "global_regression_build_overviews"),
    (["global_regression", "merge"], "merge_rasters_build_overviews"),
])
def test_overviews_choose_last_eligible_step(steps, expected):
    args = worldview._build_parser().parse_args(["--calculate-overviews-spectralmatch"])
    args.match_steps = steps
    args.overview_scales = [2, 4]

    kwargs = worldview._build_spectralmatch_kwargs(args)

    assert {key: value for key, value in kwargs.items() if key.endswith("_build_overviews")} == {expected: True}


@pytest.mark.parametrize("source", ["match", "json"])
@pytest.mark.parametrize("option", list(worldview.SPECTRALMATCH_OVERVIEW_STEPS.values()))
def test_automatic_overviews_reject_any_explicit_enabled_flag(source, option):
    args = worldview._build_parser().parse_args(["--calculate-overviews-spectralmatch"])
    args.overview_scales = [2]
    args.match_steps = ["global_regression"]
    if source == "match":
        setattr(args, "match_" + option, True)
    else:
        args.spectralmatch_kwargs_json = json.dumps({option: True})

    with pytest.raises(ValueError, match="cannot be combined.*match_" + option):
        worldview._build_spectralmatch_kwargs(args)


def test_automatic_overviews_use_json_steps_and_explicit_scales():
    args = worldview._build_parser().parse_args(["--calculate-overviews-spectralmatch"])
    args.spectralmatch_kwargs_json = '{"steps": ["global_regression", "align"]}'
    args.match_shared_window_scales = [2, 4]
    args.match_global_regression_build_overviews = False

    kwargs = worldview._build_spectralmatch_kwargs(args)

    assert kwargs["global_regression_build_overviews"] is True
    assert kwargs["shared_window_scales"] == [2, 4]


@pytest.mark.parametrize("steps", [[], ["align"], ["mask", "weighted_seamline"]])
def test_automatic_overviews_reject_pipeline_without_eligible_steps(steps):
    args = worldview._build_parser().parse_args(["--calculate-overviews-spectralmatch"])
    args.match_steps = steps
    args.overview_scales = [2]
    with pytest.raises(ValueError, match="overview-capable"):
        worldview._build_spectralmatch_kwargs(args)


def test_automatic_overviews_require_scales():
    args = worldview._build_parser().parse_args(["--calculate-overviews-spectralmatch"])
    with pytest.raises(ValueError, match="requires overview_scales"):
        worldview._build_spectralmatch_kwargs(args)


def test_overview_conflicts_fail_before_processing(monkeypatch, capsys):
    run = Mock()
    monkeypatch.setattr(worldview, "_run_workflow", run)
    with pytest.raises(SystemExit):
        worldview.main([
            "--input-file-glob", "*.tif", "--run-spectralmatch", "--calculate-overviews-spectralmatch",
            "--overview-scales", "2", "4", "--match-local-block-adjustment-build-overviews", "true",
        ])
    assert "cannot be combined" in capsys.readouterr().err
    run.assert_not_called()


@pytest.mark.parametrize("legacy", ["run_radiometric_normalization", "save_radiometric_normalization", "calculate_overviews_radiometric_normalization"])
def test_old_step_config_keys_are_rejected(legacy):
    with pytest.raises(ValueError, match="Use spectralmatch"):
        worldview._normalize_config_defaults({legacy: True})


def test_new_cli_and_null_save_config_reach_workflow(tmp_path, monkeypatch):
    config = tmp_path / "config.yml"
    config.write_text("workflow:\n  run_spectralmatch: true\n  save_spectralmatch:\n")
    run = Mock(return_value=0)
    monkeypatch.setattr(worldview, "_run_workflow", run)
    assert worldview.main([
        "--config-yaml", str(config), "--input-file-glob", "*.tif",
        "--match-shared-output-image-path", str(tmp_path / "timelapse"),
    ]) == 0
    assert run.call_args.args[0].run_spectralmatch


def test_processing_counts_do_not_inspect_spectralmatch_outputs(make_worldview_bundle, tmp_path, monkeypatch, capsys):
    bundle = make_worldview_bundle()
    scene = worldview.load_worldview_scenes_from_tif_files([str(bundle["mul_tif"]), str(bundle["pan_tif"])])[0]
    args = worldview._build_parser().parse_args(["--run-spectralmatch", "--log-to-console"])
    args.temp_dir = str(tmp_path / "temp")
    args.save_spectralmatch = str(tmp_path / "timelapse")
    checked = []
    check = worldview._existing_outputs_are_reusable

    def inspect(paths, **kwargs):
        checked.extend(paths)
        return check(paths, **kwargs)

    monkeypatch.setattr(worldview, "_existing_outputs_are_reusable", inspect)
    counts = worldview._count_processing_steps([scene], args)
    worldview._log_processing_steps(args, counts)

    assert args.save_spectralmatch not in checked
    assert counts["spectralmatch"] == {"loaded": 0, "processing": 1}
    line = next(line for line in capsys.readouterr().out.splitlines() if "spectralmatch:" in line)
    assert "handled by SpectralMatch" in line and "loaded:" not in line
