from collections.abc import Iterable
from datetime import date, datetime, timezone
import importlib
import json
from pathlib import Path
import pytest
import yaml
from vhrharmonize.cli.functions import build_parser, function_cli
from vhrharmonize.workflow.api import load_workflow, run_workflow


def test_signature_drives_cli_types_defaults_and_new_options(tmp_path, capsys):
    calls = []

    def calculate(
        name: str,
        images: Iterable[str],
        *,
        gain: float = 2.5,
        enabled: bool = True,
        metadata: dict | None = None,
        acquired: date | None = None,
    ):
        """Calculate a result.

        Args:
            gain: Multiplier applied to each image.
        """
        calls.append((name, images, gain, enabled, metadata, acquired))

    metadata = tmp_path / "metadata.json"
    metadata.write_text('{"scale": 10000}')
    parser = build_parser(calculate)
    assert "Multiplier applied to each image." in parser.format_help()
    assert "default: 2.5" in parser.format_help()
    function_cli(
        calculate,
        [
            "--name",
            "001",
            "--images",
            '["a.tif", "b.tif"]',
            "--metadata",
            "@" + str(metadata),
            "--acquired",
            "2026-01-02",
            "--no-enabled",
        ],
    )
    assert calls == [("001", ["a.tif", "b.tif"], 2.5, False, {"scale": 10000}, date(2026, 1, 2))]
    function_cli(
        calculate, ["--name", "x", "--images", "[]", "--gain", "3.5", "--metadata", "null"]
    )
    assert calls[-1] == ("x", [], 3.5, True, None, None)


@pytest.mark.parametrize(
    "module,function,args,flags",
    [
        (
            "alignment",
            "align_image_pair",
            {
                "moving_image_path": "in.tif",
                "fixed_image_path": "ref.tif",
                "output_image_path": "out.tif",
                "min_valid_fraction": 0,
            },
            [
                "--moving-image-path",
                "in.tif",
                "--fixed-image-path",
                "ref.tif",
                "--output-image-path",
                "out.tif",
                "--min-valid-fraction",
                "0",
            ],
        ),
        (
            "cloud_mask",
            "cloudmask_raster",
            {
                "input_image_path": "in.tif",
                "output_raster_path": "out.tif",
                "output_mask_path": "mask.tif",
                "red_band_index": 1,
                "green_band_index": 2,
                "nir_band_index": 3,
                "buffer_pixels": -1,
            },
            [
                "--input-image-path",
                "in.tif",
                "--output-raster-path",
                "out.tif",
                "--output-mask-path",
                "mask.tif",
                "--red-band-index",
                "1",
                "--green-band-index",
                "2",
                "--nir-band-index",
                "3",
                "--buffer-pixels",
                "-1",
            ],
        ),
        (
            "atmospheric_correction",
            "run_py6s",
            {"input_raster": "in.tif", "output_raster": "out.tif"},
            ["--input-raster", "in.tif", "--output-raster", "out.tif"],
        ),
    ],
)
def test_plugin_domain_errors_are_identical_in_python_and_cli(module, function, args, flags):
    function = getattr(importlib.import_module("vhrharmonize.plugins." + module), function)
    with pytest.raises(ValueError) as api:
        function(**args)
    with pytest.raises(type(api.value)) as cli:
        function_cli(function, flags)
    assert str(cli.value) == str(api.value)


def test_required_arguments_are_validated_by_the_python_function():

    def compute(required: str):
        return required

    with pytest.raises(TypeError) as direct:
        compute()
    with pytest.raises(TypeError) as cli:
        function_cli(compute, [])
    assert str(cli.value) == str(direct.value)


def test_decorated_api_retains_datetime_types_in_generated_cli():
    from vhrharmonize.plugins.fetch_atmosphere import fetch_modis_water_vapor_for_bbox

    values = vars(
        build_parser(fetch_modis_water_vapor_for_bbox).parse_args(
            ["--scene-datetime-utc", "2026-01-02T12:00:00Z"]
        )
    )
    assert values["scene_datetime_utc"] == datetime(2026, 1, 2, 12, tzinfo=timezone.utc)


def test_workflow_api_and_generated_cli_share_planning_and_execution(
    tmp_path, make_test_raster, capsys
):
    from vhrharmonize.cli.main import main

    source = make_test_raster(tmp_path / "source.tif")
    config = {
        "shared": {"plugin": "shared", "core:run": True, "core:log_to_console": False},
        "import_files": {
            "plugin": "import_files",
            "core:run": True,
            "param:search_glob": str(source),
            "var:raw": "returned:file_path",
            "var:input": "returned:file_path",
            "var:input_name": "expr:$split(var.file_path, '/')[-1]",
            "var:input_stem": r"expr:$replace($split(var.file_path, '/')[-1], /\.[^.]*$/, '')",
            "var:input_dir": r"expr:$replace(var.file_path, /[^\/]+$/, '')",
            "const:output_dir": "out",
        },
        "file_source": {
            "plugin": "file_source",
            "core:run": True, "core:require_outputs": True,
            "param:input_path": "var:raw",
            "var:copied": "expr:const.output_dir & '/copied.tif'",
            "param:output_path": "var:copied",
        },
    }
    recipe = tmp_path / "recipe.yml"
    recipe.write_text(yaml.safe_dump(config, sort_keys=False))
    expected = run_workflow(str(recipe), dry_run=True)
    assert main(["workflow", "--config", str(recipe), "--dry-run"]) == 0
    assert json.loads(capsys.readouterr().out) == expected
    assert not (tmp_path / "out").exists()
    assert run_workflow(config, config_dir=str(tmp_path)) == expected
    assert (tmp_path / "out/copied.tif").read_bytes() == source.read_bytes()
    assert load_workflow(recipe).counts()["file_source"]["loaded"] == 1


def test_workflow_configuration_errors_are_identical(tmp_path):
    from vhrharmonize.cli.main import main

    recipe = tmp_path / "bad.yml"
    recipe.write_text("workflow: {alignment: true}\n")
    with pytest.raises(ValueError) as api:
        run_workflow(str(recipe))
    with pytest.raises(ValueError) as cli:
        main(["workflow", "--config", str(recipe)])
    assert str(cli.value) == str(api.value)


def test_hpc_errors_come_from_the_python_api(tmp_path):
    from vhrharmonize import slurm
    from vhrharmonize.cli.main import main

    config = tmp_path / "hpc.yml"
    config.write_text("workflow_config: missing.yml\n")
    with pytest.raises(ValueError) as api:
        slurm.prepare_slurm_plan(str(config), overrides={"run_id": "test"})
    with pytest.raises(ValueError) as cli:
        main(["hpc-prepare", "--config", str(config), "--overrides", '{"run_id":"test"}'])
    assert str(cli.value) == str(api.value)


def test_hpc_download_failure_is_not_swallowed(tmp_path, monkeypatch):
    from vhrharmonize import slurm
    from vhrharmonize.cli.main import main

    config = tmp_path / "staged.yml"
    config.write_text(
        yaml.safe_dump(
            {"download_output_paths": {str(tmp_path / "out.tif"): "/remote/out.tif"}},
            sort_keys=False,
        )
    )
    monkeypatch.setattr(slurm, "_remote_is_directory", lambda *args: False)

    def fail(*args):
        raise OSError("transfer failed")

    monkeypatch.setattr(slurm, "_scp_download", fail)
    with pytest.raises(RuntimeError) as api:
        slurm.download_slurm_outputs(str(config))
    with pytest.raises(RuntimeError) as cli:
        main(["hpc-download", "--config", str(config)])
    assert str(cli.value) == str(api.value)
    assert str(cli.value.__cause__) == "transfer failed"


def test_spectralmatch_cli_discovers_native_function_options(tmp_path, monkeypatch):
    from vhrharmonize.cli.functions import function_cli
    from vhrharmonize.plugins.global_regression import global_regression

    module = importlib.import_module("vhrharmonize.plugins.global_regression")
    calls = []
    monkeypatch.setattr(module, "_upstream", lambda: lambda **kwargs: calls.append(kwargs))
    function_cli(
        global_regression,
        [
            "--input-images",
            '["input.tif"]',
            "--output-images",
            '["output.tif"]',
            "--window-scales",
            "null",
            "--pif-max-samples",
            "null",
            "--cache",
            "0.5",
            "--vector-mask",
            '["include", "mask.gpkg", "image"]',
        ],
    )
    assert calls[0]["input_images"] == ["input.tif"]
    assert calls[0]["output_images"] == ["output.tif"]
    assert calls[0]["window_scales"] is None
    assert calls[0]["pif_max_samples"] is None
    assert calls[0]["cache"] == 0.5
    assert calls[0]["vector_mask"] == ("include", "mask.gpkg", "image")


def test_qgis_csv_file_output_is_a_python_api_option(tmp_path):
    from vhrharmonize.plugins.orthorectification import qgis_gcps_to_csv

    source = tmp_path / "points.txt"
    source.write_text("1,2,3,4,1,0,0,0\n")
    destination = tmp_path / "result" / "points.csv"
    result = qgis_gcps_to_csv(str(source), output_csv_path=str(destination))
    assert destination.read_text() == result
