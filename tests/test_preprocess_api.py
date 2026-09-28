from __future__ import annotations

import sys
import importlib
from pathlib import Path
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest
import rasterio
from shapely.geometry import box
import yaml

import vhrharmonize.plugins.atmospheric_correction as atmos_mod
from vhrharmonize.plugins.fetch_dem import download_opentopography_dem_for_bbox
import vhrharmonize.plugins.fetch_atmosphere as fetch_mod
from vhrharmonize.plugins.alignment import align_image_pair
from vhrharmonize.plugins.atmospheric_correction import (
    atmospheric_correction,
    _build_flaash_kwargs_from_variables,
    _convert_flaash_params_paths_for_windows,
    _init_envi_engine,
    parallel_flaash,
    run_flaash,
    _run_flaash_wrapper,
    run_py6s,
    _validate_flaash_params,
    _wsl_path_to_windows_for_envi,
)
from vhrharmonize.plugins.cloud_mask import (
    apply_binary_cloud_mask_to_image,
    cloudmask_raster,
    create_cloud_mask_with_omnicloudmask,
)
from vhrharmonize.plugins.fetch_atmosphere import (
    fetch_modis_water_vapor_for_bbox,
    fetch_power_atmosphere_for_bbox,
    _init_ee_client,
)


def _metadata():
    return dict(
        solar_zenith=17.0,
        solar_azimuth=75.0,
        view_zenith=22.0,
        view_azimuth=80.0,
        line_of_sight_azimuth=80.0,
        line_of_sight_zenith=158.0,
        day=5,
        month=7,
        band_wavelengths_um=[0.5, 0.8],
        dn_to_radiance_factors=[1.0, 3.0],
        dn_to_radiance_offsets=[-1.0, -2.0],
    )


@pytest.mark.parametrize("delete_temp", [None, False, True])
def test_align_image_pair(monkeypatch, tmp_path: Path, delete_temp) -> None:
    calls = []
    alignment = importlib.import_module("vhrharmonize.plugins.alignment")
    monkeypatch.setattr(
        alignment,
        "coregix_align_image_pair",
        lambda **kwargs: calls.append(kwargs)
        or SimpleNamespace(output_image_path=kwargs["output_image_path"]),
    )
    options = {} if delete_temp is None else {"delete_temp_dir": delete_temp}
    output = align_image_pair("moving.tif", "fixed.tif", str(tmp_path / "out.tif"), **options)
    assert output.output_image_path.endswith("out.tif")
    assert calls[0]["keep_temp_dir"] is (delete_temp is False)
    assert calls[0]["split_factor"] == 0
    assert calls[0]["output_on_moving_grid"] is False


def test_align_image_pair_rejects_removed_keep_temp_dir() -> None:
    with pytest.raises(TypeError, match="keep_temp_dir"):
        align_image_pair("moving.tif", "fixed.tif", "out.tif", keep_temp_dir=True)


@pytest.mark.parametrize(
    "flags, delete_temp",
    [([], True), (["--delete-temp-dir"], True), (["--no-delete-temp-dir"], False)],
)
def test_align_image_pair_cli_cleanup_option(
    monkeypatch, tmp_path: Path, flags, delete_temp
) -> None:
    from vhrharmonize.cli.functions import function_cli

    image = tmp_path / "input.tif"
    image.touch()
    calls = []
    monkeypatch.setattr(
        importlib.import_module("vhrharmonize.plugins.alignment"),
        "coregix_align_image_pair",
        lambda **kwargs: calls.append(kwargs)
        or SimpleNamespace(output_image_path=kwargs["output_image_path"]),
    )

    assert (
        function_cli(
            align_image_pair,
            [
                "--moving-image-path",
                str(image),
                "--fixed-image-path",
                str(image),
                "--output-image-path",
                str(tmp_path / "out.tif"),
                *flags,
            ],
        )
        == 0
    )

    assert calls[0]["keep_temp_dir"] is not delete_temp


@pytest.mark.parametrize("flag", ["--keep-temp-dir", "--no-keep-temp-dir"])
def test_align_image_pair_cli_rejects_removed_flag(flag, capsys) -> None:
    from vhrharmonize.cli.functions import function_cli
    from vhrharmonize.cli.functions import build_parser

    with pytest.raises(SystemExit):
        build_parser(align_image_pair).parse_args(
            [
                "--moving-image-path",
                "moving.tif",
                "--fixed-image-path",
                "fixed.tif",
                "--output-image-path",
                "out.tif",
                flag,
            ]
        )
    assert f"unrecognized arguments: {flag}" in capsys.readouterr().err


def test_py6s_helpers_and_run(monkeypatch, tmp_path: Path) -> None:
    calls = []
    monkeypatch.setattr(atmos_mod.Py6SCorrector, "run", lambda self, **kwargs: calls.append(kwargs))
    result = run_py6s(
        "in.tif", str(tmp_path / "out.tif"), variables=_metadata(), ground_elevation_km=1.0
    )
    assert result.output_raster.endswith("out.tif")
    assert calls[0]["dn_to_radiance_factors"] == [1.0, 3.0]
    assert calls[0]["ground_elevation_km"] == 1.0


def test_flaash_helpers_and_dispatch(monkeypatch, tmp_path: Path) -> None:
    fake_geospatial = ModuleType("vhrharmonize.io.geospatial")
    fake_geospatial.get_image_percentile_value = lambda *args, **kwargs: 100.0
    monkeypatch.setattr(
        atmos_mod, "get_image_percentile_value", fake_geospatial.get_image_percentile_value
    )
    kwargs = _build_flaash_kwargs_from_variables(
        "in.tif",
        "dem.tif",
        box(0, 0, 1, 1),
        _metadata(),
        "out.tif",
        dem_ground_percentile=50.0,
        modtran_atm="Mid-Latitude Summer",
        modtran_aer="Rural",
        use_aerosol="AUTO",
        default_visibility=10.0,
    )
    assert kwargs["GROUND_ELEVATION"] == 0.1
    assert _validate_flaash_params({"MODTRAN_ATM": "x"})["MODTRAN_ATM"] == "x"
    assert _wsl_path_to_windows_for_envi("/mnt/c/test/file") == "C:\\test\\file"
    assert (
        _convert_flaash_params_paths_for_windows(
            {
                "INPUT_RASTER": {"url": "/mnt/c/x", "factory": "URLRaster"},
                "OUTPUT_RASTER_URI": "/mnt/c/y",
            }
        )["OUTPUT_RASTER_URI"]
        == "C:\\y"
    )
    fake_envipyengine = ModuleType("envipyengine")
    fake_envipyengine.Engine = lambda *_: SimpleNamespace(tasks=lambda: None)
    fake_envipyengine_config = ModuleType("envipyengine.config")
    fake_envipyengine_config.set = lambda *args, **kwargs: None
    fake_envipyengine.config = fake_envipyengine_config
    monkeypatch.setitem(sys.modules, "envipyengine", fake_envipyengine)
    monkeypatch.setitem(sys.modules, "envipyengine.config", fake_envipyengine_config)
    assert _init_envi_engine("engine").tasks() is None
    monkeypatch.setattr(atmos_mod, "_execute_flaash_task", lambda *args, **kwargs: None)
    assert (
        _run_flaash_wrapper(({"OUTPUT_RASTER_URI": "out.tif"}, "params.txt", object())) == "out.tif"
    )
    monkeypatch.setattr(atmos_mod, "_execute_flaash_task", lambda *args, **kwargs: None)
    result = run_flaash(
        "in.tif",
        str(tmp_path / "out.tif"),
        params={
            "INPUT_RASTER": {"url": "in.tif", "factory": "URLRaster"},
            "OUTPUT_RASTER_URI": str(tmp_path / "out.tif"),
        },
        envi_engine=object(),
    )
    assert result.output_raster.endswith("out.tif")
    monkeypatch.setattr(atmos_mod, "_run_flaash_wrapper", lambda task: task[0]["OUTPUT_RASTER_URI"])

    class _Future:
        def __init__(self, value):
            self._value = value

        def result(self):
            return self._value

    class _Executor:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def submit(self, fn, task):
            return _Future(fn(task))

    monkeypatch.setattr(atmos_mod, "ProcessPoolExecutor", _Executor)
    monkeypatch.setattr(atmos_mod, "as_completed", lambda futures: futures)
    monkeypatch.setattr(atmos_mod, "tqdm", lambda it, **kwargs: it)
    assert parallel_flaash([({"OUTPUT_RASTER_URI": "a.tif"}, "p.txt")], object()) == ["a.tif"]
    monkeypatch.setattr(
        atmos_mod.Py6SCorrector,
        "run",
        lambda self, input_raster, output_raster, **kw: output_raster,
    )
    assert (
        atmospheric_correction(
            "in.tif",
            "out.tif",
            method="py6s",
            solar_zenith=1,
            solar_azimuth=1,
            view_zenith=1,
            view_azimuth=1,
            day=1,
            month=1,
        )
        == "out.tif"
    )


def test_cloudmask_functions(monkeypatch, tmp_path: Path, make_test_raster) -> None:
    input_path = make_test_raster(
        tmp_path / "input.tif", count=3, data=np.arange(48, dtype=np.uint16).reshape(3, 4, 4)
    )
    mask_path = tmp_path / "mask.tif"
    output_path = tmp_path / "masked.tif"
    monkeypatch.setitem(
        sys.modules,
        "omnicloudmask",
        SimpleNamespace(
            predict_from_array=lambda arr, **kwargs: (
                np.array([[0, 1], [1, 0]])
                if arr.shape[-1] == 2
                else np.array([[0, 1, 0, 0], [0, 0, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0]])
            )
        ),
    )
    create_cloud_mask_with_omnicloudmask(
        str(input_path), str(mask_path), 1, 2, 3, inference_resolution_m=None
    )
    apply_binary_cloud_mask_to_image(
        str(input_path), str(mask_path), str(output_path), custom_nodata_value=999
    )
    result = cloudmask_raster(
        str(input_path),
        str(tmp_path / "masked2.tif"),
        str(tmp_path / "mask2.tif"),
        red_band_index=1,
        green_band_index=2,
        nir_band_index=3,
        inference_resolution_m=None,
    )
    assert result.mask_pixel_count > 0


def test_fetch_functions(monkeypatch, tmp_path: Path) -> None:
    class _Response:
        headers = {"content-type": "image/tiff"}
        content = b"abc"

        def raise_for_status(self):
            return None

        def json(self):
            return {}

        text = ""

    monkeypatch.setattr(fetch_mod.requests, "get", lambda *args, **kwargs: _Response())
    monkeypatch.setenv("OPENTOPOGRAPHY_API_KEY", "key")
    assert download_opentopography_dem_for_bbox(
        min_lon=0, min_lat=0, max_lon=1, max_lat=1, output_tif_path=str(tmp_path / "dem.tif")
    ).endswith("dem.tif")
    monkeypatch.setattr(
        fetch_mod,
        "_fetch_power_daily_point",
        lambda *args, **kwargs: {"aot550": 0.2, "water_vapor": 2.5, "ozone_cm_atm": 0.3},
    )
    power = fetch_power_atmosphere_for_bbox(
        day_utc=__import__("datetime").date(2020, 1, 1), min_lon=0, min_lat=0, max_lon=1, max_lat=1
    )
    assert power.aot550 == pytest.approx(0.2)
    fake_ee = SimpleNamespace(
        Authenticate=lambda: None,
        Initialize=lambda **kwargs: None,
    )
    monkeypatch.setitem(sys.modules, "ee", fake_ee)
    assert _init_ee_client(env_file=None) is fake_ee
    monkeypatch.setattr(
        fetch_mod,
        "_fetch_collection_value",
        lambda *args, **kwargs: {
            "collection": "terra",
            "band_found": True,
            "raw_value": 2000.0,
            "abs_time_diff_hours": 1.0,
            "image_time_utc": "2020-01-01T00:00:00+00:00",
        },
    )
    monkeypatch.setattr(
        fetch_mod,
        "_fetch_first_available_band_value",
        lambda *args, **kwargs: {
            "collection": "terra",
            "band_found": True,
            "raw_value": 100.0,
            "abs_time_diff_hours": 1.0,
            "image_time_utc": "2020-01-01T00:00:00+00:00",
        },
    )
    modis = fetch_modis_water_vapor_for_bbox(
        scene_datetime_utc=__import__("datetime").datetime(
            2020, 1, 1, tzinfo=__import__("datetime").timezone.utc
        ),
        min_lon=0,
        min_lat=0,
        max_lon=1,
        max_lat=1,
        ee=object(),
    )
    assert modis.status == "ok"


def test_pansharpen(monkeypatch, tmp_path: Path) -> None:
    class _PanSharpen:
        def __init__(self, pan, mul):
            self.pan, self.mul = pan, mul

        def process(self, output_image_path, write_mask=False, overwrite=True):
            with rasterio.open(self.mul) as src:
                profile = src.profile.copy()
                data = src.read()
            with rasterio.open(output_image_path, "w", **profile) as dst:
                dst.write(data)

    fake_orthority = ModuleType("orthority")
    fake_orthority.PanSharpen = _PanSharpen
    monkeypatch.setitem(sys.modules, "orthority", fake_orthority)
    pansharpen_mod = importlib.import_module("vhrharmonize.plugins.pansharpen")
    mul = tmp_path / "mul.tif"
    pan = tmp_path / "pan.tif"
    data = np.ones((1, 4, 4), dtype=np.int16)
    with rasterio.open(
        mul,
        "w",
        driver="GTiff",
        width=4,
        height=4,
        count=1,
        dtype="int16",
        transform=rasterio.transform.from_origin(0, 4, 1, 1),
        crs="EPSG:4326",
        nodata=-32768,
    ) as dst:
        dst.write(data)
    with rasterio.open(
        pan,
        "w",
        driver="GTiff",
        width=4,
        height=4,
        count=1,
        dtype="int16",
        transform=rasterio.transform.from_origin(0, 4, 1, 1),
        crs="EPSG:4326",
        nodata=-32768,
    ) as dst:
        dst.write(data)
    out = tmp_path / "ps.tif"
    pansharpen_mod.pansharpen_image(str(mul), str(pan), str(out), custom_nodata_value=0)
    assert out.exists()



@pytest.mark.parametrize(
    "mode, expected",
    [
        ("no", {"missing.tif"}),
        ("yes", {"missing.tif", "valid.tif", "broken.tif", "run.log"}),
        ("validate", {"missing.tif", "broken.tif"}),
        (None, {"missing.tif", "broken.tif"}),
    ],
)
def test_slurm_download_conflicts(tmp_path: Path, monkeypatch, mode, expected) -> None:
    slurm_mod = importlib.import_module("vhrharmonize.slurm")
    valid = tmp_path / "valid.tif"
    with rasterio.open(
        valid,
        "w",
        driver="GTiff",
        width=4,
        height=4,
        count=1,
        dtype="uint8",
        transform=rasterio.transform.from_origin(0, 4, 1, 1),
        crs="EPSG:4326",
    ) as dst:
        dst.write(np.ones((1, 4, 4), dtype=np.uint8))
    (tmp_path / "broken.tif").write_text("incomplete raster")
    (tmp_path / "run.log").write_text("existing log")
    config = {
        "download_output_paths": {
            str(tmp_path / name): f"/remote/{name}"
            for name in ("missing.tif", "valid.tif", "broken.tif")
        },
        "download_log_paths": {str(tmp_path / "run.log"): "/remote/run.log"},
    }
    config_path = tmp_path / "staged.hpc.yml"
    # Write literal yes/no to cover YAML boolean normalization.
    config_path.write_text(
        yaml.safe_dump(config)
        + (f"override_download_conflict: {mode}\n" if mode is not None else "")
    )
    monkeypatch.setattr(slurm_mod, "_remote_is_directory", lambda *args: False)
    downloaded = []
    monkeypatch.setattr(
        slurm_mod, "_scp_download", lambda data, remote, local: downloaded.append(Path(local).name)
    )
    slurm_mod.download_slurm_outputs(str(config_path))
    assert set(downloaded) == expected


def test_slurm_download_rejects_invalid_conflict_policy(tmp_path: Path, monkeypatch) -> None:
    slurm_mod = importlib.import_module("vhrharmonize.slurm")
    config_path = tmp_path / "staged.hpc.yml"
    config_path.write_text("override_download_conflict: invalid\n")
    with pytest.raises(ValueError, match="override_download_conflict"):
        slurm_mod.download_slurm_outputs(str(config_path))


def test_start_slurm_yaml_writer_and_remote_quote(tmp_path: Path) -> None:
    from vhrharmonize.slurm import (
        _parse_sbatch_log_templates,
        _remote_quote,
        _resolve_remote_sbatch_log_templates,
        _status_command,
        _status_from_text,
        _resolve_staged_hpc_file,
        _resolve_staged_workflow_file,
        _resolve_staged_slurm_start_file,
        _write_sectioned_yaml_file,
    )

    output_path = tmp_path / "RUN123.staged.hpc.yml"
    long_path = "/" + "/".join(["very_long_path_segment"] * 12) + "/image.tif"
    _write_sectioned_yaml_file(
        str(output_path),
        {"uploaded_input_paths": {long_path: "~/koa_scratch/run/output/image.tif"}},
        header_by_key={"uploaded_input_paths": "# All mappings are local file: remote file."},
    )

    text = output_path.read_text(encoding="utf-8")
    assert (
        yaml.safe_load(text)["uploaded_input_paths"][long_path]
        == "~/koa_scratch/run/output/image.tif"
    )
    assert _remote_quote("~/koa_scratch/run/output") == "~/koa_scratch/run/output"
    assert _resolve_staged_hpc_file(
        {}, config_path=str(tmp_path / "example.hpc.yml"), run_id="RUN123"
    ) == str(tmp_path / "RUN123.staged.hpc.yml")
    assert _resolve_staged_workflow_file(
        {},
        config_path=str(tmp_path / "example.hpc.yml"),
        workflow_config=str(tmp_path / "example.worldview.yml"),
        run_id="RUN123",
    ) == str(tmp_path / "RUN123.staged.worldview.yml")
    assert _resolve_staged_slurm_start_file(
        {},
        slurm_start_file=str(tmp_path / "example.slurm.sbatch"),
        run_id="RUN123",
    ) == str(tmp_path / "RUN123.staged.slurm.sbatch")
    sbatch_path = tmp_path / "job.sbatch"
    sbatch_path.write_text(
        "#!/bin/bash\n"
        "#SBATCH --output ../logs/slurm-%j.out\n"
        "#SBATCH -e ../logs/slurm-%j.err\n",
        encoding="utf-8",
    )
    assert _resolve_remote_sbatch_log_templates(
        _parse_sbatch_log_templates(str(sbatch_path)),
        remote_slurm_start_file="~/remote/references/job.sbatch",
    ) == {
        "output": "~/remote/logs/slurm-%j.out",
        "error": "~/remote/logs/slurm-%j.err",
    }
    status_command = _status_command("123")
    assert status_command == "scontrol show job 123 -dd"
    assert (
        _status_from_text("JobId=13590014 JobName=vhr-workflow\n   JobState=FAILED Reason=None")
        == "failed"
    )
    assert (
        _status_from_text("JobId=123 JobName=vhr-workflow\n   JobState=RUNNING Reason=None")
        == "running"
    )


@pytest.mark.parametrize("debug_logs", [False, True])
def test_slurm_upload_uses_rsync(tmp_path, monkeypatch, debug_logs):
    import vhrharmonize.slurm as slurm_mod

    calls = []

    def run(command, **kwargs):
        calls.append(command)
        return SimpleNamespace(stdout="", stderr="", returncode=0)

    monkeypatch.setattr(slurm_mod, "_run_local_command", run)
    result = slurm_mod._rsync_upload_tree(
        {
            "ssh_user": "user",
            "ssh_host": "host",
            "ssh_private_key": "~/.ssh/id_ed25519",
            "debug_logs": debug_logs,
        },
        stage_root=str(tmp_path),
        remote_root="~/remote/output",
    )
    assert result.returncode == 0
    command = calls[0]
    assert command[:3] == ["rsync", "-aL", "--itemize-changes"]
    assert ("--info=progress2" in command) is debug_logs
    assert str(Path.home() / ".ssh/id_ed25519") in command[command.index("-e") + 1]
    assert command[-2:] == [str(tmp_path) + "/", "user@host:~/remote/output/"]


def test_slurm_upload_tree_can_enable_rsync_checksum(tmp_path: Path, monkeypatch) -> None:
    import vhrharmonize.slurm as slurm_mod

    calls = []

    def fake_run_local_command(command, *, check=True, capture_output=True, stream_output=False):
        calls.append(command)
        return SimpleNamespace(stdout="", stderr="", returncode=0)

    monkeypatch.setattr(slurm_mod, "_run_local_command", fake_run_local_command)
    slurm_mod._rsync_upload_tree(
        {"ssh_user": "user", "ssh_host": "host", "enable_rsync_checksum": True},
        stage_root=str(tmp_path),
        remote_root="~/remote/output",
    )

    assert calls == [
        [
            "rsync",
            "-aL",
            "--itemize-changes",
            "--checksum",
            "--rsync-path",
            "mkdir -p ~/remote/output && rsync",
            f"{tmp_path}/",
            "user@host:~/remote/output/",
        ]
    ]


def test_slurm_upload_and_start_are_separate(tmp_path: Path, monkeypatch) -> None:
    import vhrharmonize.slurm as slurm_mod

    local_path = tmp_path / "input.tif"
    local_path.write_text("data", encoding="utf-8")
    staged_slurm = tmp_path / "RUN123.staged.hpc.yml"
    staged_slurm.write_text(
        yaml.safe_dump(
            {
                "run_id": "RUN123",
                "workflow": "vhr-workflow",
                "staged_hpc_file": str(staged_slurm),
                "ssh_host": "host",
                "ssh_user": "user",
                "remote_workflow_config": "~/remote/workflow.yml",
                "remote_slurm_start_file": "~/remote/worldview.sbatch",
                "remote_slurm_log_templates": {
                    "output": "~/remote/slurm-%j.out",
                    "error": "~/remote/slurm-%j.err",
                },
                "remote_slurm_log_paths": {},
                "uploaded_input_paths": {str(local_path): "~/remote/input.tif"},
                "uploaded_reference_paths": {},
                "download_output_paths": {},
                "download_log_paths": {},
                "debug_logs": False,
                "status": "prepared",
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        slurm_mod,
        "_upload_required_files",
        lambda slurm_data: {str(local_path): {"status": "synced"}},
    )
    uploaded = slurm_mod.upload_slurm_files(str(staged_slurm))
    assert uploaded["status"] == "uploaded"
    assert (
        slurm_mod._load_yaml_file(str(staged_slurm))["upload_results"][str(local_path)]["status"]
        == "synced"
    )

    def fail_upload(_slurm_data):
        raise AssertionError("start should not upload files")

    ssh_calls = []

    def fake_run_ssh(
        slurm_data, remote_command, *, check=True, capture_output=True, stream_output=False
    ):
        ssh_calls.append(remote_command)
        if "sbatch" in remote_command:
            return SimpleNamespace(stdout="Submitted batch job 123\n", stderr="", returncode=0)
        return SimpleNamespace(
            stdout="JobId=123 JobName=vhr-workflow\n   JobState=PENDING Reason=Priority\n",
            stderr="",
            returncode=0,
        )

    monkeypatch.setattr(slurm_mod, "_upload_required_files", fail_upload)
    monkeypatch.setattr(slurm_mod, "_run_ssh", fake_run_ssh)
    started = slurm_mod.start_slurm_job(str(staged_slurm))
    assert started["submitted_job_id"] == "123"
    assert started["status"] == "running"
    assert started["remote_slurm_log_paths"] == {
        "output": "~/remote/slurm-123.out",
        "error": "~/remote/slurm-123.err",
    }
    assert len(ssh_calls) == 2


def test_slurm_status_reads_sbatch_logs(tmp_path: Path, monkeypatch, capsys) -> None:
    import vhrharmonize.slurm as slurm_mod

    staged_slurm = tmp_path / "RUN123.staged.hpc.yml"
    staged_slurm.write_text(
        yaml.safe_dump(
            {
                "run_id": "RUN123",
                "ssh_host": "host",
                "ssh_user": "user",
                "submitted_job_id": "123",
                "remote_slurm_log_templates": {
                    "output": "~/remote/slurm-%j.out",
                    "error": "~/remote/slurm-%j.err",
                },
                "status": "submitted",
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    def fake_run_ssh(
        slurm_data, remote_command, *, check=True, capture_output=True, stream_output=False
    ):
        if "scontrol show job" in remote_command:
            return SimpleNamespace(
                stdout="JobId=123 JobName=vhr-workflow\n   JobState=COMPLETED Reason=None\n",
                stderr="",
                returncode=0,
            )
        if "slurm-123.err" in remote_command:
            return SimpleNamespace(stdout="stderr text\n", stderr="", returncode=0)
        if "slurm-123.out" in remote_command:
            return SimpleNamespace(stdout="stdout text\n", stderr="", returncode=0)
        raise AssertionError(remote_command)

    monkeypatch.setattr(slurm_mod, "_run_ssh", fake_run_ssh)
    updated = slurm_mod.update_status_slurm_file(str(staged_slurm))
    captured = capsys.readouterr().out

    assert updated["status"] == "completed"
    assert updated["remote_slurm_log_paths"] == {
        "output": "~/remote/slurm-123.out",
        "error": "~/remote/slurm-123.err",
    }
    assert updated["raw_slurm_log_text"] == {"error": "stderr text\n", "output": "stdout text\n"}
    assert "JobState=COMPLETED" in updated["raw_status_text"]
    assert "=== Slurm error log: ~/remote/slurm-123.err ===" in captured
    assert "stderr text" in captured
    assert "=== Slurm output log: ~/remote/slurm-123.out ===" in captured
    assert "stdout text" in captured
    assert "=== Slurm status ===" in captured
    assert "JobState=COMPLETED" in captured
    assert captured.index("=== Slurm output log") < captured.index("=== Slurm error log")
    assert captured.index("=== Slurm error log") < captured.index("=== Slurm status ===")


def test_slurm_stop_and_close_are_narrow(tmp_path: Path, monkeypatch) -> None:
    import pytest

    import vhrharmonize.slurm as slurm_mod

    staged_slurm = tmp_path / "RUN123.staged.hpc.yml"
    staged_slurm.write_text(
        yaml.safe_dump(
            {
                "ssh_host": "host",
                "ssh_user": "user",
                "submitted_job_id": "123",
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    ssh_calls = []
    local_calls = []

    def fake_run_ssh(
        slurm_data, remote_command, *, check=True, capture_output=True, stream_output=False
    ):
        ssh_calls.append(remote_command)
        return SimpleNamespace(stdout="", stderr="", returncode=0)

    def fake_run_local(command, *, check=True, capture_output=True, stream_output=False):
        local_calls.append(command)
        return SimpleNamespace(stdout="", stderr="", returncode=0)

    monkeypatch.setattr(slurm_mod, "_run_ssh", fake_run_ssh)
    monkeypatch.setattr(slurm_mod, "_run_local_command", fake_run_local)

    slurm_mod.stop_slurm_job(str(staged_slurm))
    slurm_mod.close_hpc_connection(str(staged_slurm))

    assert ssh_calls == ["scancel 123"]
    assert local_calls == [["ssh", "-O", "exit", "user@host"]]

    staged_slurm.write_text(
        yaml.safe_dump({"ssh_host": "host", "ssh_user": "user", "submitted_job_id": None}),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="submitted_job_id"):
        slurm_mod.stop_slurm_job(str(staged_slurm))


def test_py6s_samples_dem_using_mapped_geometry(monkeypatch, tmp_path):
    from shapely.geometry import mapping

    calls = []
    monkeypatch.setattr(atmos_mod.Py6SCorrector, "run", lambda self, **kwargs: calls.append(kwargs))

    def sample(path, percentile, mask):
        assert path == "dem.tif"
        assert percentile == 60
        assert mask.equals(box(0, 0, 1, 1))
        return 500

    monkeypatch.setattr(atmos_mod, "get_image_percentile_value", sample)
    run_py6s(
        "in.tif",
        str(tmp_path / "out.tif"),
        variables={**_metadata(), "geometry": mapping(box(0, 0, 1, 1))},
        dem_file_path="dem.tif",
        dem_ground_percentile=60,
    )
    assert calls[0]["ground_elevation_km"] == 0.5
