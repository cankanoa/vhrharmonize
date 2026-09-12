from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

from osgeo import gdal
import pytest
import rasterio

from vhrharmonize.cli import worldview
from vhrharmonize.providers.worldview import load_worldview_scenes_from_tif_files


@pytest.mark.parametrize("use_exceptions", [False, True])
def test_raster_validity_handles_unreadable_tiff_directory(tmp_path, make_test_raster, use_exceptions):
    corrupt_path = tmp_path / "corrupt.tif"
    # TIFF header pointing to a directory that was never written.
    corrupt_path.write_bytes(b"II\x2a\x00\x08\x02\x00\x00")
    valid_path = make_test_raster(tmp_path / "valid.tif")

    with gdal.ExceptionMgr(useExceptions=use_exceptions):
        is_valid, reason = worldview._gdal_raster_is_valid(str(corrupt_path))
        assert not is_valid
        assert reason
        assert worldview._gdal_raster_is_valid(str(valid_path)) == (True, None)


@pytest.mark.parametrize("skip_existing", [False, True])
def test_processing_counts_schedule_corrupt_raster_for_regeneration(make_worldview_bundle, tmp_path, skip_existing):
    bundle = make_worldview_bundle()
    scene = load_worldview_scenes_from_tif_files([str(bundle["mul_tif"])])[0]
    args = worldview._build_parser().parse_args([])
    args.temp_dir = str(tmp_path / "temp")
    args.output_dir = str(tmp_path / "output")
    args.skip_existing = skip_existing
    args.max_cloud_cover_to_process = None
    args.run_fetch_atmosphere = False
    args.run_orthorectification = False
    args.run_pansharpen = False
    args.run_cloud_mask = False
    state = worldview._initialize_scene_state(scene, args)
    output_path = Path(worldview._get_expected_scene_step_outputs(state, args)["atmospheric_correction"][0])
    output_path.write_bytes(b"II\x2a\x00\x08\x02\x00\x00")

    with gdal.ExceptionMgr(useExceptions=True):
        counts = worldview._count_processing_steps([scene], args)

    assert counts["atmospheric_correction"] == {"loaded": 0, "processing": 1}
    assert output_path.exists()  # Planning must not clean outputs.


@pytest.mark.parametrize("reuse", [False, True])
@pytest.mark.parametrize("check_validity", [False, True])
def test_prepare_outputs_cleans_all_invalid_files_in_partial_sets(tmp_path, make_test_raster, reuse, check_validity):
    args = worldview._build_parser().parse_args([])
    args.run_from_existing = reuse
    args.run_from_existing_check_validity = check_validity
    valid = make_test_raster(tmp_path / "valid.tif")
    valid_bytes = valid.read_bytes()
    invalid = tmp_path / "invalid.tif"
    invalid.write_bytes(b"II\x2a\x00\x08\x02\x00\x00")
    invalid_json = tmp_path / "invalid.json"
    invalid_json.write_text('{"water_vapor":')
    sidecars = [Path(f"{invalid}{suffix}") for suffix in (".aux.xml", ".msk", ".ovr")]
    for path in sidecars:
        path.write_text("stale")
    paths = [str(tmp_path / "missing.tif"), str(valid), str(invalid), str(invalid_json)]

    with gdal.ExceptionMgr(useExceptions=True):
        assert not worldview._existing_outputs_are_reusable(
            paths, check_validity=check_validity, validity_check_grid_size=0,
            log_to_console=False, step="cloud_mask",
        )
        assert all(path.exists() for path in [invalid, invalid_json, *sidecars])
        assert not worldview._prepare_step_outputs(paths, input_paths=[], args=args, step="cloud_mask")

    assert valid.read_bytes() == valid_bytes
    assert all(path.exists() is (not check_validity) for path in [invalid, invalid_json])
    assert all(path.read_text() == "stale" for path in sidecars)


@pytest.mark.parametrize("reuse", [False, True])
def test_prepare_outputs_preserves_valid_files(tmp_path, make_test_raster, reuse):
    args = worldview._build_parser().parse_args([])
    args.run_from_existing = reuse
    output = make_test_raster(tmp_path / "valid.tif")
    before = output.read_bytes()
    assert worldview._prepare_step_outputs([str(output)], input_paths=[], args=args, step="alignment") is reuse
    assert output.read_bytes() == before


@pytest.mark.parametrize("named", [False, True])
def test_spectralmatch_outputs_are_left_to_pipeline(tmp_path, make_test_raster, monkeypatch, named):
    args = worldview._build_parser().parse_args([])
    source = make_test_raster(tmp_path / "source.tif")
    target = tmp_path / "result.tif"
    target.write_bytes(b"II\x2a\x00\x08\x02\x00\x00")

    def regenerate(**kwargs):
        assert target.read_bytes() == b"II\x2a\x00\x08\x02\x00\x00"
        target.unlink()  # The pipeline owns replacement of its invalid outputs.
        make_test_raster(target)
        return str(target)

    monkeypatch.setattr(worldview, "spectralmatch", regenerate)
    with gdal.ExceptionMgr(useExceptions=True):
        if named:
            result = worldview._run_named_spectralmatch_group(
                "result.tif", "auto:*.tif", available_paths=[str(source)], args=args,
                temp_root=str(tmp_path), output_root=str(tmp_path),
            )
        else:
            result = worldview._run_default_spectralmatch(
                [str(source)], args=args, output_path=str(target), temp_root=str(tmp_path),
            )
    assert result == str(target)
    assert target.exists()


def test_all_scene_steps_regenerate_corrupt_outputs(make_worldview_bundle, tmp_path, make_test_raster, monkeypatch):
    bundle = make_worldview_bundle()
    source_paths = [bundle["mul_tif"], bundle["pan_tif"]]
    originals = {path: path.read_bytes() for path in source_paths}
    scene = load_worldview_scenes_from_tif_files([str(path) for path in source_paths])[0]
    args = worldview._build_parser().parse_args([])
    args.temp_dir = str(tmp_path / "temp")
    args.output_dir = str(tmp_path / "output")
    args.run_file_source = True
    args.run_alignment = True
    args.alignment_fixed_image = str(bundle["pan_tif"])
    args.dem_file_path = str(bundle["pan_tif"])
    state = worldview._initialize_scene_state(scene, args)
    expected = worldview._get_expected_scene_step_outputs(state, args)
    corrupted = set()
    for step, paths in expected.items():
        if step == "raw":
            continue
        for path in paths:
            if Path(path).suffix.lower() not in {".tif", ".json"}:
                continue
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            Path(path).write_bytes(b"II\x2a\x00\x08\x02\x00\x00")
            corrupted.add(path)

    def regenerate(path):
        assert not Path(path).exists(), f"Step did not clean corrupt output: {path}"
        make_test_raster(Path(path), count=8)

    def correct(**kwargs):
        regenerate(kwargs["output_raster"])
        return SimpleNamespace(effective_params={}, auto_atmos_estimate=None)

    def cloudmask(**kwargs):
        regenerate(kwargs["output_raster_path"])
        regenerate(kwargs["output_mask_path"])
        return SimpleNamespace(mask_pixel_count=0, output_mask_path=kwargs["output_mask_path"])

    def align(**kwargs):
        regenerate(kwargs["output_image_path"])
        return SimpleNamespace(output_image_path=kwargs["output_image_path"])

    monkeypatch.setattr(worldview, "fetch_power_atmosphere_for_bbox", lambda **kwargs: SimpleNamespace(
        source="nasa_power", date_used="2017-07-05", sample_count=1,
        aot550=0.2, water_vapor=2.5, ozone_cm_atm=0.3,
    ))
    monkeypatch.setattr(worldview, "get_image_percentile_value", lambda *args, **kwargs: 1000.0)
    monkeypatch.setattr(worldview, "run_py6s", correct)
    monkeypatch.setattr(worldview, "gcp_refined_rpc_orthorectification", lambda src, dst, *args, **kwargs: regenerate(dst))
    monkeypatch.setattr(worldview, "pansharpen_image", lambda src, pan, dst, **kwargs: regenerate(dst))
    monkeypatch.setattr(worldview, "cloudmask_raster", cloudmask)
    monkeypatch.setattr(worldview, "align_image_pair", align)

    with gdal.ExceptionMgr(useExceptions=True):
        for step in ("file_source", "fetch_atmosphere", "atmospheric_correction", "orthorectification", "pansharpen", "cloud_mask", "alignment"):
            getattr(worldview, f"_run_{step}_step")(state, args)
        assert worldview._existing_outputs_are_reusable(
            list(corrupted), check_validity=True, validity_check_grid_size=0,
            log_to_console=False, step="workflow",
        )
    assert all(path.read_bytes() == contents for path, contents in originals.items())


def test_py6s_regenerates_corrupt_output_after_shared_preparation(tmp_path, make_test_raster, monkeypatch):
    from vhrharmonize.preprocess.atmospheric_correction import Py6SCorrector, SixS

    source = make_test_raster(tmp_path / "source.tif")
    target = tmp_path / "corrected.tif"
    target.write_bytes(b"II\x2a\x00\x08\x02\x00\x00")
    args = worldview._build_parser().parse_args([])
    def coefficients(sixs):
        sixs.outputs = SimpleNamespace(coef_xa=0.1, coef_xb=0.0, coef_xc=0.0)
    monkeypatch.setattr(SixS, "run", coefficients)

    with gdal.ExceptionMgr(useExceptions=True):
        assert not worldview._prepare_step_outputs(
            [str(target)], input_paths=[str(source)], args=args, step="atmospheric_correction",
        )
    Py6SCorrector().run(
        str(source), str(target), solar_zenith=17, solar_azimuth=75,
        view_zenith=22, view_azimuth=80, day=5, month=7, sixs_executable="unused-sixs",
    )
    with rasterio.open(target) as dst:
        assert dst.read(1)[0, 0] == pytest.approx(0.1)


@pytest.mark.parametrize("check_validity", [False, True])
@pytest.mark.parametrize(
    "contents, valid",
    [
        (b'{"water_vapor": 2.5}', True),
        (b"{}", True),
        (b"[]", True),
        (b"null", True),
        (b"", False),
        (b'{"water_vapor":', False),
        (b"{} trailing text", False),
        (b'"\xff"', False),
    ],
)
def test_json_output_reuse_checks_syntax_only_when_enabled(tmp_path, contents, valid, check_validity):
    output_path = tmp_path / "atmosphere.JSON"
    output_path.write_bytes(contents)

    assert worldview._existing_outputs_are_reusable(
        [str(output_path)],
        check_validity=check_validity,
        validity_check_grid_size=0,
        log_to_console=False,
        step="fetch_atmosphere",
    ) is (valid or not check_validity)


@pytest.mark.parametrize("contents, should_fetch", [("{}", False), ('{"water_vapor":', True)])
def test_fetch_atmosphere_refetches_invalid_json(make_worldview_bundle, tmp_path, monkeypatch, contents, should_fetch):
    bundle = make_worldview_bundle()
    scene = load_worldview_scenes_from_tif_files([str(bundle["mul_tif"])])[0]
    args = worldview._build_parser().parse_args([])
    args.run_fetch_atmosphere = True
    args.run_from_existing = True
    args.run_from_existing_check_validity = True
    args.fetch_atmosphere_source = "nasa_power"
    output_dir = tmp_path / "fetch_atmosphere"
    output_dir.mkdir()
    output_path = output_dir / f'{bundle["basename"]}{args.fetch_atmosphere_output_suffix}.json'
    output_path.write_text(contents, encoding="utf-8")
    state = worldview.SceneWorkflowState(
        scene=scene,
        step_dirs={"fetch_atmosphere": str(output_dir)},
        current_files=[str(bundle["mul_tif"])],
    )
    result = dict(source="nasa_power", date_used="2017-07-05", sample_count=1,
                  aot550=0.2, water_vapor=2.5, ozone_cm_atm=0.3)
    fetch = Mock(return_value=SimpleNamespace(**result))
    monkeypatch.setattr(worldview, "fetch_power_atmosphere_for_bbox", fetch)

    worldview._run_fetch_atmosphere_step(state, args)

    assert fetch.call_count == int(should_fetch)
    assert state.fetch_atmosphere_result == (result if should_fetch else {})
    assert json.loads(output_path.read_text(encoding="utf-8")) == state.fetch_atmosphere_result
    assert scene.step_outputs["fetch_atmosphere"] == [str(output_path)]
