import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from vhrharmonize.cli import worldview
from vhrharmonize.providers.worldview import load_worldview_scenes_from_tif_files


@pytest.fixture
def cleanup_scene(make_worldview_bundle, tmp_path):
    bundle = make_worldview_bundle()
    scene = load_worldview_scenes_from_tif_files([str(bundle["mul_tif"]), str(bundle["pan_tif"])])[0]
    args = worldview._build_parser().parse_args(["--delete-temp-steps-proactively"])
    args.temp_dir = str(tmp_path / "temp")
    args.output_dir = str(tmp_path / "output")
    args.delete_temp_dir = False  # Retain the directory while allowing individual files to be cleaned.
    args.run_file_source = True
    args.save_fetch_atmosphere = "$temp/atmosphere"
    args.save_atmospheric_correction = "$temp/correction"
    args.save_orthorectification = "$temp/ortho"
    args.save_pansharpen = "$temp/pansharpen"
    args.save_cloud_mask = "$output"
    args.max_cloud_cover_to_process = None
    state = worldview._initialize_scene_state(scene, args)
    return args, state, bundle


def write_outputs(paths, make_test_raster):
    for path in paths:
        if Path(path).suffix.lower() == ".json":
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            Path(path).write_text("{}")
        else:
            make_test_raster(Path(path), count=8)


def save_scene_outputs(state, args, make_test_raster):
    paths = worldview._scene_skip_required_outputs(state, args)
    write_outputs(paths, make_test_raster)
    return paths


@pytest.mark.parametrize("enabled", [False, True])
@pytest.mark.parametrize("delete_temp", [False, True])
def test_skipped_scene_cleans_partial_temp_steps_and_sidecars(
    cleanup_scene, make_test_raster, enabled, delete_temp, capsys,
):
    args, state, bundle = cleanup_scene
    args.delete_temp_steps_proactively = enabled
    args.delete_temp_dir = delete_temp
    args.log_to_console = True
    permanent = save_scene_outputs(state, args, make_test_raster)
    expected = state.cleanup_step_outputs
    # Only some intermediate files remain; no temp outputs are registered on the scene.
    leftovers = [expected["fetch_atmosphere"][0], expected["atmospheric_correction"][0],
                 expected["orthorectification_pan"][0], *expected["file_source"]]
    for path in leftovers:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text("old temporary output")
    missing_raster = Path(expected["pansharpen"][0])
    leftovers.extend([str(missing_raster) + ".ovr", str(missing_raster.with_suffix(".hdr")),
                      expected["atmospheric_correction"][0] + ".params.txt"])
    for path in leftovers[-3:]:
        Path(path).write_text("orphaned sidecar")
    dem = Path(state.step_dirs["temp_root"]) / "dem" / f'{bundle["basename"]}_dem.tif'
    make_test_raster(dem)
    leftovers.append(str(dem))
    unrelated = make_test_raster(missing_raster.parent / "another_scene.tif")
    source_bytes = {path: path.read_bytes() for path in (bundle["mul_tif"], bundle["pan_tif"], bundle["mul_imd"])}

    result = worldview._process_scene(state.scene, args)

    assert result.current_step == "cloud_mask"
    assert all(Path(path).exists() is (not enabled) for path in leftovers)
    assert all(Path(path).exists() for path in permanent)
    assert unrelated.exists()
    assert Path(state.step_dirs["temp_root"]).is_dir()
    assert all(path.read_bytes() == contents for path, contents in source_bytes.items())
    assert "Reason: desired scene level outputs exist" in capsys.readouterr().out


@pytest.mark.parametrize("bad_step", ["cloud_mask", "cloud_mask_mask", "fetch_atmosphere"])
@pytest.mark.parametrize("failure", ["missing", "invalid"])
def test_incomplete_saved_outputs_preserve_all_temp_steps(
    cleanup_scene, make_test_raster, bad_step, failure,
):
    args, state, _ = cleanup_scene
    args.save_fetch_atmosphere = "$output"
    state = worldview._initialize_scene_state(state.scene, args)
    save_scene_outputs(state, args, make_test_raster)
    bad_path = Path(state.cleanup_step_outputs[bad_step][0])
    if failure == "missing":
        bad_path.unlink()
    else:
        bad_path.write_bytes(b"invalid")
    temp = Path(state.cleanup_step_outputs["atmospheric_correction"][0])
    temp.write_bytes(b"partial temp file")

    worldview._cleanup_completed_scene_temp_steps(state, args)

    assert temp.read_bytes() == b"partial temp file"
    assert bad_path.exists() is (failure != "missing")  # Completion inspection is read-only.


@pytest.mark.parametrize("skip_existing", [False, True])
def test_preflight_cleans_every_completed_scene_before_processing_starts(
    cleanup_scene, make_test_raster, make_worldview_bundle, monkeypatch, skip_existing,
):
    args, state, _ = cleanup_scene
    args.skip_existing = skip_existing
    other = make_worldview_bundle(
        basename="17JUL05211635-M1BS-016445286010_01_P002",
        pan_basename="17JUL05211635-P1BS-016445286010_01_P002",
    )
    second = load_worldview_scenes_from_tif_files([str(other["mul_tif"]), str(other["pan_tif"])])[0]
    states = [state, worldview._initialize_scene_state(second, args)]
    leftovers = []
    for item in states:
        save_scene_outputs(item, args, make_test_raster)
        leftovers.append(make_test_raster(Path(item.cleanup_step_outputs["pansharpen"][0])))

    def fail_first_scene(scene, args):
        assert not any(path.exists() for path in leftovers)
        raise RuntimeError("processing interrupted")

    monkeypatch.setattr(worldview, "_process_scene", fail_first_scene)
    with pytest.raises(RuntimeError, match="processing interrupted"):
        worldview._process_scenes([item.scene for item in states], args)


def test_newly_processed_scene_cleans_only_after_saved_outputs_exist(
    cleanup_scene, make_test_raster, monkeypatch,
):
    args, state, bundle = cleanup_scene
    args.run_file_source = False
    args.run_fetch_atmosphere = False
    args.run_atmospheric_correction = False
    args.run_orthorectification = False
    args.run_pansharpen = False
    args.run_alignment = True
    args.alignment_fixed_image = str(bundle["pan_tif"])
    args.save_cloud_mask = "$temp/cloud"
    args.save_alignment = "$output"
    state = worldview._initialize_scene_state(state.scene, args)
    temp_raster = Path(state.cleanup_step_outputs["cloud_mask"][0])
    temp_mask = Path(state.cleanup_step_outputs["cloud_mask_mask"][0])

    def cloudmask(**kwargs):
        write_outputs([kwargs["output_raster_path"], kwargs["output_mask_path"]], make_test_raster)
        return SimpleNamespace(mask_pixel_count=0, output_mask_path=kwargs["output_mask_path"])

    def align(**kwargs):
        assert temp_raster.exists() and temp_mask.exists()
        make_test_raster(Path(kwargs["output_image_path"]))
        return SimpleNamespace(output_image_path=kwargs["output_image_path"])

    monkeypatch.setattr(worldview, "cloudmask_raster", cloudmask)
    monkeypatch.setattr(worldview, "align_image_pair", align)
    result = worldview._process_scene(state.scene, args)

    assert not temp_raster.exists() and not temp_mask.exists()
    assert Path(result.current_files[0]).exists()
    assert Path(result.scene.metadata_report_path).exists()


@pytest.mark.parametrize("alias", ["same_path", "symlink", "hardlink"])
def test_proactive_cleanup_protects_source_aliases(cleanup_scene, make_test_raster, alias):
    args, state, bundle = cleanup_scene
    save_scene_outputs(state, args, make_test_raster)
    source = bundle["mul_tif"]
    before = source.read_bytes()
    target = Path(state.cleanup_step_outputs["atmospheric_correction"][0])
    if alias == "same_path":
        state.cleanup_step_outputs["atmospheric_correction"] = [str(source)]
    elif alias == "symlink":
        target.symlink_to(source)
    else:
        target.hardlink_to(source)
    sidecar = Path(str(source if alias == "same_path" else target) + ".ovr")
    sidecar.write_bytes(b"source overview")

    worldview._cleanup_completed_scene_temp_steps(state, args)

    assert source.read_bytes() == before
    assert alias == "same_path" or target.exists()
    assert sidecar.exists()


@pytest.mark.parametrize("consumer", ["alignment_fixed_image", "dem_file_path", "group_by_basename", "spectralmatch_kwargs_json"])
def test_proactive_cleanup_preserves_explicit_file_inputs(cleanup_scene, make_test_raster, consumer):
    args, state, _ = cleanup_scene
    save_scene_outputs(state, args, make_test_raster)
    temp = make_test_raster(Path(state.cleanup_step_outputs["atmospheric_correction"][0]))
    value = {"mosaic.tif": [f"file:{temp}"]} if consumer == "group_by_basename" else str(temp)
    if consumer == "spectralmatch_kwargs_json":
        value = json.dumps({"shared_input_images": [str(temp)]})
    setattr(args, consumer, value)

    worldview._cleanup_completed_scene_temp_steps(state, args)

    assert temp.exists()


@pytest.mark.parametrize("aggregate", ["seamline_metadata", "spectralmatch"])
@pytest.mark.parametrize("succeeds", [False, True])
@pytest.mark.parametrize("delete_temp", [False, True])
def test_aggregate_inputs_survive_until_consumers_finish(
    cleanup_scene, make_test_raster, monkeypatch, aggregate, succeeds, delete_temp,
):
    args, state, _ = cleanup_scene
    args.delete_temp_dir = delete_temp
    setattr(args, f"run_{aggregate}", True)
    setattr(args, f"save_{aggregate}", str(Path(args.output_dir) / "aggregate.tif"))
    args.save_fetch_atmosphere = "$output"
    args.save_cloud_mask = "$temp/cloud"
    state = worldview._initialize_scene_state(state.scene, args)
    permanent = save_scene_outputs(state, args, make_test_raster)
    final = make_test_raster(Path(state.cleanup_step_outputs["final_raster"][0]))
    overview = Path(str(final) + ".ovr")
    overview.write_bytes(b"overview")
    intermediate = make_test_raster(Path(state.cleanup_step_outputs["pansharpen"][0]))
    monkeypatch.setattr(worldview, "_collect_input_files_by_stage", lambda value: {"file_source": state.source_files})
    monkeypatch.setattr(worldview, "_load_worldview_scenes_from_stage_paths", lambda *a, **k: [state.scene])

    def consume(*a, **k):
        assert final.exists() and overview.exists()
        assert not intermediate.exists()
        if not succeeds:
            raise RuntimeError("aggregate failed")
        return str(make_test_raster(Path(getattr(args, f"save_{aggregate}"))))

    monkeypatch.setattr(worldview, f"_run_{aggregate}_workflow", consume)
    if succeeds:
        assert worldview._run_workflow(args) == 0
    else:
        with pytest.raises(RuntimeError, match="aggregate failed"):
            worldview._run_workflow(args)

    assert final.exists() is (not succeeds)
    assert overview.exists() is (not succeeds)
    assert all(Path(path).exists() for path in permanent)
    assert Path(args.temp_dir).exists() is (not delete_temp or not succeeds)


def test_temp_only_scene_is_retained_without_a_saved_aggregate(cleanup_scene, make_test_raster):
    args, state, _ = cleanup_scene
    args.save_cloud_mask = "$temp/cloud"
    state = worldview._initialize_scene_state(state.scene, args)
    final = make_test_raster(Path(state.cleanup_step_outputs["final_raster"][0]))

    worldview._cleanup_completed_scene_temp_steps(state, args)

    assert final.exists()


def test_temp_only_scene_can_be_cleaned_after_a_persistent_aggregate(cleanup_scene, make_test_raster):
    args, state, _ = cleanup_scene
    args.save_cloud_mask = "$temp/cloud"
    args.run_spectralmatch = True
    state = worldview._initialize_scene_state(state.scene, args)
    final = make_test_raster(Path(state.cleanup_step_outputs["final_raster"][0]))
    mosaic = make_test_raster(Path(args.output_dir) / "mosaic.tif")

    worldview._cleanup_completed_scene_temp_steps(state, args, aggregate_outputs=[str(mosaic)])

    assert not final.exists()
    assert mosaic.exists()


@pytest.mark.parametrize("succeeds", [False, True])
def test_spectralmatch_success_controls_cleanup_without_output_validation(
    cleanup_scene, make_test_raster, monkeypatch, succeeds,
):
    args, state, _ = cleanup_scene
    args.run_spectralmatch = True
    args.save_spectralmatch = str(Path(args.output_dir) / "timelapse")
    args.save_cloud_mask = "$temp/cloud"
    state = worldview._initialize_scene_state(state.scene, args)
    final = make_test_raster(Path(state.cleanup_step_outputs["final_raster"][0]))
    make_test_raster(Path(state.cleanup_step_outputs["cloud_mask_mask"][0]))
    monkeypatch.setattr(worldview, "_collect_input_files_by_stage", lambda value: {"file_source": state.source_files})
    monkeypatch.setattr(worldview, "_load_worldview_scenes_from_stage_paths", lambda *a, **k: [state.scene])
    inspect = worldview._existing_outputs_are_reusable

    def check(paths, **kwargs):
        assert args.save_spectralmatch not in paths
        return inspect(paths, **kwargs)

    def consume(*a, **k):
        assert final.exists()
        if not succeeds:
            raise RuntimeError("SpectralMatch failed")
        # Successful return is the completion signal; the wrapper must not inspect it.
        return args.save_spectralmatch

    monkeypatch.setattr(worldview, "_existing_outputs_are_reusable", check)
    monkeypatch.setattr(worldview, "_run_spectralmatch_workflow", consume)
    if succeeds:
        assert worldview._run_workflow(args) == 0
    else:
        with pytest.raises(RuntimeError, match="SpectralMatch failed"):
            worldview._run_workflow(args)

    assert final.exists() is (not succeeds)


def test_staged_permanent_source_bundle_must_be_complete(cleanup_scene, make_test_raster):
    args, state, _ = cleanup_scene
    args.save_file_source = "$output/source"
    # The workflow loader registers raw files under file_source before staging.
    state.scene.step_outputs["file_source"] = [str(state.scene.mul_image.tif_file)]
    state = worldview._initialize_scene_state(state.scene, args)
    write_outputs(state.cleanup_step_outputs["cloud_mask"] + state.cleanup_step_outputs["cloud_mask_mask"], make_test_raster)
    temp = make_test_raster(Path(state.cleanup_step_outputs["pansharpen"][0]))

    worldview._cleanup_completed_scene_temp_steps(state, args)

    assert temp.exists()  # Raw inputs cannot stand in for the requested staged outputs.


@pytest.mark.parametrize("outputs_complete", [False, True])
def test_cloud_filtered_scenes_only_clean_when_saved_outputs_are_complete(
    cleanup_scene, make_test_raster, outputs_complete,
):
    args, state, _ = cleanup_scene
    args.max_cloud_cover_to_process = 10
    if outputs_complete:
        save_scene_outputs(state, args, make_test_raster)
    temp = make_test_raster(Path(state.cleanup_step_outputs["pansharpen"][0]))

    worldview._process_scene(state.scene, args)

    assert temp.exists() is (not outputs_complete)


def test_proactive_option_yaml_and_cli_override(tmp_path):
    config = tmp_path / "config.yml"
    config.write_text("shared:\n  delete_temp_steps_proactively: true\n")
    parser = worldview._build_parser()
    assert not parser.parse_args([]).delete_temp_steps_proactively
    parser.set_defaults(**worldview._normalize_config_defaults(worldview._load_worldview_yaml_config(str(config))))
    assert parser.parse_args([]).delete_temp_steps_proactively
    assert not parser.parse_args(["--no-delete-temp-steps-proactively"]).delete_temp_steps_proactively


def test_directory_cleanup_option_yaml_and_cli_override(tmp_path):
    config = tmp_path / "config.yml"
    config.write_text("shared:\n  delete_temp_dir: false\n")
    parser = worldview._build_parser()
    assert parser.parse_args([]).delete_temp_dir
    assert not parser.parse_args(["--no-delete-temp-dir"]).delete_temp_dir
    parser.set_defaults(**worldview._normalize_config_defaults(worldview._load_worldview_yaml_config(str(config))))
    args = parser.parse_args([])
    assert not args.delete_temp_dir
    assert "keep_temp_dir" not in vars(args)
    assert parser.parse_args(["--delete-temp-dir"]).delete_temp_dir


@pytest.mark.parametrize("key", ["keep_temp_dir", "keep-temp-dir"])
@pytest.mark.parametrize("value", ["true", "false"])
def test_removed_directory_cleanup_config_is_rejected(tmp_path, key, value):
    config = tmp_path / "config.yml"
    config.write_text(f"shared:\n  {key}: {value}\n")
    with pytest.raises(ValueError, match="Unrecognized config key: keep_temp_dir"):
        worldview.main(["--config-yaml", str(config)])


@pytest.mark.parametrize("flag", ["--keep-temp-dir", "--no-keep-temp-dir"])
def test_removed_directory_cleanup_flags_are_rejected(flag):
    with pytest.raises(SystemExit, match=f"Unrecognized argument: {flag}"):
        worldview.main([flag])


@pytest.mark.parametrize("delete_temp", [False, True])
@pytest.mark.parametrize("match_override", [None, False, True])
def test_spectralmatch_cleanup_defaults_follow_workflow_flag(delete_temp, match_override):
    args = worldview._build_parser().parse_args([])
    args.delete_temp_dir = delete_temp
    args.match_delete_temp_dir = match_override

    kwargs = worldview._build_spectralmatch_kwargs(args)

    assert kwargs["delete_temp_dir"] is (delete_temp if match_override is None else match_override)


@pytest.mark.parametrize("delete_temp", [False, True])
@pytest.mark.parametrize("proactive", [False, True])
def test_directory_and_per_file_cleanup_are_independent(
    cleanup_scene, make_test_raster, monkeypatch, delete_temp, proactive,
):
    args, state, bundle = cleanup_scene
    args.delete_temp_dir = delete_temp
    args.delete_temp_steps_proactively = proactive
    args.run_seamline_metadata = True
    args.save_seamline_metadata = str(Path(args.output_dir) / "footprints.gpkg")
    permanent = save_scene_outputs(state, args, make_test_raster)
    intermediate = make_test_raster(Path(state.cleanup_step_outputs["pansharpen"][0]))
    scratch = Path(args.temp_dir) / "unregistered_scratch.txt"
    scratch.write_text("workflow scratch")
    monkeypatch.setattr(worldview, "_collect_input_files_by_stage", lambda value: {"file_source": state.source_files})
    monkeypatch.setattr(worldview, "_load_worldview_scenes_from_stage_paths", lambda *a, **k: [state.scene])

    def aggregate(*a, **k):
        assert intermediate.exists() is (not proactive)
        assert scratch.exists()  # Only directory cleanup removes unregistered scratch.
        assert args.delete_temp_dir is delete_temp
        output = Path(args.save_seamline_metadata)
        output.write_text("completed footprints")
        return str(output)

    monkeypatch.setattr(worldview, "_run_seamline_metadata_workflow", aggregate)
    assert worldview._run_workflow(args) == 0

    assert Path(args.temp_dir).exists() is (not delete_temp)
    assert scratch.exists() is (not delete_temp)
    assert intermediate.exists() is (not delete_temp and not proactive)
    assert all(Path(path).exists() for path in permanent)
    assert bundle["mul_tif"].exists() and bundle["mul_imd"].exists()


@pytest.mark.parametrize("protected_kind", ["source", "configured_input", "saved_output", "source_symlink"])
def test_directory_cleanup_preserves_roots_containing_protected_files(
    cleanup_scene, make_test_raster, protected_kind,
):
    args, state, bundle = cleanup_scene
    args.delete_temp_dir = True
    args.delete_temp_steps_proactively = False
    if protected_kind == "source":
        args.temp_dir = str(bundle["scene_root"])
    elif protected_kind == "configured_input":
        args.alignment_fixed_image = str(make_test_raster(Path(args.temp_dir) / "reference.tif"))
    elif protected_kind == "saved_output":
        args.output_dir = str(Path(args.temp_dir) / "saved")
    state = worldview._initialize_scene_state(state.scene, args)
    if protected_kind == "source_symlink":
        source_link = Path(args.temp_dir) / "original.tif"
        source_link.symlink_to(bundle["mul_tif"])
        state.source_files.append(str(source_link))
    permanent = save_scene_outputs(state, args, make_test_raster)

    worldview._cleanup_workflow_temp_dirs([state], args, [])

    assert Path(args.temp_dir).is_dir()
    assert all(Path(path).exists() for path in permanent)
    assert all(Path(path).exists() for path in state.source_files)
    if protected_kind == "configured_input":
        assert Path(args.alignment_fixed_image).exists()


def test_directory_cleanup_preserves_parent_of_working_directory(cleanup_scene, monkeypatch):
    args, state, _ = cleanup_scene
    args.delete_temp_dir = True
    working_dir = Path(args.temp_dir) / "active_work"
    working_dir.mkdir()
    monkeypatch.chdir(working_dir)

    worldview._cleanup_workflow_temp_dirs([state], args, [])

    assert working_dir.is_dir()
