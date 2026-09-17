from copy import deepcopy
from pathlib import Path

import pytest

from vhrharmonize.cli import worldview


@pytest.fixture
def partial_run(make_worldview_bundle, tmp_path):
    args = worldview._build_parser().parse_args([])
    for step in ("fetch_atmosphere", "atmospheric_correction", "orthorectification",
                 "pansharpen", "cloud_mask", "alignment", "seamline_metadata", "spectralmatch"):
        setattr(args, f"run_{step}", False)
    args.run_file_source = True
    args.skip_existing = False
    args.run_from_existing = True
    args.temp_dir = str(tmp_path / "temp")
    args.output_dir = str(tmp_path / "output")
    args.save_file_source = "$output/source"
    scenes = []
    for second in (35, 36, 37):
        bundle = make_worldview_bundle(
            basename=f"17JUL052116{second}-M1BS-016445286010_01_P001",
            pan_basename=f"17JUL052116{second}-P1BS-016445286010_01_P001",
        )
        scenes.extend(worldview.load_worldview_scenes_from_tif_files([
            str(bundle["mul_tif"]), str(bundle["pan_tif"]),
        ]))
    state = worldview._initialize_scene_state(deepcopy(scenes[1]), args)
    worldview._run_file_source_step(state, args)
    args.log_to_console = True
    return scenes, args


@pytest.mark.parametrize("workers", [1, 2])
def test_progress_counts_only_pending_scenes(partial_run, workers, capfd):
    scenes, args = partial_run
    args.concurrent_processing = workers
    counts = worldview._count_processing_steps(scenes, args)
    assert counts["file_source"] == {"loaded": 1, "processing": 2}

    states = worldview._process_scenes(scenes, args)

    lines = [line for line in capfd.readouterr().out.splitlines()
             if "file_source] Completed" in line]
    assert len(lines) == 2
    assert lines[0].endswith("Completed 1/2/3")
    assert lines[1].endswith("Completed 2/2/3")
    assert all(scenes[1].primary_basename not in line for line in lines)
    assert all(Path(state.current_files[0]).exists() for state in states)


def test_dask_workers_share_step_progress(partial_run, monkeypatch, capfd):
    distributed = pytest.importorskip("dask.distributed")
    scenes, args = partial_run
    args.concurrent_processing_backend = "dask"
    with distributed.LocalCluster(
        n_workers=2, threads_per_worker=1, processes=False,
        protocol="inproc", dashboard_address=None,
    ) as cluster:
        monkeypatch.setattr(worldview, "_make_dask_client", lambda args: distributed.Client(cluster))
        worldview._process_scenes(scenes, args)

    lines = [line for line in capfd.readouterr().out.splitlines()
             if "file_source] Completed" in line]
    assert len(lines) == 2
    assert lines[0].endswith("Completed 1/2/3")
    assert lines[1].endswith("Completed 2/2/3")


def test_failed_step_does_not_advance_progress(partial_run, monkeypatch, capsys):
    scenes, args = partial_run
    worldview._count_processing_steps(scenes, args)
    args.scene_total = len(scenes)
    args._completed_steps = {}
    state = worldview._initialize_scene_state(scenes[0], args)

    def fail(path_map):
        raise RuntimeError("copy failed")

    monkeypatch.setattr(worldview, "_copy_file_source_bundle", fail)
    with pytest.raises(RuntimeError, match="copy failed"):
        worldview._run_file_source_step(state, args)

    assert "Completed" not in capsys.readouterr().out
    assert args._completed_steps == {}


def test_completion_counts_follow_finish_order(partial_run, capsys):
    scenes, args = partial_run
    worldview._count_processing_steps(scenes, args)
    args.scene_total = len(scenes)
    args._completed_steps = {}
    for scene in reversed(scenes):
        worldview._run_file_source_step(worldview._initialize_scene_state(scene, args), args)

    lines = [line for line in capsys.readouterr().out.splitlines()
             if "file_source] Completed" in line]
    assert lines == [
        f"[{scenes[2].primary_basename} file_source] Completed 1/2/3",
        f"[{scenes[0].primary_basename} file_source] Completed 2/2/3",
    ]
