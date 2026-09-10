from contextlib import contextmanager
import errno
import json
import stat
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from vhrharmonize.cli import worldview
from vhrharmonize.providers.worldview import load_worldview_scenes_from_tif_files


ATMOSPHERE = {
    "source": "nasa_power",
    "date_used": "2017-07-05",
    "sample_count": 9,
    "aot550": 0.2,
    "water_vapor": 2.5,
    "ozone_cm_atm": 0.3,
}


@pytest.fixture
def atmosphere_scene(make_worldview_bundle, tmp_path):
    bundle = make_worldview_bundle()
    scene = load_worldview_scenes_from_tif_files([str(bundle["mul_tif"])])[0]
    args = worldview._build_parser().parse_args([])
    args.temp_dir = str(tmp_path / "scratch")
    args.output_dir = str(tmp_path / "output")
    args.log_to_console = True
    args.run_from_existing = True
    args.run_fetch_atmosphere = True
    args.run_atmospheric_correction = True
    args.fetch_atmosphere_source = "nasa_power"
    args.max_cloud_cover_to_process = None
    for step in (
        "file_source", "orthorectification", "pansharpen", "cloud_mask",
        "alignment", "seamline_metadata", "radiometric_normalization",
    ):
        setattr(args, f"run_{step}", False)
    state = worldview._initialize_scene_state(scene, args)
    cache = Path(state.step_dirs["fetch_atmosphere"]) / (
        f"{bundle['basename']}{args.fetch_atmosphere_output_suffix}.json"
    )
    return state, args, cache


@pytest.mark.parametrize("check_validity", [True, False])
@pytest.mark.parametrize("contents", [b"", b'{"source":', b"[]", b"\xff"])
def test_invalid_atmosphere_cache_is_refetched(
    atmosphere_scene, monkeypatch, capsys, check_validity, contents,
):
    state, args, cache = atmosphere_scene
    args.run_from_existing_check_validity = check_validity
    cache.write_bytes(contents)
    fetch = Mock(return_value=SimpleNamespace(**ATMOSPHERE))
    monkeypatch.setattr(worldview, "fetch_power_atmosphere_for_bbox", fetch)

    worldview._run_fetch_atmosphere_step(state, args)

    fetch.assert_called_once()
    assert json.loads(cache.read_text()) == ATMOSPHERE
    assert state.fetch_atmosphere_result == ATMOSPHERE
    assert state.scene.step_outputs["fetch_atmosphere"] == [str(cache)]
    output = capsys.readouterr().out
    assert "Existing output invalid; rerunning" in output
    assert "Skipping because output exists" not in output


@pytest.mark.parametrize("check_validity", [True, False])
def test_valid_atmosphere_cache_is_reused(
    atmosphere_scene, monkeypatch, capsys, check_validity,
):
    state, args, cache = atmosphere_scene
    args.run_from_existing_check_validity = check_validity
    contents = json.dumps(ATMOSPHERE)
    cache.write_text(contents)
    fetch = Mock(side_effect=AssertionError("valid cache must not be fetched"))
    bounds = Mock(side_effect=AssertionError("valid cache must not require bounds"))
    write = Mock(side_effect=AssertionError("valid cache must not be rewritten"))
    monkeypatch.setattr(worldview, "fetch_power_atmosphere_for_bbox", fetch)
    monkeypatch.setattr(worldview, "materialize_scene_bounds", bounds)
    monkeypatch.setattr(worldview, "_write_json", write)

    worldview._run_fetch_atmosphere_step(state, args)

    fetch.assert_not_called()
    bounds.assert_not_called()
    write.assert_not_called()
    assert cache.read_text() == contents
    assert state.fetch_atmosphere_result == ATMOSPHERE
    assert state.scene.step_outputs["fetch_atmosphere"] == [str(cache)]
    assert "Skipping because output exists" in capsys.readouterr().out


@pytest.mark.parametrize("contents, loaded", [
    (b"", 0),
    (b'{"source":', 0),
    (b"[]", 0),
    (json.dumps(ATMOSPHERE).encode(), 1),
])
def test_processing_counts_validate_atmosphere_json(atmosphere_scene, contents, loaded):
    state, args, cache = atmosphere_scene
    args.run_from_existing_check_validity = True
    cache.write_bytes(contents)

    counts = worldview._count_processing_steps([state.scene], args)

    assert counts["fetch_atmosphere"] == {"loaded": loaded, "processing": 1 - loaded}


@pytest.mark.parametrize("existing", [True, False])
@pytest.mark.parametrize("relative", [True, False])
def test_json_write_publishes_complete_document(tmp_path, monkeypatch, existing, relative):
    monkeypatch.chdir(tmp_path)
    path = Path("atmosphere.json") if relative else tmp_path / "nested" / "atmosphere.json"
    if existing:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('{"source": "previous"}')

    worldview._write_json(str(path), ATMOSPHERE)

    assert json.loads(path.read_text()) == ATMOSPHERE
    assert list(path.parent.iterdir()) == [path]


@pytest.mark.parametrize("existing", [True, False])
@pytest.mark.parametrize("failure", ["write", "flush", "close", "fsync", "replace"])
def test_failed_json_write_preserves_destination_and_cleans_temp(
    tmp_path, monkeypatch, existing, failure,
):
    path = tmp_path / "atmosphere.json"
    previous = b'{"source": "previous"}'
    if existing:
        path.write_bytes(previous)
    disk_full = OSError(errno.ENOSPC, "No space left on device")

    if failure == "write":
        def interrupted_dump(payload, handle, **kwargs):
            handle.write('{"source":')
            raise disk_full
        monkeypatch.setattr(worldview.json, "dump", interrupted_dump)
    elif failure in {"flush", "close"}:
        @contextmanager
        def failing_open(*args, **kwargs):
            with open(*args, **kwargs) as handle:
                if failure == "flush":
                    monkeypatch.setattr(handle, "flush", Mock(side_effect=disk_full))
                yield handle
            # Simulate an error reported only when the stream is closed.
            if failure == "close":
                raise disk_full

        monkeypatch.setattr(worldview, "open", failing_open, raising=False)
    else:
        monkeypatch.setattr(worldview.os, failure, Mock(side_effect=disk_full))

    with pytest.raises(OSError) as caught:
        worldview._write_json(str(path), ATMOSPHERE)

    assert caught.value.errno == errno.ENOSPC
    if existing:
        assert path.read_bytes() == previous
    else:
        assert not path.exists()
    assert list(tmp_path.iterdir()) == ([path] if existing else [])


@pytest.mark.parametrize("existing", [True, False])
def test_atomic_json_write_preserves_file_permissions(tmp_path, existing):
    reference = tmp_path / "reference.json"
    reference.write_text("{}")
    path = tmp_path / "atmosphere.json"
    if existing:
        path.write_text("{}")
        path.chmod(0o640)
    expected_mode = stat.S_IMODE((path if existing else reference).stat().st_mode)

    worldview._write_json(str(path), ATMOSPHERE)

    assert stat.S_IMODE(path.stat().st_mode) == expected_mode
