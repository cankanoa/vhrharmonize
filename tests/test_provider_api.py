from __future__ import annotations

from pathlib import Path

from vhrharmonize.providers.standardized import StandardizedMetadata
from vhrharmonize.providers.worldview import (
    discover_worldview_scene_tree_from_tif_files,
    enrich_worldview_scenes_with_metadata,
    find_files,
    iter_worldview_scenes,
    load_worldview_metadata,
    load_worldview_scenes_from_tif_files,
    parse_worldview_basename,
    parse_worldview_imd_file,
    parse_worldview_imd_text,
)


def test_parse_worldview_basename() -> None:
    parts = parse_worldview_basename("17JUL05211635-M1BS-016445286010_01_P001.TIF")
    assert parts is not None
    assert parts.image_role == "mul"
    assert parts.scene_id == "P001"


def test_parse_worldview_imd_text_and_file(tmp_path: Path) -> None:
    imd_text = 'version = "1";\nBEGIN_GROUP = IMAGE_1\nsatId = "WV03";\nEND_GROUP = IMAGE_1\nEND;'
    imd_path = tmp_path / "scene.IMD"
    imd_path.write_text(imd_text, encoding="utf-8")
    assert parse_worldview_imd_text(imd_text)["version"] == "1"
    assert parse_worldview_imd_file(str(imd_path))["IMAGE_1"]["satId"] == "WV03"


def test_worldview_scene_discovery_and_metadata(make_worldview_bundle) -> None:
    bundle = make_worldview_bundle()
    tif_files = [str(bundle["mul_tif"]), str(bundle["pan_tif"])]
    tree = discover_worldview_scene_tree_from_tif_files(tif_files)
    scenes = iter_worldview_scenes(tree)
    assert len(scenes) == 1
    enriched = enrich_worldview_scenes_with_metadata(scenes)
    assert enriched[0].mul_image is not None
    assert enriched[0].mul_image.standardized_metadata is not None
    loaded = load_worldview_scenes_from_tif_files(tif_files)
    assert loaded[0].primary_basename == bundle["basename"]


def test_find_files_and_standardized_metadata(make_worldview_bundle) -> None:
    bundle = make_worldview_bundle()
    found = find_files(str(bundle["scene_root"]))
    scene = next(iter(found.values()))
    metadata = load_worldview_metadata(str(bundle["mul_imd"]), photo_basename=bundle["basename"])
    standardized = StandardizedMetadata.from_worldview_metadata(metadata)
    assert scene["mul_photo_basename"] == bundle["basename"]
    assert round(float(standardized.cloud_cover or 0.0), 1) == 81.4


def test_source_bundle_uses_only_image_directory_and_raster_companions(make_worldview_bundle):
    from vhrharmonize.providers.worldview.core import _worldview_image_source_files
    from vhrharmonize.cli.worldview import _image_source_file_map
    from vhrharmonize.slurm import _worldview_file_source_upload_inputs
    from vhrharmonize.providers.standardized import materialize_scene_bounds

    bundle = make_worldview_bundle()
    image_path = bundle["mul_tif"]
    image_path.with_suffix(".RPB").write_text("RPC metadata")
    for extension in (".shp", ".shx", ".dbf", ".prj", ".cpg", ".shp.xml"):
        image_path.with_suffix(extension).write_text("unneeded vector companion")
    scene = load_worldview_scenes_from_tif_files([str(image_path)])[0]
    image = scene.mul_image
    paths = _worldview_image_source_files(image)
    assert {Path(path).suffix for path in paths} == {".TIF", ".IMD", ".RPB"}
    assert set(_image_source_file_map(image, "/staging")) == set(paths)
    uploads, inputs = _worldview_file_source_upload_inputs(scene, None, source_step="file_source", state=None)
    assert {path for step, path in uploads} == set(paths)
    assert inputs == [("file_source", str(image_path))]
    assert not hasattr(image, "shp_file")
    assert not hasattr(image, "worldview_metadata")
    assert materialize_scene_bounds(image.standardized_metadata.source_metadata).bounds == (0, 0, 4, 4)


def test_scene_initialization_and_reused_steps_do_not_materialize_bounds(make_worldview_bundle, tmp_path, monkeypatch):
    from unittest.mock import Mock
    from vhrharmonize.cli import worldview

    bundle = make_worldview_bundle()
    scene = load_worldview_scenes_from_tif_files([str(bundle["mul_tif"])])[0]
    args = worldview._build_parser().parse_args([])
    args.log_to_console = False
    args.run_atmospheric_correction = True
    args.run_fetch_atmosphere = True
    args.run_orthorectification = False
    args.run_pansharpen = False
    args.run_from_existing = True
    step_dirs = {step: str(tmp_path / step) for step in ("atmospheric_correction", "fetch_atmosphere")}
    monkeypatch.setattr(worldview, "_resolve_scene_step_dirs", lambda *args: step_dirs)
    materialize = Mock(side_effect=AssertionError("unused bounds were materialized"))
    monkeypatch.setattr(worldview, "materialize_scene_bounds", materialize)
    state = worldview._initialize_scene_state(scene, args)
    monkeypatch.setattr(worldview, "_prepare_step_outputs", lambda *args, **kwargs: True)
    monkeypatch.setattr(worldview, "_read_json", lambda path: {"source": "existing"})
    worldview._run_fetch_atmosphere_step(state, args)
    worldview._run_atmospheric_correction_step(state, args)
    materialize.assert_not_called()
    assert state.fetch_atmosphere_result == {"source": "existing"}
    assert state.current_step == "atmospheric_correction"
