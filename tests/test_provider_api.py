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
