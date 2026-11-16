"""Raw footprints and ranked coverage use core reuse and completeness checks."""
from pathlib import Path

import pytest

from vhrharmonize.plugins.create_footprints import CreateFootprints
from vhrharmonize.workflow.engine import Workflow


def recipe(tmp_path, make_test_raster, *, scope="aggregate", validity=False):
    for name in ("a", "b"):
        make_test_raster(tmp_path / f"{name}.tif", crs="EPSG:3857")
    return {
        "shared": {"plugin": "shared", "core:run": True,
                   "core:log_to_console": False, "core:show_progress": False,
                   "core:save_statistics_path": None, "core:load_statistics_path": None,
                   "core:concurrent_processing": 1, "core:check_validity": validity},
        "input": {"plugin": "import_files", "core:run": True,
                  "param:search_glob": str(tmp_path / "*.tif"),
                  "var:image": "returned:file_path"},
        "footprints": {"plugin": "create_footprints", "core:run": True,
                       "core:scope": scope, "core:require_outputs": True,
                       "param:input_images": "collect:image" if scope == "aggregate" else "expr:[var.image]",
                       "param:output_polygons": str(tmp_path / "footprints.gpkg") if scope == "aggregate" else "expr:var.image & '.gpkg'",
                       "param:footprint_source": "calculate_bounds"},
    }


@pytest.mark.parametrize("scope", ["aggregate", "var"])
@pytest.mark.parametrize("validity", [False, "before_reuse"])
def test_completed_raw_footprints_reuse_without_calling_polygonizer(
    tmp_path, make_test_raster, monkeypatch, scope, validity,
):
    config = recipe(tmp_path, make_test_raster, scope=scope, validity=validity)
    Workflow(config).run()

    def unexpected_run(*args, **kwargs):
        pytest.fail("Core should reuse the footprint without invoking SpectralMatch")

    monkeypatch.setattr(CreateFootprints, "run", unexpected_run)
    workflow = Workflow(config)
    assert workflow.counts()["footprints"]["loaded"] == (2 if scope == "var" else 1)
    assert len(workflow.run()) == 2


@pytest.mark.parametrize("change", [
    {"param:image_rank_function": "{stat:area_crs}"},
    {"param:extrude_distance": 500},
    {"param:outline_path": "outline.gpkg"},
])
@pytest.mark.parametrize("validity", [False, "before_reuse"])
def test_coverage_reuses_without_calling_function(tmp_path, make_test_raster, monkeypatch, change, validity):
    config = recipe(tmp_path, make_test_raster, validity=validity)
    Workflow(config).run()
    config["footprints"].update(change)
    if "param:outline_path" in change:
        config["footprints"]["param:outline_path"] = str(tmp_path / "outline.gpkg")
        (tmp_path / "outline.gpkg").write_bytes((tmp_path / "footprints.gpkg").read_bytes())
    # Build real coverage once, then ensure the next run restores it directly.
    config["footprints"]["core:reuse"] = False
    retained = Workflow(config).run()
    config["footprints"].pop("core:reuse")

    def unexpected_run(*args, **kwargs):
        pytest.fail("Completed coverage should be reused without rerunning SpectralMatch")

    monkeypatch.setattr(CreateFootprints, "run", unexpected_run)
    workflow = Workflow(config)
    assert workflow.counts()["footprints"]["loaded"] == 1
    assert [r["id"] for r in workflow.run()] == [r["id"] for r in retained]


def test_explicit_rerun_does_not_reuse(tmp_path, make_test_raster):
    config = recipe(tmp_path, make_test_raster)
    config["footprints"]["param:image_rank_function"] = "{stat:area_crs}"
    Workflow(config).run()
    config["footprints"]["core:reuse"] = False
    workflow = Workflow(config)
    assert workflow.counts()["footprints"]["processing"] == 1


@pytest.mark.parametrize("postfix", [None, "_writing"])
@pytest.mark.parametrize("ranked", [False, True])
def test_incomplete_cache_recomputed_even_with_validation_disabled(tmp_path, make_test_raster, postfix, ranked):
    config = recipe(tmp_path, make_test_raster)
    if ranked:
        config["footprints"]["param:image_rank_function"] = "{stat:area_crs}"
    Workflow(config).run()
    output = tmp_path / "footprints.gpkg"
    marker = Path(str(output) + ".incomplete")
    marker.write_text("Interrupted streaming output")
    config["shared"]["core:validity_writing_postfix"] = postfix
    workflow = Workflow(config)
    assert workflow.counts()["footprints"]["processing"] == 1
    assert len(workflow.run()) == (1 if ranked else 2)
    assert not marker.exists()


def test_invalid_geopackage_recomputed_when_validation_enabled(tmp_path, make_test_raster):
    config = recipe(tmp_path, make_test_raster, validity="before_reuse")
    (tmp_path / "footprints.gpkg").write_bytes(b"not a geopackage")
    workflow = Workflow(config)
    assert workflow.counts()["footprints"]["processing"] == 1
    assert len(workflow.run()) == 2


def test_reuse_hook_only_called_once_per_node_in_planning(tmp_path, make_test_raster, monkeypatch):
    config = recipe(tmp_path, make_test_raster)
    Workflow(config).run()
    calls = []
    original = CreateFootprints.can_reuse_outputs

    def check(self, *, params):
        calls.append(params["output_polygons"])
        return original(self, params=params)

    monkeypatch.setattr(CreateFootprints, "can_reuse_outputs", check)
    Workflow(config).counts()
    assert calls == [str(tmp_path / "footprints.gpkg")]
