from copy import deepcopy
import json
from pathlib import Path
import shutil
import pytest
import yaml
from vhrharmonize.workflow.config import load_config, validate_config
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.staging import stage_workflow
from vhrharmonize.workflow.values import expression, resolve
from workflow_helpers import import_settings, copy_step, install_function


@pytest.fixture
def pipeline(tmp_path, make_test_raster):
    raw = make_test_raster(tmp_path / "source/image.tif")
    config = {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False},
        "import_files": import_settings(raw, tmp_path),
        'file_source_1': {**(copy_step("corrected", "mul", "expr:const.temp_dir & '/corrected.tif'")), "plugin": 'file_source'}, 'file_source_2': {**(copy_step("cloudmasked", "corrected", "expr:const.output_dir & '/cloudmasked.tif'")), "plugin": 'file_source'}, 'file_source_3': {**(copy_step("aligned", "cloudmasked", "expr:const.output_dir & '/aligned.tif'")), "plugin": 'file_source'},
    }
    return (config, raw)


def test_resume_skips_missing_upstream_files_and_processes_only_alignment(
    pipeline, make_test_raster
):
    config, raw = pipeline
    masked = make_test_raster(
        Path(config["import_files"]["const:output_dir"]) / "cloudmasked.tif", fill=7
    )
    before = masked.read_bytes()
    workflow = Workflow(config)
    counts = workflow.counts()
    assert counts["file_source_1"] == {"loaded": 0, "processing": 0, "unused": 1}
    assert counts["file_source_2"]["loaded"] == 1
    assert counts["file_source_3"]["processing"] == 1
    assert not Path(config["import_files"]["const:temp_dir"]).exists()
    workflow.run()
    assert Path(config["import_files"]["const:output_dir"], "aligned.tif").read_bytes() == before
    assert masked.read_bytes() == before
    assert not Path(config["import_files"]["const:temp_dir"], "corrected.tif").exists()


def test_saved_intermediate_remains_a_required_output(pipeline, make_test_raster):
    config, raw = pipeline
    config["file_source_1"]["var:corrected"] = "expr:const.output_dir & '/corrected.tif'"
    make_test_raster(Path(config["import_files"]["const:output_dir"]) / "cloudmasked.tif")
    counts = Workflow(config).counts()
    assert counts["file_source_1"]["processing"] == 1
    assert counts["file_source_2"]["processing"] == 0


def test_each_record_resumes_at_its_own_available_output(pipeline, make_test_raster):
    config, raw = pipeline
    make_test_raster(raw.with_name("second.tif"))
    config["import_files"]["param:search_glob"] = str(raw.parent / "*.tif")
    for step, name in zip([config[f"file_source_{i}"] for i in (1, 2, 3)], ["corrected", "cloudmasked", "aligned"]):
        step["var:" + name] = "expr:const.temp_dir & '/' & var.basename & '_" + name + ".tif'"
    root = Path(config["import_files"]["const:temp_dir"])
    make_test_raster(root / "image_cloudmasked.tif")
    make_test_raster(root / "second_corrected.tif")
    workflow = Workflow(config)
    assert [r["processing"] for r in workflow.counts().values()] == [0, 1, 2]
    workflow.run()
    assert len(workflow.records) == 2


def test_corrupt_output_is_regenerated_and_sources_preserved(pipeline):
    config, raw = pipeline
    original = raw.read_bytes()
    corrupt = Path(config["import_files"]["const:output_dir"], "cloudmasked.tif")
    corrupt.parent.mkdir()
    corrupt.write_bytes(b"II*\x00\x08\x02\x00\x00")
    workflow = Workflow(config)
    assert workflow.counts()["file_source_1"]["processing"] == 1
    assert workflow.counts()["file_source_2"]["processing"] == 1
    assert corrupt.read_bytes() != original
    workflow.run()
    assert corrupt.read_bytes() == original
    assert raw.read_bytes() == original


def test_reuse_can_be_disabled(pipeline, make_test_raster):
    config, raw = pipeline
    make_test_raster(Path(config["import_files"]["const:output_dir"]) / "cloudmasked.tif")
    config["shared"]["core:run_from_existing"] = False
    assert all((row["processing"] == 1 for row in Workflow(config).counts().values()))


def test_yaml_order_is_authoritative_and_forward_references_fail(pipeline):
    config, raw = pipeline
    config["file_source_1"], config["file_source_2"] = (
        config["file_source_2"],
        config["file_source_1"],
    )
    with pytest.raises(ValueError, match="Undefined variable field: corrected"):
        Workflow(config)


def test_disabled_step_does_not_publish_an_input_alias(pipeline):
    config, raw = pipeline
    config["file_source_1"]["core:run"] = False
    config["file_source_1"]["var:corrected"] = "var:mul"
    with pytest.raises(ValueError, match="Undefined variable field: corrected"):
        Workflow(config)


def test_explicit_input_link_skips_a_disabled_step(pipeline):
    config, raw = pipeline
    config["file_source_1"]["core:run"] = False
    config["file_source_2"]["param:input_path"] = "var:mul"
    workflow = Workflow(config)
    workflow.run()
    assert "corrected" not in workflow.contexts[0]["var"]
    assert (
        Path(config["import_files"]["const:output_dir"], "aligned.tif").read_bytes()
        == raw.read_bytes()
    )
    assert not Path(config["import_files"]["const:temp_dir"], "corrected.tif").exists()


def test_failed_consumer_preserves_intermediate(pipeline, monkeypatch):
    from vhrharmonize.plugins.file_source import FileSource

    config, raw = pipeline
    config["shared"]["core:delete_temp_steps_proactively"] = True
    run = FileSource.run

    def fail(self, **kwargs):
        if kwargs["params"]["output_path"].endswith("cloudmasked.tif"):
            raise RuntimeError("failed")
        return run(self, **kwargs)

    monkeypatch.setattr(FileSource, "run", fail)
    with pytest.raises(RuntimeError, match="failed"):
        Workflow(config).run()
    assert Path(config["import_files"]["const:temp_dir"], "corrected.tif").exists()


def test_metadata_mapping_and_json_handoff(pipeline, monkeypatch):
    from vhrharmonize.plugins.file_source import FileSource

    config, raw = pipeline
    raw.with_suffix(".IMD").write_text("BEGIN_GROUP = IMAGE\nangle = 20;\nEND_GROUP = IMAGE\nEND;")
    config["import_files"].update(
        {
            "param:create_metadata_json": {"source": {"to_json": "literal:expr:$replace(var.file_path, '.tif', '.IMD')"}},
            "var:solar_zenith": "expr:90 - var.source.IMAGE.angle",
            "var:custom": "sensor-independent",
        }
    )
    config["file_source_1"]["var:processed"] = "expr:var.solar_zenith + 1"
    config["file_source_2"]["var:copied"] = "var:processed"
    for step in [config[f"file_source_{i}"] for i in (1, 2, 3)]:
        step["param:snapshot"] = "var:$"
    observed = []
    run = FileSource.run

    def observe(self, **kwargs):
        observed.append(deepcopy(kwargs["params"].pop("snapshot")))
        return run(self, **kwargs)

    monkeypatch.setattr(FileSource, "run", observe)
    workflow = Workflow(config)
    workflow.run()
    assert observed[0]["solar_zenith"] == 70
    assert observed[1]["processed"] == 71
    assert workflow.records[0]["context"]["var"]["copied"] == 71
    Path(config["import_files"]["const:output_dir"], "aligned.tif").unlink()
    observed.clear()
    Workflow(config).run()
    assert len(observed) == 1 and observed[0]["copied"] == 71


def test_function_plugins_use_shared_and_explicit_arguments(monkeypatch):
    from vhrharmonize.plugins.base import FunctionPlugin

    calls = []

    def function(input_path, output_path, gain, custom_nodata_value, variables):
        calls.append((input_path, output_path, gain, custom_nodata_value, variables))

    plugin = FunctionPlugin()
    monkeypatch.setattr(plugin, "function", lambda: function)
    plugin.run(
        params={"input_path": "in", "output_path": "out", "gain": 3, "variables": {"gain": 1}},
        shared={"gain": 2, "custom_nodata_value": -9999, "concurrent_processing": 2},
    )
    assert calls == [("in", "out", 3, -9999, {"gain": 1})]
    with pytest.raises(ValueError, match="Unsupported"):
        plugin.run(params={"typo": 1}, shared={})


def test_staging_resumes_from_cloudmasked_output(pipeline, make_test_raster, tmp_path):
    config, raw = pipeline
    masked = make_test_raster(Path(config["import_files"]["const:output_dir"]) / "cloudmasked.tif")
    remote_root = tmp_path / "remote"
    staged, uploads, downloads = stage_workflow(
        config,
        config_dir=str(tmp_path),
        remote_output_dir=str(remote_root / "out"),
        remote_temp_dir=str(remote_root / "temp"),
        remote_reference_dir=str(remote_root / "ref"),
    )
    assert str(raw) not in uploads
    assert str(masked) in uploads
    assert str(Path(config["import_files"]["const:output_dir"], "aligned.tif")) in downloads
    for source, target in uploads.items():
        Path(target).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    remote = Workflow(staged)
    assert remote.counts()["file_source_1"]["processing"] == 0
    assert remote.counts()["file_source_3"]["processing"] == 1
    remote.run()
    assert (remote_root / "out" / "aligned.tif").exists()


def test_path_collisions_are_rejected_before_processing(pipeline):
    config, raw = pipeline
    config["file_source_2"]["var:cloudmasked"] = "expr:const.temp_dir & '/corrected.tif'"
    with pytest.raises(ValueError, match="collision"):
        Workflow(config)


def test_old_flat_workflow_and_duplicate_keys_rejected(tmp_path):
    for bad in [{"workflow": []}, {"typo": {"plugin": "unknown"}}, {"output_dir": "./output"}]:
        with pytest.raises(ValueError):
            validate_config(bad)
    filename = tmp_path / "bad.yml"
    filename.write_text("import_files:\n  core:run: true\nimport_files: {}\n")
    with pytest.raises(ValueError, match="Duplicate YAML key"):
        load_config(filename)


def test_expressions_preserve_types_and_disallow_python_execution():
    assert resolve("var:values", {"const": {}, "var": {"values": [1, 2]}}) == [1, 2]
    assert expression(
        "$map(var.source, function($v){$v * 2})", {"const": {}, "var": {"source": [1, 2]}}
    ) == [2, 4]
    with pytest.raises(ValueError):
        expression("__import__('os').system('true')", {"const": {}, "var": {}})


def test_plugin_entry_point_registration(monkeypatch):
    from vhrharmonize.workflow import registry
    from vhrharmonize.plugins.file_source import FileSource

    class Entries(list):

        def select(self, **kwargs):
            assert kwargs == {"group": "vhrharmonize.plugins"}
            return self

    class Entry:
        name = "custom"

        def load(self):
            return FileSource

    monkeypatch.setattr(registry.metadata, "entry_points", lambda: Entries([Entry()]))
    assert isinstance(registry.load_plugin("custom"), FileSource)
    with pytest.raises(ValueError, match="Unknown workflow plugin"):
        registry.load_plugin("absent")


def test_actual_alignment_adapter_receives_cached_cloudmasked_image(
    pipeline, make_test_raster, monkeypatch
):
    from types import SimpleNamespace
    from vhrharmonize.plugins import alignment
    from vhrharmonize.plugins.cloud_mask import CloudMask

    config, raw = pipeline
    root = Path(config["import_files"]["const:output_dir"])
    masked = make_test_raster(root / "cloudmasked.tif", fill=9)
    mask = make_test_raster(root / "mask.tif")
    reference = make_test_raster(raw.parent / "reference.tif")
    config["file_source"] = config.pop("file_source_1")
    config.pop("file_source_2")
    config.pop("file_source_3")
    config["cloud_mask"] = {"plugin": 'cloud_mask', 
        "core:run": True,
        "param:input_image_path": "var:corrected",
        "var:cloudmasked": str(masked),
        "param:output_raster_path": "var:cloudmasked",
        "param:output_mask_path": str(mask),
    }
    config["alignment"] = {"plugin": 'alignment', 
        "core:run": True,
        "param:moving_image_path": "var:cloudmasked",
        "param:fixed_image_path": str(reference),
        "param:output_image_path": str(root / "aligned.tif"),
        "param:moving_band_index": 0,
        "param:fixed_band_index": 0,
    }
    calls = []

    def align(**kwargs):
        calls.append(kwargs)
        shutil.copy2(kwargs["moving_image_path"], kwargs["output_image_path"])
        return SimpleNamespace(output_image_path=kwargs["output_image_path"])

    monkeypatch.setattr(alignment, "coregix_align_image_pair", align)
    monkeypatch.setattr(
        CloudMask, "run", lambda *a, **k: pytest.fail("Cloud mask should be reused")
    )
    Workflow(config).run()
    assert calls[0]["moving_image_path"] == str(masked)
    assert (root / "aligned.tif").read_bytes() == masked.read_bytes()
    assert not Path(config["import_files"]["const:temp_dir"], "corrected.tif").exists()


def test_dask_executes_required_records_in_yaml_order(pipeline, make_test_raster, monkeypatch):
    from concurrent.futures import Future
    from types import ModuleType
    import sys

    config, raw = pipeline
    make_test_raster(raw.with_name("second.tif"))
    config["import_files"]["param:search_glob"] = str(raw.parent / "*.tif")
    for step, name in zip([config[f"file_source_{i}"] for i in (1, 2, 3)], ["corrected", "cloudmasked", "aligned"]):
        step["var:" + name] = "expr:const.output_dir & '/' & var.basename & '_" + name + ".tif'"
    config["shared"].update(
        {
            "core:concurrent_processing_backend": "dask",
            "core:dask_scheduler_address": "tcp://scheduler:8786",
        }
    )
    submitted = []

    class Client:

        def __init__(self, address):
            assert address == "tcp://scheduler:8786"

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def submit(self, function, payload, pure):
            assert pure is False
            submitted.append(Path(payload[1]["output_path"]).name)
            future = Future()
            future.set_result(function(payload))
            return future

        def cancel(self, futures):
            pass

    distributed = ModuleType("dask.distributed")
    distributed.Client = Client
    distributed.as_completed = lambda futures: reversed(list(futures))
    monkeypatch.setitem(sys.modules, "dask.distributed", distributed)
    Workflow(config).run()
    assert submitted == [
        "image_corrected.tif",
        "second_corrected.tif",
        "image_cloudmasked.tif",
        "second_cloudmasked.tif",
        "image_aligned.tif",
        "second_aligned.tif",
    ]


@pytest.mark.parametrize("alias", ["symlink", "hardlink"])
def test_outputs_never_overwrite_source_aliases(pipeline, alias):
    import os

    config, raw = pipeline
    output = Path(config["import_files"]["const:output_dir"], "aligned.tif")
    output.parent.mkdir()
    if alias == "symlink":
        output.symlink_to(raw)
    else:
        os.link(raw, output)
    original = raw.read_bytes()
    config["shared"]["core:run_from_existing"] = False
    with pytest.raises(ValueError, match="protected input"):
        Workflow(config).run()
    assert raw.read_bytes() == original


def test_aggregate_metadata_reaches_following_scene_steps(pipeline, monkeypatch):
    config, raw = pipeline
    install_function(monkeypatch, "summary", lambda: {"total": 42}, scope="aggregate")
    for i in (1, 2, 3):
        config.pop(f"file_source_{i}")
    config["summary"] = {"plugin": 'summary', "core:run": True, "const:total": "returned:total"}
    config["file_source"] = copy_step("copied", "mul", "expr:const.output_dir & '/copied.tif'")
    config["file_source"]["var:seen"] = "const:total"
    workflow = Workflow(config)
    workflow.run()
    assert workflow.records[0]["context"]["var"]["seen"] == 42


def test_cached_atmosphere_restores_metadata_when_checkpoint_is_corrupt(pipeline, monkeypatch):
    from vhrharmonize.plugins.fetch_atmosphere import FetchAtmosphere

    config, raw = pipeline
    atmosphere = raw.parent / "atmosphere.json"
    atmosphere.write_text('{"water_vapor":2.5}')
    Path(str(atmosphere) + ".context.json").write_text("{broken")
    for i in (1, 2, 3):
        config.pop(f"file_source_{i}")
    config["fetch_atmosphere"] = {"plugin": 'fetch_atmosphere', 
        "core:run": True,
        "param:output_path": str(atmosphere),
        "var:atmosphere": "returned:$",
    }
    config["file_source"] = copy_step("copied", "mul", "expr:const.output_dir & '/copied.tif'")
    config["file_source"]["var:water_vapor"] = "var:atmosphere.water_vapor"
    monkeypatch.setattr(FetchAtmosphere, "run", lambda *a, **k: pytest.fail("Must reuse JSON"))
    workflow = Workflow(config)
    workflow.run()
    assert workflow.records[0]["context"]["var"]["water_vapor"] == 2.5


def test_cached_raster_does_not_require_unused_temporary_mask(
    pipeline, make_test_raster, monkeypatch
):
    config, raw = pipeline
    masked = make_test_raster(Path(config["import_files"]["const:output_dir"]) / "cloudmasked.tif")
    config["file_source"] = config.pop("file_source_1")
    config.pop("file_source_2")
    config.pop("file_source_3")
    config["cloud_mask"] = {"plugin": 'cloud_mask', 
        "core:run": True,
        "param:input_image_path": "var:corrected",
        "var:cloudmasked": str(masked),
        "param:output_raster_path": "var:cloudmasked",
        "param:output_mask_path": "expr:const.temp_dir & '/mask.tif'",
    }
    install_function(
        monkeypatch,
        "finish",
        lambda input_path, output_path: shutil.copy2(input_path, output_path),
        input_paths={"input_path"},
        output_paths={"output_path"},
    )
    config["finish"] = copy_step("aligned", "cloudmasked", "expr:const.output_dir & '/aligned.tif'")
    config["finish"]["plugin"] = "finish"
    workflow = Workflow(config)
    workflow.plan()
    assert workflow.counts()["file_source"]["processing"] == 0
    assert workflow.counts()["cloud_mask"]["processing"] == 0
    workflow.run()
    assert masked.exists()
