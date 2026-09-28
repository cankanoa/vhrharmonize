"""Individual SpectralMatch adapters delegate algorithms and expose file products."""

import importlib
import inspect
import json
from pathlib import Path
from unittest.mock import Mock
import pytest
from osgeo import gdal
from vhrharmonize.parameters import function_parameters
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.registry import load_plugin, plugin_names
from workflow_helpers import import_settings, install_function, stage, transfer

FUNCTIONS = (
    "joint_coregistration",
    "global_regression",
    "local_block_adjustment",
    "align_rasters",
    "create_footprints",
    "postprocess_footprints",
    "voronoi_center_seamline",
    "weighted_seamline",
    "markov_triangles",
    "mask_rasters",
    "merge_rasters",
    "merge_vectors",
    "band_math",
    "create_cloud_mask_with_omnicloudmask",
    "process_raster_values_to_vector_polygons",
    "compute_overviews",
    "search_paths",
    "create_paths",
    "match_paths",
)


@pytest.mark.parametrize("name", FUNCTIONS)
def test_every_public_function_has_a_delegating_plugin(name, monkeypatch):
    module = importlib.import_module("vhrharmonize.plugins." + name)
    native = module._upstream()
    function = getattr(module, name)
    params = function_parameters(function)
    assert set(inspect.signature(native).parameters) == set(params)
    supplied = {
        key: "example" for key, value in params.items() if value.default is inspect.Parameter.empty
    }
    mock = Mock(return_value=["result.tif"])
    monkeypatch.setattr(module, "_upstream", lambda: mock)
    result = load_plugin(name).run(params=supplied, shared={})
    assert result == ["result.tif"]
    mock.assert_called_once_with(**supplied)


def test_old_pipeline_and_statistics_are_not_registered():
    names = plugin_names()
    assert "spectralmatch" not in names
    assert not any(name.startswith("compare_") for name in names)
    with pytest.raises(ValueError, match="Unknown workflow plugin"):
        load_plugin("spectralmatch")


def test_json_sequences_use_native_tuples_without_implicit_parameter_mapping(monkeypatch):
    module = importlib.import_module("vhrharmonize.plugins.global_regression")
    call = Mock(return_value=[])
    monkeypatch.setattr(module, "_upstream", lambda: call)
    load_plugin("global_regression").run(
        params={
            "input_images": ["a.tif"],
            "output_images": ["out.tif"],
            "window_scales": [2, 4],
            "vector_mask": ["include", "mask.gpkg", "image"],
            "custom_nodata_value": -9999,
            "output_dtype": "int16",
            "image_threads": 1,
            "dask_scheduler": None,
        },
        shared={"concurrent_processing": 8, "dask_scheduler_file": "not-implicitly-forwarded.json"},
    )
    args = call.call_args.kwargs
    assert args["window_scales"] == (2, 4) and args["vector_mask"] == (
        "include",
        "mask.gpkg",
        "image",
    )
    assert args["custom_nodata_value"] == -9999 and args["output_dtype"] == "int16"
    assert args["image_threads"] == 1 and args["dask_scheduler"] is None
    assert "resume_from_outputs" not in args and "shared_input_images" not in args


def test_standalone_merge_function_writes_real_raster(tmp_path, make_test_raster):
    from vhrharmonize.plugins.merge_rasters import merge_rasters

    source = make_test_raster(tmp_path / "input.tif", width=32, height=32)
    output = tmp_path / "mosaic.tif"
    assert merge_rasters([str(source)], str(output), io_threads=1, tile_threads=1) == str(output)
    with gdal.Open(str(output)) as dataset:
        assert dataset.RasterXSize == dataset.RasterYSize == 32
        assert (dataset.ReadAsArray() == 1).all()


def test_worldview_selects_one_chain_and_inherits_function_settings():
    from vhrharmonize.workflow.config import load_config

    config = load_config(Path(__file__).parents[1] / "configs/example.worldview.yml")
    setup = config["setup_spectralmatch"]
    assert "plugin" not in setup and setup["core:scope"] == "aggregate"
    assert list(config)[-1] == "merge_rasters"
    selected = [name for name in config if name in FUNCTIONS]
    assert selected == [
        "global_regression",
        "local_block_adjustment",
        "markov_triangles",
        "mask_rasters",
        "merge_rasters",
    ]
    assert all(config[name]["core:run"] for name in selected)
    assert (
        list(config).index("local_block_adjustment")
        < list(config).index("seamline_metadata")
        < list(config).index("markov_triangles")
    )
    assert config["seamline_metadata"]["param:image_paths"] == "collect:current_image_paths"
    assert config["shared"]["param:window_scales"] == [2, 4, 8, 16, 32, 64, 128, 256, 512]
    assert not any(
        key.startswith("const:")
        for key in config["shared"]
        if any(p in key for p in ("output_dtype", "custom_nodata_value", "window_scales"))
    )
    for name in selected:
        assert config[name]["plugin"] == name
        for parameter in (
            "output_dtype",
            "custom_nodata_value",
            "window_scales",
            "build_overviews",
            "window_size",
        ):
            assert "param:" + parameter not in config[name]
            assert "const:spectralmatch." + parameter not in setup
        if name != "markov_triangles":
            assert config[name]["core:calculate_overviews"] is True
    assert config["merge_rasters"]["param:output_image_path"].startswith("expr:")


def test_batch_outputs_checkpoint_reuse_cleanup_and_hpc(tmp_path, monkeypatch):
    source = tmp_path / "input.txt"
    source.write_text("input")
    outputs = [str(tmp_path / "temp/a.txt"), str(tmp_path / "temp/b.txt")]
    final = str(tmp_path / "output/done.txt")
    calls = []

    def batch(input_images, output_images):
        calls.append("batch")
        for filename in output_images:
            Path(filename).write_text("batch")
        return {"images": output_images, "gain": 2}

    def finish(input_images, output_path, gain):
        calls.append("finish")
        assert gain == 2
        Path(output_path).write_text("|".join(Path(p).read_text() for p in input_images))

    install_function(
        monkeypatch,
        "batch",
        batch,
        scope="aggregate",
        input_paths={"input_images"},
        output_paths={"output_images"},
    )
    install_function(
        monkeypatch,
        "finish",
        finish,
        scope="aggregate",
        input_paths={"input_images"},
        output_paths={"output_path"},
    )
    config = {
        "defaults": {
            "plugin": "shared", "core:run": True, "core:log_to_console": False,
            # Retain intermediates for the checkpoint/staging checks; enable cleanup below.
            "core:delete_temp_steps_proactively": False,
        },
        "files": import_settings(source, tmp_path),
        "matching": {
            "plugin": "batch",
            "core:run": True,
            "param:input_images": [str(source)],
            "param:output_images": outputs,
            "const:result": "returned:$",
        },
        "mosaic": {
            "plugin": "finish",
            "core:run": True,
            "param:input_images": outputs,
            "param:output_path": final,
            "param:gain": "const:result.gain",
        },
    }
    Workflow(config).run()
    calls.clear()
    Workflow(config).run()
    assert not calls
    Path(final).unlink()
    staged, uploads, downloads = stage(config, tmp_path)
    assert set(outputs) <= uploads.keys() and str(source) not in uploads
    assert outputs[0] + ".context.json" in uploads
    transfer(uploads)
    remote = Workflow(staged)
    remote.run()
    assert calls == ["finish"]
    assert remote.context["const"]["result"]["images"] == [uploads[p] for p in outputs]
    assert Path(downloads[final]).read_text() == "batch|batch"
    config["defaults"]["core:delete_temp_steps_proactively"] = True
    Workflow(config).run()
    assert not any(Path(p).exists() for p in outputs)


def test_partial_batch_is_rebuilt_and_output_collisions_are_detected(tmp_path, monkeypatch):
    source = tmp_path / "input.txt"
    source.write_text("source")
    outputs = [str(tmp_path / "output/a.txt"), str(tmp_path / "output/b.txt")]
    calls = []

    def batch(output_images):
        calls.append(1)
        for name in output_images:
            Path(name).write_text("complete")

    install_function(monkeypatch, "batch", batch, scope="aggregate", output_paths={"output_images"})
    recipe = {
        "files": import_settings(source, tmp_path),
        "batch": {"plugin": "batch", "core:run": True, "param:output_images": outputs},
    }
    Path(outputs[0]).parent.mkdir()
    Path(outputs[0]).write_text("partial")
    Workflow(recipe).run()
    assert calls == [1]
    recipe["duplicate"] = {"plugin": "batch", "core:run": True, "param:output_images": [outputs[1]]}
    with pytest.raises(ValueError, match="collision"):
        Workflow(recipe)


def test_worldview_function_chain_runs_stages_and_reuses_final_product(
    tmp_path, monkeypatch, make_test_raster
):
    """Exercise the editable recipe's actual links, ordering and native merge options."""
    from functools import wraps
    from vhrharmonize.workflow.config import load_config

    source = make_test_raster(tmp_path / "scene.tif", width=32, height=32)
    example = load_config(Path(__file__).parents[1] / "configs/example.worldview.yml")
    recipe = {
        "defaults": {
            **example["shared"],
            "core:log_to_console": False,
            "core:validity_check_grid_size": 0,
            "param:custom_nodata_value": -9999,
            "param:output_dtype": "float32",
            "param:epsg": 4326,
            "param:window_scales": [2, 4],
        },
        "files": {
            **import_settings(source, tmp_path),
            "var:current_image_paths": "returned:file_path",
            "var:footprint_metadata": {"quality": 10},
        },
        "setup": {**example["setup_spectralmatch"], "core:run": True},
        **{
            name: {**example[name], "core:run": True}
            for name in (
                "global_regression",
                "local_block_adjustment",
                "seamline_metadata",
                "markov_triangles",
                "mask_rasters",
                "merge_rasters",
            )
        },
    }
    calls = []
    for name in ("global_regression", "local_block_adjustment"):
        module = importlib.import_module("vhrharmonize.plugins." + name)
        native = module._upstream()

        def fake(*, _name=name, **kwargs):
            calls.append((_name, kwargs["input_images"], kwargs["output_images"]))
            assert kwargs["custom_nodata_value"] == -9999
            assert kwargs["output_dtype"] == "float32"
            assert "build_overviews" not in kwargs
            for source, destination in zip(kwargs["input_images"], kwargs["output_images"]):
                assert Path(source).is_file()
                gdal.Translate(
                    destination,
                    source,
                    outputType=gdal.GDT_Float32,
                    noData=kwargs["custom_nodata_value"],
                )
            return kwargs["output_images"]

        fake = wraps(native)(fake)
        monkeypatch.setattr(module, "_upstream", lambda f=fake: f)

    # The solver is outside this recipe test; check its actual footprint/image binding.
    import geopandas as gpd

    module = importlib.import_module("vhrharmonize.plugins.markov_triangles")
    native = module._upstream()

    @wraps(native)
    def markov(**kwargs):
        footprints = gpd.read_file(kwargs["input_polygons"], layer=kwargs["input_layer"])
        assert list(footprints[kwargs["image_field_name"]]) == kwargs["input_images"]
        footprints.to_file(kwargs["output_mask"], layer="seamlines", driver="GPKG")
        return kwargs["output_mask"]

    monkeypatch.setattr(module, "_upstream", lambda: markov)

    # Verify the final invocation separately from the native mosaic algorithm
    # (covered by test_standalone_merge_function_writes_real_raster).
    module = importlib.import_module("vhrharmonize.plugins.merge_rasters")
    native = module._upstream()

    @wraps(native)
    def merge(**kwargs):
        import rasterio

        assert kwargs["image_threads"] is None
        assert kwargs["concurrent_processing_backend"] is None
        assert kwargs["dask_scheduler"] is None
        assert kwargs["output_dtype"] == "float32"
        assert "build_overviews" not in kwargs
        assert kwargs["input_images"][0].endswith("_global_local_masked.tif")
        with rasterio.open(kwargs["input_images"][0]) as source:
            with rasterio.open(kwargs["output_image_path"], "w", **source.profile) as target:
                target.write(source.read())
        return kwargs["output_image_path"]

    monkeypatch.setattr(module, "_upstream", lambda: merge)

    staged, uploads, downloads = stage(recipe, tmp_path)
    transfer(uploads)
    remote = Workflow(staged)
    remote.run()
    final = str(tmp_path / "output/mosaic.tif")
    assert Path(downloads[final]).is_file()
    assert [call[0] for call in calls] == ["global_regression", "local_block_adjustment"]
    assert calls[0][2] == calls[1][1]
    assert calls[1][2][0].endswith("_global_local.tif")
    assert remote.records[0]["context"]["var"]["current_image_paths"] == downloads[final]
    with gdal.Open(downloads[final]) as mosaic:
        assert mosaic.GetRasterBand(1).GetOverviewCount() == 2
    calls.clear()
    # Simulate downloading only the final product; missing intermediates must not force recomputation.
    Path(final).parent.mkdir(parents=True, exist_ok=True)
    import shutil

    shutil.copy2(downloads[final], final)
    local = Workflow(recipe)
    local.run()
    assert not calls
    assert local.counts()["global_regression"]["processing"] == 0
    assert local.counts()["local_block_adjustment"]["processing"] == 0


def test_tiled_output_folder_never_counts_as_complete_cache(tmp_path, monkeypatch):
    directory = tmp_path / "output/tiles"
    directory.mkdir(parents=True)
    (directory / "one.txt").write_text("incomplete")
    calls = []

    def batch(output_path):
        calls.append(output_path)
        (Path(output_path) / "two.txt").write_text("remaining")

    install_function(monkeypatch, "tiles", batch, output_paths={"output_path"})
    source = tmp_path / "source.txt"
    source.write_text("source")
    recipe = {
        "files": import_settings(source, tmp_path),
        "tile_merge": {"plugin": "tiles", "core:run": True, "param:output_path": str(directory)},
    }
    Workflow(recipe).run()
    assert calls == [str(directory)]
    assert (directory / "two.txt").is_file()
