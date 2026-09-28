"""Processing functions use explicit destinations; naming is ordinary YAML state."""

from importlib import import_module

import pytest

from vhrharmonize.parameters import function_parameters
from vhrharmonize.plugins.file_source import copy_file
from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import import_settings, copy_step


@pytest.mark.parametrize(
    "module,function,arguments,output",
    [
        ("file_source", "copy_file", {"input_path": "input.tif"}, "output_path"),
        (
            "alignment",
            "align_image_pair",
            {"moving_image_path": "input.tif", "fixed_image_path": "ref.tif"},
            "output_image_path",
        ),
        ("atmospheric_correction", "run_py6s", {"input_raster": "input.tif"}, "output_raster"),
        ("atmospheric_correction", "run_flaash", {"input_raster": "input.tif"}, "output_raster"),
        (
            "orthorectification",
            "gcp_refined_rpc_orthorectification",
            {
                "input_image_path": "input.tif",
                "dem_image_path": "dem.tif",
                "output_epsg": 4326,
            },
            "output_image_path",
        ),
        (
            "pansharpen",
            "pansharpen_image",
            {
                "input_low_resolution_path": "mul.tif",
                "input_high_resolution_path": "pan.tif",
            },
            "output_image_path",
        ),
        (
            "cloud_mask",
            "cloudmask_raster",
            {
                "input_image_path": "input.tif",
                "output_mask_path": "mask.tif",
                "red_band_index": 1,
                "green_band_index": 2,
                "nir_band_index": 3,
            },
            "output_raster_path",
        ),
    ],
)
def test_missing_output_is_rejected_before_processing(module, function, arguments, output):
    target = getattr(import_module("vhrharmonize.plugins." + module), function)
    assert "append_to_name" not in function_parameters(target)
    for value in (None, ""):
        with pytest.raises(ValueError, match=output + " must be an explicit output path"):
            target(**arguments, **{output: value})


def test_removed_naming_option_is_rejected_in_python_and_yaml(tmp_path):
    source = tmp_path / "input.txt"
    source.write_text("input")
    destination = tmp_path / "output.txt"
    with pytest.raises(TypeError, match="append_to_name"):
        copy_file(source, destination, append_to_name="_old")
    recipe = {
        "import": import_settings(source, tmp_path),
        "copy": {**copy_step("copied", "mul", str(destination)), "param:append_to_name": "_old"},
    }
    with pytest.raises(ValueError, match="Unsupported FileSource options: append_to_name"):
        Workflow(recipe).run()
    assert not destination.exists()
