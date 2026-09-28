"""The editable recipe documents every public function parameter and its defaults."""

import importlib
import inspect
from pathlib import Path
import re

import pytest
import yaml

from vhrharmonize.parameters import function_parameters
from test_spectralmatch_workflow import FUNCTIONS

RECIPE = Path(__file__).resolve().parents[1] / "configs/example.worldview.yml"


def plugin_sections():
    result = {}
    text = re.sub(r"^# ", "", RECIPE.read_text(), flags=re.MULTILINE)
    for section in re.split(r"(?=^\w+:$)", text, flags=re.MULTILINE):
        match = re.search(r"^  plugin: (\w+)$", section, re.MULTILINE)
        if match:
            result.setdefault(match[1], []).append(section)
    return result


def documented_plugin_keys():
    return {
        name: set(re.findall(r"^  (?:#\s*)?param:(\w+):", "\n".join(sections), re.MULTILINE))
        for name, sections in plugin_sections().items()
    }


@pytest.mark.parametrize(
    "plugin,functions",
    [
        ("import_files", ["import_files"]),
        ("file_source", ["copy_file"]),
        ("fetch_dem", ["download_opentopography_dem_for_bbox"]),
        (
            "fetch_atmosphere",
            ["fetch_power_atmosphere_for_bbox", "fetch_modis_water_vapor_for_bbox"],
        ),
        ("atmospheric_correction", ["run_py6s", "run_flaash"]),
        ("orthorectification", ["gcp_refined_rpc_orthorectification"]),
        ("pansharpen", ["pansharpen_image"]),
        ("cloud_mask", ["cloudmask_raster"]),
        ("alignment", ["align_image_pair"]),
        ("seamline_metadata", ["write_seamline_metadata_gpkg"]),
        *[(name, [name]) for name in FUNCTIONS],
    ],
)
def test_worldview_documents_every_public_plugin_parameter(plugin, functions):
    module = importlib.import_module("vhrharmonize.plugins." + plugin)
    documented = documented_plugin_keys()[plugin]
    for name in functions:
        assert set(function_parameters(getattr(module, name))) <= documented


def test_worldview_retains_all_supported_flaash_task_parameters():
    from vhrharmonize.plugins.atmospheric_correction import FLAASH_ALLOWED_PARAMS

    text = "\n".join(plugin_sections()["atmospheric_correction"])
    assert FLAASH_ALLOWED_PARAMS <= set(
        re.findall(r"^\s+#\s+([A-Z][A-Z_0-9]+):", text, re.MULTILINE)
    )


@pytest.mark.parametrize("name", FUNCTIONS)
def test_spectralmatch_commented_defaults_match_the_installed_api(name):
    text = "\n".join(plugin_sections()[name])
    function = getattr(importlib.import_module("vhrharmonize.plugins." + name), name)
    for key, parameter in function_parameters(function).items():
        line = re.search(r"^  (#\s*)?param:" + key + r": (.*?) # (.*)$", text, re.MULTILINE)
        assert line, f"{name}.{key} needs a value and a concise options/default comment"
        if parameter.default is inspect.Parameter.empty:
            assert "required" in line[3]
            continue
        assert "default:" in line[3], (name, key)
        if line[1]:
            default = parameter.default
            if isinstance(default, tuple):
                default = list(default)
            assert yaml.safe_load(line[2]) == default, (name, key)


def test_every_example_step_is_explicitly_enabled_or_disabled():
    for filename in (RECIPE, RECIPE.with_name("example.planet.yml")):
        recipe = yaml.safe_load(filename.read_text())
        assert all(isinstance(step.get("core:run"), bool) for step in recipe.values())
