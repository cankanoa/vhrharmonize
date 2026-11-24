"""The editable recipe uses supported options and accurate commented defaults."""

import importlib
import inspect
from pathlib import Path
import re

import pytest
import yaml

from vhrharmonize.parameters import function_parameters
from vhrharmonize.workflow.registry import load_plugin
from test_spectralmatch_workflow import FUNCTIONS

RECIPE = Path(__file__).resolve().parents[1] / "configs/example.worldview.yml"


def plugin_sections():
    result = {}
    text = re.sub(r"^# ", "", RECIPE.read_text(), flags=re.MULTILINE)
    for section in re.split(r"(?=^\w+:$)", text, flags=re.MULTILINE):
        match = re.search(r"^  plugin: (\w+)(?:\s+#.*)?\s*$", section, re.MULTILINE)
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
        ("py6s_atmospheric_correction", ["run_py6s"]),
        ("flaash_atmospheric_correction", ["run_flaash"]),
        ("orthorectification", ["gcp_refined_rpc_orthorectification"]),
        ("pansharpen", ["pansharpen_image"]),
        ("cloud_mask", ["cloudmask_raster"]),
        ("alignment", ["align_image_pair"]),
        *[(name, [name]) for name in FUNCTIONS],
    ],
)
def test_worldview_parameters_are_supported_by_their_plugin(plugin, functions):
    module = importlib.import_module("vhrharmonize.plugins." + plugin)
    documented = documented_plugin_keys().get(plugin, set())
    accepted = set(load_plugin(plugin).options)
    if plugin == 'fetch_atmosphere':
        accepted.update({'source', 'output_path'})  # Consumed by the adapter.
    for name in functions:
        accepted.update(function_parameters(getattr(module, name)))
    # This is an executable workflow, not an exhaustive catalog of every plugin.
    # Each function's full public API is covered by its adapter contract tests.
    assert documented <= accepted


def test_worldview_retains_all_supported_flaash_task_parameters():
    from vhrharmonize.plugins.flaash_atmospheric_correction import FLAASH_ALLOWED_PARAMS

    text = "\n".join(plugin_sections()["flaash_atmospheric_correction"])
    assert FLAASH_ALLOWED_PARAMS <= set(
        re.findall(r"^\s+(?:#\s+)?([A-Z][A-Z_0-9]+):", text, re.MULTILINE)
    )


@pytest.mark.parametrize("name", FUNCTIONS)
def test_spectralmatch_commented_defaults_match_the_installed_api(name):
    text = "\n".join(plugin_sections().get(name, []))
    function = getattr(importlib.import_module("vhrharmonize.plugins." + name), name)
    parameters = function_parameters(function)
    for key, value, comment in re.findall(r"^  #\s*param:(\w+): (.*?) # (.*)$", text, re.MULTILINE):
        parameter = parameters[key]
        assert parameter.default is not inspect.Parameter.empty, (name, key)
        # Expressions show optional wiring, and inherited settings document the
        # shared override. Literal commented defaults must match the native API.
        if value.startswith(('var:', 'const:', 'expr:')) or 'default: shared.' in comment.lower():
            continue
        default = parameter.default
        if isinstance(default, tuple):
            default = list(default)
        assert yaml.safe_load(value) == default, (name, key)


def test_every_example_step_is_explicitly_enabled_or_disabled():
    for filename in (RECIPE, RECIPE.with_name("example.planet.yml")):
        recipe = yaml.safe_load(filename.read_text())
        assert all(isinstance(step.get("core:run"), bool) for step in recipe.values())
