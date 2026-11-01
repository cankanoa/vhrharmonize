"""Atmospheric correction algorithms have independent workflow plugins."""

import importlib
from pathlib import Path
import re

import pytest
import yaml

from vhrharmonize.parameters import function_parameters
from vhrharmonize.workflow.registry import BUILTINS, load_plugin, plugin_names


@pytest.mark.parametrize("algorithm,foreign_option", [("py6s", "envi_engine"), ("flaash", "aot550")])
def test_algorithm_plugins_expose_only_their_own_parameters(algorithm, foreign_option):
    name = f"{algorithm}_atmospheric_correction"
    plugin = load_plugin(name)
    assert name in plugin_names()
    assert plugin.function().__module__ == f"vhrharmonize.plugins.{name}"
    assert "method" not in function_parameters(plugin.function())
    for option in (foreign_option, "method"):
        with pytest.raises(ValueError, match="Unsupported .* options: " + option):
            plugin.arguments(plugin.function(), {option: "unused"}, {})
    features = plugin.file_features()
    assert ("wavelengths_response_csv" in features["input_hpc_staging_paths"]) == (algorithm == "py6s")
    assert ("output_params_path" in features["output_hpc_download_paths"]) == (algorithm == "flaash")
    assert "atmospheric_correction" not in BUILTINS


@pytest.mark.parametrize("algorithm", ["py6s", "flaash"])
def test_plugins_run_the_selected_algorithm_and_preserve_public_python_exports(monkeypatch, tmp_path, algorithm):
    import vhrharmonize

    name = f"{algorithm}_atmospheric_correction"
    module = importlib.import_module("vhrharmonize.plugins." + name)
    assert getattr(vhrharmonize, "run_" + algorithm) is getattr(module, "run_" + algorithm)
    calls = []
    params = {"input_raster": "input.tif", "output_raster": str(tmp_path / "out.tif")}
    if algorithm == "py6s":
        monkeypatch.setattr(module.Py6SCorrector, "run", lambda self, **kw: calls.append(kw))
        params.update(wavelengths_center_um=[0.55], ground_elevation_km=0)
    else:
        monkeypatch.setattr(module, "_execute_flaash_task", lambda *args, **kw: calls.append(args))
        params.update(params={}, envi_engine=object())
    result = load_plugin(name).run(params=params, shared={})
    assert result["output_raster"] == params["output_raster"]
    assert len(calls) == 1


def test_example_offers_flaash_as_a_separate_commented_alternative():
    text = (Path(__file__).parents[1] / "configs/example.worldview.yml").read_text()
    recipe = yaml.safe_load(text)
    assert recipe["py6s_atmospheric_correction"]["plugin"] == "py6s_atmospheric_correction"
    assert "flaash_atmospheric_correction" not in recipe
    start = text.index("# flaash_atmospheric_correction:")
    end = text.index("# Orthorectify multispectral imagery", start)
    alternative = yaml.safe_load(re.sub(r"^# ", "", text[start:end], flags=re.MULTILINE))
    step = alternative["flaash_atmospheric_correction"]
    assert step["plugin"] == "flaash_atmospheric_correction"
    assert step["core:run"] is True
    assert "param:envi_engine_path" in step and "param:variables" in step
    assert "param:method" not in step and "param:custom_params" not in step
    assert "param:wavelengths_response_csv" not in step
