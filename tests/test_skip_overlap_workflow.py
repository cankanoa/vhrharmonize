"""VHR exposes the native opt-in skip without adding error handling in core."""
import importlib
from pathlib import Path
from unittest.mock import Mock

import pytest
import yaml

from vhrharmonize.parameters import function_parameters
from vhrharmonize.workflow.registry import load_plugin


@pytest.mark.parametrize("name", ["joint_coregistration", "global_regression"])
@pytest.mark.parametrize("enabled", [False, True])
def test_native_overlap_skip_option_is_exposed_and_forwarded(name, enabled, monkeypatch):
    module = importlib.import_module("vhrharmonize.plugins." + name)
    parameters = function_parameters(getattr(module, name))
    assert parameters["skip_overlap_when_not_enough_matches"].default is False
    native = Mock(return_value=["out.tif"])
    monkeypatch.setattr(module, "_upstream", lambda: native)
    result = load_plugin(name).run(
        params={"input_images": ["a.tif", "b.tif"], "output_images": ["out.tif"],
                "skip_overlap_when_not_enough_matches": enabled}, shared={})
    assert result == ["out.tif"]
    assert native.call_args.kwargs["skip_overlap_when_not_enough_matches"] is enabled
