"""Native overlap placement and pair quotas survive VHR argument resolution."""
import importlib
from pathlib import Path
from unittest.mock import Mock

import pytest
import yaml

from vhrharmonize.parameters import function_parameters
from vhrharmonize.workflow.registry import load_plugin


@pytest.mark.parametrize('name,default', [('joint_coregistration', 'grid'), ('global_regression', 'overlap_spread')])
def test_native_sampling_defaults_and_quota_forwarding(name, default, monkeypatch):
    module = importlib.import_module('vhrharmonize.plugins.' + name)
    parameters = function_parameters(getattr(module, name))
    assert parameters['tie_sampling_method'].default == default
    assert parameters['tie_max_samples'].default == '1/30crs'
    assert parameters['tie_min_samples'].default == '1/500crs'
    native = Mock(return_value=['out.tif'])
    monkeypatch.setattr(module, '_upstream', lambda: native)
    params = dict(input_images=['a.tif', 'b.tif'], output_images=['out.tif'],
                  tie_sampling_method='overlap_spread', tie_max_samples='5/3crs',
                  tie_min_samples=4, tie_max_attempts=100)
    assert load_plugin(name).run(params=params, shared={}) == ['out.tif']
    native.assert_called_once_with(**params)


def test_worldview_pair_quotas_are_separate_from_pif_pixels():
    root = Path(__file__).parents[1]
    paths = [root/'configs/example.worldview.yml', *sorted((root/'configs/drylands_groups').glob('*/*.worldview.yml'))]
    for path in paths:
        config = yaml.safe_load(path.read_text())
        assert config['joint_coregistration']['param:tie_sampling_method'] == 'grid', path
        matching = config['global_regression']
        assert matching['param:tie_sampling_method'] == 'overlap_spread', path
        assert matching['param:tie_max_samples'] == '1/30crs', path
        assert matching['param:tie_min_samples'] == '1/500crs', path
        assert matching['param:pif_region_radius'] == 30, path
