"""Native footprint meshing controls are exposed and forwarded by VHR."""
import importlib
from pathlib import Path
from unittest.mock import Mock

import yaml

from vhrharmonize.parameters import function_parameters
from vhrharmonize.workflow.registry import load_plugin


def test_native_mesh_controls_and_optional_footprints(monkeypatch):
    module = importlib.import_module('vhrharmonize.plugins.markov_triangles')
    parameters = function_parameters(module.markov_triangles)
    assert 'mesh_spacing' not in parameters
    assert parameters['mesh_max_triangle_area'].default == 10000
    assert parameters['mesh_min_triangle_angle'].default == 20
    assert parameters['input_polygons'].default is None
    native = Mock(return_value='out.gpkg')
    monkeypatch.setattr(module, '_upstream', lambda: native)
    params = dict(input_images=['a.tif', 'b.tif'], input_polygons='footprints.gpkg',
                  output_mask='out.gpkg', mesh_max_triangle_area=10000,
                  mesh_min_triangle_angle=20)
    assert load_plugin('markov_triangles').run(params=params, shared={}) == 'out.gpkg'
    native.assert_called_once_with(**params)


def test_worldview_mesh_defaults():
    root = Path(__file__).parents[1]
    paths = [root/'configs/example.worldview.yml',
             *sorted((root/'configs/drylands_groups').glob('*/*.worldview.yml'))]
    for path in paths:
        params = yaml.safe_load(path.read_text())['markov_triangles']
        assert 'param:mesh_spacing' not in params, path
        assert params['param:mesh_max_triangle_area'] == 10000, path
        assert params['param:mesh_min_triangle_angle'] == 20, path
        assert params['param:input_polygons'] == 'const:footprints', path
