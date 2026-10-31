import inspect
import json
from pathlib import Path
import numpy as np
import pytest
import yaml
from osgeo import gdal

from vhrharmonize.parameters import function_parameters
from vhrharmonize.plugins.joint_coregistration import joint_coregistration
from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import import_settings


def test_native_matching_and_fixed_parameters_exposed():
    parameters = function_parameters(joint_coregistration)
    assert 'fixed_images' in parameters
    assert 'elastix' in str(parameters['moving_tie_method'].annotation)
    assert 'elastix' in str(parameters['fixed_tie_method'].annotation)
    assert parameters['moving_tie_method'].default == 'orb'
    assert parameters['fixed_tie_method'].default == 'elastix'
    assert 'tie_feature_method' not in parameters
    assert 'tie_fixed_feature_method' not in parameters
    assert parameters['input_images_bands'].default == [1]
    assert parameters['fixed_images_bands'].default == [1]


def test_lidar_recipe_uses_separate_reference_paths_and_only_moving_outputs(tmp_path, make_test_raster):
    recipe = yaml.safe_load((Path(__file__).parents[1] / 'configs/example.worldview.yml').read_text())
    step = recipe['lidar_joint_coregistration']
    assert step['param:input_images'] == 'var:current_image_paths'
    assert step['param:output_images'] == 'var:current_image_paths'
    assert step['param:fixed_images'] == 'expr:[const.lidar_reference]'
    source = make_test_raster(tmp_path / 'source_cloudmasked.tif', crs='EPSG:3857')
    reference = make_test_raster(tmp_path / 'intensity.tif', crs='EPSG:3857')
    discover = import_settings(source,tmp_path)
    discover.update({'var:current_image_paths':'returned:file_path','const:reference_path':str(reference),
                     'var:relative_output_dir':str(tmp_path/'processed')})
    config = {'shared': {'plugin':'shared','core:run':True,'core:save_statistics_path':None},
              'discover':discover,'lidar':step}
    workflow = Workflow(config)
    workflow.plan()
    assert len(workflow.records) == 1
    node = next(n for n in workflow.nodes if n.step['name']=='lidar')
    assert node.params['input_images'] == str(source)
    assert node.params['fixed_images'] == [str(reference)]
    assert node.params['input_images_bands'] == node.params['fixed_images_bands'] == [1]
    assert node.params['output_images'] == str(tmp_path/'processed'/f"source_cloudmasked{step['var:suffix']}.tif")
    assert node.params['tie_save_path'] == str(tmp_path/'processed'/f"source_cloudmasked{step['var:suffix']}_ties.json")


@pytest.mark.parametrize('explicit_lists', [False, True])
def test_recipe_pairs_each_scene_with_fixed_lidar(tmp_path, explicit_lists):
    from spectralmatch.joint_coregistration.joint_coregistration import _save_ties
    source_dir = tmp_path/'scenes'
    source_dir.mkdir()
    reference = tmp_path/'lidar_intensity.tif'
    def write(path, offset):
        ds = gdal.GetDriverByName('GTiff').Create(str(path),24,24,1,gdal.GDT_UInt16)
        ds.SetGeoTransform((offset,1,0,24,0,-1))
        from osgeo import osr
        crs = osr.SpatialReference()
        crs.ImportFromEPSG(3857)
        ds.SetProjection(crs.ExportToWkt())
        ds.GetRasterBand(1).WriteArray(np.ones((24,24),dtype=np.uint16))
        ds = None
    write(reference,0)
    reference_bytes = reference.read_bytes()
    for name in ('a','b'):
        write(source_dir/f'{name}.tif',2)
        _save_ties(str(tmp_path/f'{name}.json'),{(name,'lidar_intensity'):np.array([[10,10,10,10]])})
    recipe = yaml.safe_load((Path(__file__).parents[1]/'configs/example.worldview.yml').read_text())
    step = dict(recipe['lidar_joint_coregistration'])
    step.update({'core:run':True, 'core:require_outputs':'param:output_images', 'param:skip_single_image':True,
                 'param:local_model':'none',
                 'param:tie_load_path':"expr:'"+str(tmp_path)+"/' & var.basename & '.json'"})
    if explicit_lists:
        step['param:input_images'] = ['var:current_image_paths']
        step['param:output_images'] = ['var:current_image_paths']
    discover = import_settings(source_dir/'*.tif',tmp_path)
    discover.update({'var:current_image_paths':'returned:file_path','const:reference_path':str(reference),
                     'var:relative_output_dir':str(tmp_path/'processed')})
    workflow = Workflow({'shared':{'plugin':'shared','core:run':True,'core:log_to_console':False,
                                  'core:save_statistics_path':None},
                         'discover':discover,'lidar':step})
    workflow.run()
    for name in ('a','b'):
        output = tmp_path/'processed'/f"{name}{step['var:suffix']}.tif"
        assert output.exists()
        assert gdal.Open(str(output)).GetGeoTransform()[0] == 0
        assert not (tmp_path/'temp/lidar_joint_coregistration').exists()
        assert reference.read_bytes() == reference_bytes
        ties = json.loads((tmp_path/'processed'/f"{name}{step['var:suffix']}_ties.json").read_text())
        assert len(ties['tie_points']) == 1
