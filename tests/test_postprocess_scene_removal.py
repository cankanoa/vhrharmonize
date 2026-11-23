"""Footprint smoothing selects downstream scenes on fresh and cached runs."""
from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
from osgeo import gdal, osr
from shapely.geometry import box

from spectralmatch import postprocess_footprints
from vhrharmonize.plugins.postprocess_footprints import PostprocessFootprints
from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import install_function


def recipe(monkeypatch, tmp_path, survivors):
    records, geometries = [], []
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3857)
    for i in range(3):
        path = str(tmp_path / f'image{i}_coreg_global_local.tif')
        raster = gdal.GetDriverByName('GTiff').Create(path, 64, 64, 1, gdal.GDT_Int16)
        raster.SetProjection(srs.ExportToWkt())
        raster.SetGeoTransform((0, 1, 0, 64, 0, -1))
        raster.GetRasterBand(1).SetNoDataValue(-9999)
        raster.GetRasterBand(1).Fill(i + 1)
        raster = None
        # Stable IDs intentionally differ from raster names. Footprints predate
        # radiometric matching, and multiple features may identify one scene.
        records.append({'scene_id': f'record-{i}', 'image': path})
        geometries.append(box(4 + i, 4, 60, 60) if i < survivors else box(4, 4, 8, 60))
    source, output = tmp_path / 'aligned.gpkg', tmp_path / 'processed.gpkg'
    frame = gpd.GeoDataFrame({
        'scene_id': [r['scene_id'] for r in records],
        'image': [f'image{i}_coreg' for i in range(3)],
        'image_path': [r['image'] for r in records],
    }, geometry=geometries, crs=3857)
    frame.to_file(source, layer='aligned')
    install_function(monkeypatch, 'postprocess_inputs', lambda: records, scope='aggregate',
                     var_records_return='$', var_id_return='scene_id')
    params = dict(input_polygons=str(source), output_polygons=str(output),
                  input_layer='aligned', output_layer='processed',
                  hole_edge_distance=0, hole_to_hole_distance=0,
                  simplify_smoothing_radius=4, simplify_tolerance=0)
    config = {
        'shared': {'plugin': 'shared', 'core:run': True, 'core:log_to_console': False,
                   'core:show_progress': False, 'core:save_statistics_path': None,
                   'core:load_statistics_path': None, 'core:concurrent_processing': 1,
                   'core:check_validity': None, 'param:tile_threads': 1,
                   'param:io_threads': 1},
        'input': {'plugin': 'postprocess_inputs', 'core:run': True},
        'postprocess': {'plugin': 'postprocess_footprints', 'core:run': True,
                        **{'param:' + key: value for key, value in params.items()}},
        'seams': {'plugin': 'markov_triangles', 'core:run': True,
                  'param:input_images': 'collect:image',
                  'param:input_footprints': str(output), 'param:input_layer': 'processed',
                  'param:output_mask': str(tmp_path / 'seams.gpkg')},
        'mask': {'plugin': 'mask_rasters', 'core:run': True, 'core:scope': 'var',
                 'core:requires': str(tmp_path / 'seams.gpkg'),
                 'param:input_images': 'expr:[var.image]',
                 'var:masked': "expr:var.image & '_masked.tif'",
                 'param:output_images': 'expr:[var.masked]',
                 'param:vector_mask': ['include', str(tmp_path / 'seams.gpkg'), 'image'],
                 'param:include_touched_pixels': True},
        'merge': {'plugin': 'merge_rasters', 'core:run': True, 'core:require_outputs': True,
                  'param:input_images': 'collect:masked',
                  'param:output_image_path': str(tmp_path / 'mosaic.tif')},
    }
    return config, params


@pytest.mark.parametrize('cached', [False, True])
@pytest.mark.parametrize('survivors', [1, 2])
def test_postprocess_selects_surviving_scenes_through_mosaic(monkeypatch, tmp_path, cached, survivors):
    config, params = recipe(monkeypatch, tmp_path, survivors)
    if cached:
        postprocess_footprints(**params)
        Path(params['input_polygons']).unlink()  # Typical proactive cleanup.
        original = Path(params['output_polygons']).read_bytes()
    workflow = Workflow(config)
    workflow.counts()
    assert len(workflow.records) == 3  # Planning must not remove scenes.
    result = workflow.run()
    assert [r['id'] for r in result] == [f'record-{i}' for i in range(survivors)]
    assert workflow.errors == []
    seams = gpd.read_file(tmp_path / 'seams.gpkg')
    assert set(seams.image).issubset({f'image{i}_coreg_global_local' for i in range(survivors)})
    with gdal.Open(str(tmp_path / 'mosaic.tif')) as mosaic:
        assert np.any(mosaic.ReadAsArray() != -9999)
    assert not (tmp_path / 'image2_coreg_global_local.tif_masked.tif').exists()
    assert (tmp_path / 'image2_coreg_global_local.tif').exists()
    if cached:
        assert Path(params['output_polygons']).read_bytes() == original


@pytest.mark.parametrize('cached', [False, True])
def test_no_coverage_stops_before_seamlines(monkeypatch, tmp_path, cached):
    config, params = recipe(monkeypatch, tmp_path, 0)
    if cached:
        gpd.GeoDataFrame({'scene_id': []}, geometry=[], crs=3857).to_file(
            params['output_polygons'], layer='processed')
        Path(params['input_polygons']).unlink()
    with pytest.raises(ValueError, match='No coverage remains'):
        Workflow(config).run()
    assert not (tmp_path / 'seams.gpkg').exists()
    assert not (tmp_path / 'mosaic.tif').exists()


@pytest.mark.parametrize('ids,error', [(None, 'requires a scene_id'),
                                      (['unknown'], 'outside the current workflow'),
                                      ([None], 'missing scene_id')])
def test_bad_identifiers_are_not_treated_as_removed_scenes(tmp_path, ids, error):
    path = tmp_path / 'processed.gpkg'
    gpd.GeoDataFrame({} if ids is None else {'scene_id': ids},
                     geometry=[box(0, 0, 10, 10)], crs=3857).to_file(path, layer='footprints')
    with pytest.raises(ValueError, match=error):
        PostprocessFootprints().scene_ids_to_remove(
            params={'output_polygons': str(path)}, returned=None, scene_ids=['a'])


def test_surviving_component_keeps_scene_and_incomplete_output_is_not_reused(tmp_path):
    path = tmp_path / 'processed.gpkg'
    gpd.GeoDataFrame({'scene_id': ['a', 'a']},
                     geometry=[box(0, 0, 10, 10), box(20, 20, 30, 30)], crs=3857).to_file(
                         path, layer='footprints')
    plugin = PostprocessFootprints()
    params = {'output_polygons': str(path)}
    assert plugin.scene_ids_to_remove(params=params, returned=None, scene_ids=['a', 'b']) == ['b']
    assert plugin.can_reuse_outputs(params=params)
    Path(str(path) + '.incomplete').touch()
    assert not plugin.can_reuse_outputs(params=params)
