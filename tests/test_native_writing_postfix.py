"""The runner passes final paths and a string suffix to native functions."""
from pathlib import Path

import pytest
from osgeo import gdal

from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.output_files import writing_path
from workflow_helpers import import_settings


def test_polygonize_reuses_original_vector_with_core_validation_off(tmp_path, make_test_raster, monkeypatch):
    from spectralmatch.seamline import footprints
    source = make_test_raster(tmp_path / 'source.tif', width=32, height=32)
    output = tmp_path / 'footprint.gpkg'
    config = {
        'shared': {'plugin': 'shared', 'core:run': True, 'core:check_validity': False,
                   'core:log_to_console': False, 'core:show_progress': False,
                   'core:save_statistics_path': None, 'core:load_statistics_path': None},
        'files': import_settings(source, tmp_path),
        'polygonize': {'plugin': 'create_footprints', 'core:run': True, 'core:scope': 'var',
                      'core:require_outputs': True, 'param:input_images': 'expr:[var.mul]',
                      'param:output_polygons': str(output), 'param:resume_from_outputs': 'validate'},
    }
    Workflow(config).run()
    previous = output.read_bytes()
    def unexpected(*args, **kwargs):
        pytest.fail('Existing footprints must be found using the original path')
    monkeypatch.setattr(footprints, '_footprint_from_image', unexpected)
    Workflow(config).run()
    assert output.read_bytes() == previous
    assert not Path(writing_path(str(output), '_writing')).exists()


@pytest.mark.parametrize('swap', [False, True])
def test_native_publication_validates_and_builds_overviews_before_rename(tmp_path, make_test_raster, monkeypatch, swap):
    from vhrharmonize.workflow.output_files import OutputPolicy
    source = make_test_raster(tmp_path / 'source.tif', width=32, height=32)
    output = str(tmp_path / 'mosaic.tif')
    original = OutputPolicy.finish
    publications = []
    def finish(self, outputs, mapping):
        if output in mapping:
            assert mapping[output] == writing_path(output, '_writing')
            assert not Path(output).exists()
            assert self.check_validity == 'around'
            publications.append(output)
        return original(self, outputs, mapping)
    monkeypatch.setattr(OutputPolicy, 'finish', finish)
    config = {
        'shared': {'plugin': 'shared', 'core:run': True, 'core:check_validity': 'around',
                   'core:show_progress': False, 'core:save_statistics_path': None,
                   'core:load_statistics_path': None, 'param:window_scales': [2, 4]},
        'merge': {'plugin': 'merge_rasters', 'core:run': True, 'core:require_outputs': True,
                  'core:calculate_overviews': True, 'param:input_images': [str(source)],
                  'param:output_image_path': output},
    }
    if swap:
        config['shared']['core:swap_temp_dir'] = str(tmp_path / 'swap')
        config['merge']['core:swap_for_temp_file'] = ['output_image_path']
    Workflow(config).run()
    assert publications == [output]
    with gdal.Open(output) as dataset:
        assert dataset.GetRasterBand(1).GetOverviewCount() == 2
