"""Regression coverage for tiled merge outputs and directory downloads."""
import importlib
import subprocess
from argparse import Namespace

import pytest
import rasterio
import yaml


@pytest.mark.parametrize('mode', ['no', 'yes', 'validate'])
def test_directory_download_resumes_as_one_group(tmp_path, monkeypatch, make_test_raster, mode):
    slurm = importlib.import_module('vhrharmonize.slurm')
    remote = tmp_path / 'remote tiles'
    local = tmp_path / 'local tiles'
    make_test_raster(remote / 'tile[1].tif', fill=2)
    make_test_raster(remote / 'broken.tif', fill=3)
    make_test_raster(remote / '1' / 'missing.tif', fill=4)
    make_test_raster(local / 'tile[1].tif', fill=1)
    local.joinpath('broken.tif').write_text('broken')
    remote.joinpath('MergedImage.vrt').write_text('relative tile references')
    commands = []

    def run(command):
        commands.append(command)
        # Exercise real recursive rsync and exclusion semantics over local folders.
        command[-2] = str(remote) + '/'
        return subprocess.run(command, check=True, capture_output=True, text=True)

    monkeypatch.setattr(slurm, '_run_local_command', run)
    monkeypatch.setattr(slurm, '_remote_is_directory', lambda *args: True)
    config = tmp_path / 'staged.yml'
    config.write_text(yaml.safe_dump({
        'ssh_host': 'host', 'ssh_user': 'user', 'override_download_conflict': mode,
        'download_output_paths': {str(local): '/remote/tiles'},
    }))
    valid_before = local.joinpath('tile[1].tif').read_bytes()
    slurm.download_slurm_outputs(str(config))
    assert len(commands) == 1
    assert local.joinpath('1/missing.tif').read_bytes() == remote.joinpath('1/missing.tif').read_bytes()
    assert local.joinpath('MergedImage.vrt').exists()
    assert not local.joinpath('remote tiles').exists()
    assert local.joinpath('tile[1].tif').read_bytes() == (
        remote.joinpath('tile[1].tif').read_bytes() if mode == 'yes' else valid_before
    )
    assert local.joinpath('broken.tif').read_bytes() == (
        b'broken' if mode == 'no' else remote.joinpath('broken.tif').read_bytes()
    )


@pytest.mark.parametrize('steps', [['merge'], ['align', 'merge']])
def test_pipeline_runs_real_tiled_merge(tmp_path, make_test_raster, steps):
    rad = importlib.import_module('vhrharmonize.preprocess.radiometric_normalization')
    raster = make_test_raster(tmp_path / 'input.tif', width=32, height=32)
    output = tmp_path / 'tiles'
    worldview = importlib.import_module('vhrharmonize.cli.worldview')
    kwargs = worldview._build_radiometric_kwargs(Namespace(
        radiometric_normalization_kwargs_json=None, overview_scales=[2],
        match_merge_rasters_build_overviews=True,
    ))
    result = rad.radiometric_normalization(
        [str(raster)], str(output), steps=steps,
        shared_image_threads=1, shared_io_threads=1, shared_tile_threads=1,
        merge_rasters_output_tiles=True, shared_window_size=16,
        merge_rasters_custom_tiles_csv='tiles.csv', **kwargs,
    )
    assert result == str(output)
    assert len(list(output.glob('*.tif'))) == 4
    assert (output / 'MergedImage.vrt').exists()
    assert (output / 'tiles.csv').exists()
    assert (output / '1' / 'MergedImage.vrt').exists()
    with rasterio.open(output / 'MergedImage.vrt') as mosaic:
        assert mosaic.read(1).shape == (32, 32)
        assert (mosaic.read(1) == 1).all()


def test_tiled_group_plan_keeps_one_folder_mapping(tmp_path):
    slurm = importlib.import_module('vhrharmonize.slurm')
    result = slurm._collect_group_by_basename_output_downloads(
        {'merged_tiles': 'auto:*.tif'}, local_temp_root=str(tmp_path / 'temp'),
        local_output_root=str(tmp_path / 'output'), remote_output_dir='/remote/output',
    )
    assert result == {str(tmp_path / 'output/merged_tiles'): '/remote/output/merged_tiles'}


@pytest.mark.parametrize('override', [[2, 4], [], None])
def test_shared_scales_override_workflow_overviews(override):
    worldview = importlib.import_module('vhrharmonize.cli.worldview')
    kwargs = worldview._build_radiometric_kwargs(Namespace(
        radiometric_normalization_kwargs_json=None, overview_scales=[2, 4, 8],
        match_shared_window_scales=override,
    ))
    assert kwargs['shared_window_scales'] == override
    assert 'merge_rasters_window_scales' not in kwargs


def test_match_yaml_options_reach_pipeline_with_nulls_and_tuple_types(tmp_path, monkeypatch):
    rad = importlib.import_module('vhrharmonize.preprocess.radiometric_normalization')
    from spectralmatch.types_and_validation import Match
    calls = []

    def pipeline(**kwargs):
        calls.append(kwargs)
        Match._validate_local_block_adjustment(
            number_of_blocks=kwargs['local_block_adjustment_number_of_blocks'],
            load_block_maps=kwargs['local_block_adjustment_load_block_maps'],
        )
        return {'output': kwargs['shared_output_image_path']}

    monkeypatch.setattr(rad, 'spectralmatch_pipeline', pipeline)
    options = yaml.safe_load('''
global_regression_pif_method: flood_from_match_points
global_regression_pif_max_samples: null
local_block_adjustment_number_of_blocks: [2, 3]
local_block_adjustment_load_block_maps: [null, [one.tif, two.tif]]
shared_window_scales: null
''')
    worldview = importlib.import_module('vhrharmonize.cli.worldview')
    kwargs = worldview._build_radiometric_kwargs(Namespace(
        radiometric_normalization_kwargs_json=None, overview_scales=[2, 4],
        **{'match_' + key: value for key, value in options.items()},
    ))
    rad.radiometric_normalization(['input.tif'], str(tmp_path / 'out.tif'), **kwargs)
    assert len(calls) == 1
    assert calls[0]['global_regression_pif_method'] == 'flood_from_match_points'
    assert calls[0]['global_regression_pif_max_samples'] is None
    assert calls[0]['shared_window_scales'] is None
    assert calls[0]['local_block_adjustment_load_block_maps'] == (None, ['one.tif', 'two.tif'])


@pytest.mark.parametrize('final_file_overviews', [False, True])
def test_pipeline_overview_flags_are_independent(final_file_overviews):
    worldview = importlib.import_module('vhrharmonize.cli.worldview')
    args = Namespace(
        radiometric_normalization_kwargs_json=None,
        calculate_overviews_radiometric_normalization=final_file_overviews,
    )
    kwargs = worldview._build_radiometric_kwargs(args)
    assert not any(key.endswith('_build_overviews') for key in kwargs)

    flags = {
        'joint_coregistration_build_overviews': True,
        'global_regression_build_overviews': False,
        'local_block_adjustment_build_overviews': True,
        'merge_rasters_build_overviews': False,
    }
    for key, value in flags.items():
        setattr(args, 'match_' + key, value)
    kwargs = worldview._build_radiometric_kwargs(args)
    assert {key: kwargs[key] for key in flags} == flags
