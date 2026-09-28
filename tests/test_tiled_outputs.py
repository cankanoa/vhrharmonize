"""Regression coverage for tiled merge outputs and directory downloads."""

import importlib
import subprocess
from argparse import Namespace

import pytest
import rasterio
import yaml


@pytest.mark.parametrize("mode", ["no", "yes", "validate"])
def test_directory_download_resumes_as_one_group(tmp_path, monkeypatch, make_test_raster, mode):
    slurm = importlib.import_module("vhrharmonize.slurm")
    remote = tmp_path / "remote tiles"
    local = tmp_path / "local tiles"
    make_test_raster(remote / "tile[1].tif", fill=2)
    make_test_raster(remote / "broken.tif", fill=3)
    make_test_raster(remote / "1" / "missing.tif", fill=4)
    make_test_raster(local / "tile[1].tif", fill=1)
    local.joinpath("broken.tif").write_text("broken")
    remote.joinpath("MergedImage.vrt").write_text("relative tile references")
    commands = []

    def run(command):
        commands.append(command)
        # Exercise real recursive rsync and exclusion semantics over local folders.
        command[-2] = str(remote) + "/"
        return subprocess.run(command, check=True, capture_output=True, text=True)

    monkeypatch.setattr(slurm, "_run_local_command", run)
    monkeypatch.setattr(slurm, "_remote_is_directory", lambda *args: True)
    config = tmp_path / "staged.yml"
    config.write_text(
        yaml.safe_dump(
            {
                "ssh_host": "host",
                "ssh_user": "user",
                "override_download_conflict": mode,
                "download_output_paths": {str(local): "/remote/tiles"},
            }
        )
    )
    valid_before = local.joinpath("tile[1].tif").read_bytes()
    slurm.download_slurm_outputs(str(config))
    assert len(commands) == 1
    assert (
        local.joinpath("1/missing.tif").read_bytes()
        == remote.joinpath("1/missing.tif").read_bytes()
    )
    assert local.joinpath("MergedImage.vrt").exists()
    assert not local.joinpath("remote tiles").exists()
    assert local.joinpath("tile[1].tif").read_bytes() == (
        remote.joinpath("tile[1].tif").read_bytes() if mode == "yes" else valid_before
    )
    assert local.joinpath("broken.tif").read_bytes() == (
        b"broken" if mode == "no" else remote.joinpath("broken.tif").read_bytes()
    )


@pytest.mark.parametrize("align_first", [False, True])
def test_individual_functions_run_real_tiled_merge(tmp_path, make_test_raster, align_first):
    from vhrharmonize.plugins.merge_rasters import merge_rasters
    from vhrharmonize.plugins.align_rasters import align_rasters

    raster = make_test_raster(tmp_path / "input.tif", width=32, height=32)
    images = [str(raster)]
    if align_first:
        images = align_rasters(images, [str(tmp_path / "aligned.tif")], image_threads=1)
    output = tmp_path / "tiles"
    result = merge_rasters(
        images,
        str(output),
        image_threads=1,
        io_threads=1,
        tile_threads=1,
        output_tiles=True,
        window_size=16,
        custom_tiles_csv="tiles.csv",
        build_overviews=True,
        window_scales=[2],
    )
    assert result == str(output)
    assert len(list(output.glob("*.tif"))) == 4
    assert (output / "MergedImage.vrt").exists()
    assert (output / "tiles.csv").exists()
    assert (output / "1" / "MergedImage.vrt").exists()
    with rasterio.open(output / "MergedImage.vrt") as mosaic:
        assert mosaic.read(1).shape == (32, 32)
        assert (mosaic.read(1) == 1).all()
