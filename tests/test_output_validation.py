from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

from osgeo import gdal
import pytest
import rasterio

from vhrharmonize.io import validation


@pytest.mark.parametrize("use_exceptions", [False, True])
def test_raster_validity_handles_unreadable_tiff_directory(tmp_path, make_test_raster, use_exceptions):
    corrupt_path = tmp_path / "corrupt.tif"
    # TIFF header pointing to a directory that was never written.
    corrupt_path.write_bytes(b"II\x2a\x00\x08\x02\x00\x00")
    valid_path = make_test_raster(tmp_path / "valid.tif")

    with gdal.ExceptionMgr(useExceptions=use_exceptions):
        is_valid, reason = validation._gdal_raster_is_valid(str(corrupt_path))
        assert not is_valid
        assert reason
        assert validation._gdal_raster_is_valid(str(valid_path)) == (True, None)


@pytest.mark.parametrize("check_validity", [False, True])
@pytest.mark.parametrize(
    "contents, valid",
    [
        (b'{"water_vapor": 2.5}', True),
        (b"{}", True),
        (b"[]", True),
        (b"null", True),
        (b"", False),
        (b'{"water_vapor":', False),
        (b"{} trailing text", False),
        (b'"\xff"', False),
    ],
)
def test_json_output_reuse_checks_syntax_only_when_enabled(tmp_path, contents, valid, check_validity):
    output_path = tmp_path / "atmosphere.JSON"
    output_path.write_bytes(contents)

    assert validation._existing_outputs_are_reusable(
        [str(output_path)],
        check_validity=check_validity,
        validity_check_grid_size=0,
        log_to_console=False,
        step="fetch_atmosphere",
    ) is (valid or not check_validity)
