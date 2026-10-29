"""Native merge modes accept workflow-owned concurrency settings."""

import pytest
from osgeo import gdal

from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import import_settings


@pytest.mark.parametrize("output_tiles", [False, True])
def test_workflow_merge_writes_single_file_or_tiles(tmp_path, make_test_raster, output_tiles):
    source = make_test_raster(tmp_path / "source.tif", width=32, height=32)
    output = tmp_path / ("tiles" if output_tiles else "mosaic.tif")
    workflow = Workflow({
        "shared": {
            "plugin": "shared",
            "core:run": True,
            "core:concurrent_processing": 1,
            "core:tile_threads": 1,
            "core:log_to_console": False,
        },
        "files": import_settings(source, tmp_path),
        "merge": {
            "plugin": "merge_rasters",
            "core:run": True,
            "core:require_outputs": "param:output_image_path",
            "param:input_images": "collect:mul",
            "param:output_image_path": str(output),
            "param:output_tiles": output_tiles,
            "param:window_size": 16,
        },
    })
    workflow.run()
    mosaic = output / "MergedImage.vrt" if output_tiles else output
    with gdal.Open(str(mosaic)) as dataset:
        assert dataset.RasterXSize == dataset.RasterYSize == 32
        assert (dataset.ReadAsArray() == 1).all()
    if output_tiles:
        assert len(list(output.glob("*.tif"))) == 4
