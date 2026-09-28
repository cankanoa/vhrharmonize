import numpy as np
from osgeo import gdal
import pytest
from tifffile import TiffFile, TiffFileError, TiffWriter, imwrite

from vhrharmonize.io import validation


@pytest.mark.parametrize("bigtiff", [False, True])
@pytest.mark.parametrize("tiled", [False, True])
@pytest.mark.parametrize("byteorder", ["<", ">"])
def test_tiff_bounds_detect_unreadable_middle_block(tmp_path, bigtiff, tiled, byteorder):
    path = tmp_path / "middle_block.tif"
    imwrite(
        path, np.ones((1024, 1024), dtype=np.uint16), bigtiff=bigtiff,
        byteorder=byteorder, metadata=None,
        **({"tile": (128, 128)} if tiled else {"rowsperstrip": 32}),
    )
    assert validation._gdal_raster_is_valid(str(path)) == (True, None)
    with TiffFile(path, mode="r+") as tif:
        page = tif.pages[0]
        offsets = list(page.dataoffsets)
        # Keep both corner windows readable; only an interior block is damaged.
        offsets[len(offsets) // 2 + 1] = (1 << 32) + 123 if bigtiff else path.stat().st_size - 8
        page.tags["TileOffsets" if tiled else "StripOffsets"].overwrite(offsets)

    dataset = gdal.Open(str(path))
    band = dataset.GetRasterBand(1)
    assert band.ReadRaster(0, 0, 256, 256)
    assert band.ReadRaster(768, 768, 256, 256)
    band = dataset = None

    before = path.read_bytes()
    valid, reason = validation._gdal_raster_is_valid(str(path))

    assert not valid
    assert "Truncated TIFF:" in reason
    assert "byte count" in reason and "file size" in reason
    assert path.read_bytes() == before  # Validation itself never deletes or rewrites.


@pytest.mark.parametrize("subifd", [False, True])
def test_tiff_bounds_detect_truncated_internal_overview(tmp_path, subifd):
    path = tmp_path / "overview.tif"
    with TiffWriter(path, bigtiff=True) as tif:
        tif.write(np.ones((512, 512), dtype=np.uint8), metadata=None, subifds=1 if subifd else None)
        tif.write(np.ones((256, 256), dtype=np.uint8), metadata=None, subfiletype=1)
    with path.open("r+b") as stream:
        stream.truncate(path.stat().st_size - 8)

    valid, reason = validation._gdal_raster_is_valid(str(path))

    assert not valid
    assert "Truncated TIFF:" in reason


@pytest.mark.parametrize("bigtiff", [False, True])
def test_tiff_bounds_allow_sparse_blocks_and_trailing_bytes(tmp_path, bigtiff):
    path = tmp_path / "sparse.tif"
    dataset = gdal.GetDriverByName("GTiff").Create(
        str(path), 512, 512, 1, gdal.GDT_Byte,
        options=["TILED=YES", "SPARSE_OK=YES", f"BIGTIFF={'YES' if bigtiff else 'NO'}"],
    )
    dataset.GetRasterBand(1).SetNoDataValue(0)
    dataset = None
    with TiffFile(path) as tif:
        assert all(offset == count == 0 for offset, count in zip(
            tif.pages[0].dataoffsets, tif.pages[0].databytecounts,
        ))
    with path.open("ab") as stream:
        stream.write(b"extra metadata after image data")

    assert validation._gdal_raster_is_valid(str(path)) == (True, None)


def test_tiff_metadata_inspection_failure_falls_back_to_gdal(tmp_path, make_test_raster, monkeypatch):
    path = make_test_raster(tmp_path / "valid.tif")

    def unsupported(*args, **kwargs):
        raise TiffFileError("unsupported TIFF metadata")

    monkeypatch.setattr(validation, "TiffFile", unsupported)

    assert validation._gdal_raster_is_valid(str(path)) == (True, None)


def test_non_tiff_rasters_keep_existing_validation(tmp_path, make_test_raster, monkeypatch):
    source = make_test_raster(tmp_path / "source.tif")
    path = tmp_path / "image.vrt"
    dataset = gdal.Translate(str(path), str(source), format="VRT")
    dataset = None

    def unexpected(*args, **kwargs):
        raise AssertionError("TIFF inspection must not run on a VRT")

    monkeypatch.setattr(validation, "TiffFile", unexpected)

    assert validation._gdal_raster_is_valid(str(path)) == (True, None)
