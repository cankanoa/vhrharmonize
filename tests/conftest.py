from __future__ import annotations

from pathlib import Path
import sys
from typing import Any, Callable

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


@pytest.fixture
def make_test_raster() -> Callable[..., Path]:
    def _make_test_raster(
        path: Path,
        *,
        count: int = 1,
        width: int = 4,
        height: int = 4,
        dtype: str = "uint16",
        nodata: float | int | None = 0,
        crs: str = "EPSG:4326",
        fill: int | float = 1,
        data: np.ndarray | None = None,
    ) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        if data is None:
            array = np.full((count, height, width), fill, dtype=dtype)
        else:
            array = np.asarray(data).astype(dtype)
            if array.ndim == 2:
                array = array[np.newaxis, ...]
        with rasterio.open(
            path,
            "w",
            driver="GTiff",
            width=array.shape[2],
            height=array.shape[1],
            count=array.shape[0],
            dtype=str(array.dtype),
            transform=from_origin(0, 4, 1, 1),
            crs=crs,
            nodata=nodata,
        ) as dst:
            dst.write(array)
        return path

    return _make_test_raster


@pytest.fixture
def make_worldview_bundle(tmp_path: Path, make_test_raster: Callable[..., Path]) -> Callable[..., dict[str, Any]]:
    def _make_worldview_bundle(
        basename: str = "17JUL05211635-M1BS-016445286010_01_P001",
        pan_basename: str = "17JUL05211635-P1BS-016445286010_01_P001",
    ) -> dict[str, Any]:
        scene_root = tmp_path / "scene"
        mul_dir = scene_root / "MUL"
        pan_dir = scene_root / "PAN"
        mul_tif = make_test_raster(mul_dir / f"{basename}.TIF", count=8)
        pan_tif = make_test_raster(pan_dir / f"{pan_basename}.TIF", count=1)
        imd_text = """version = \"28.3\";
productOrderId = \"016445286010_01_P001\";
BEGIN_GROUP = BAND_C
absCalFactor = 1.0;
effectiveBandwidth = 2.0;
END_GROUP = BAND_C
BEGIN_GROUP = BAND_B
absCalFactor = 2.0;
effectiveBandwidth = 4.0;
END_GROUP = BAND_B
BEGIN_GROUP = IMAGE_1
satId = \"WV03\";
TLCTime = 2017-07-05T21:16:43.287350Z;
meanSunAz = 75.2;
meanSunEl = 73.0;
meanSatAz = 79.7;
meanOffNadirViewAngle = 22.0;
meanProductRowGSD = 1.342;
meanProductColGSD = 1.463;
meanProductGSD = 1.401;
cloudCover = 0.814;
END_GROUP = IMAGE_1
END;"""
        (mul_dir / f"{basename}.IMD").write_text(imd_text, encoding="utf-8")
        (pan_dir / f"{pan_basename}.IMD").write_text(imd_text.replace("M1BS", "P1BS"), encoding="utf-8")
        return {
            "scene_root": scene_root,
            "mul_tif": mul_tif,
            "pan_tif": pan_tif,
            "mul_imd": mul_dir / f"{basename}.IMD",
            "pan_imd": pan_dir / f"{pan_basename}.IMD",
            "basename": basename,
            "pan_basename": pan_basename,
        }

    return _make_worldview_bundle
