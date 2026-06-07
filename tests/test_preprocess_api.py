from __future__ import annotations

import sys
import importlib
from pathlib import Path
from types import ModuleType, SimpleNamespace

import numpy as np
import rasterio

import vhrharmonize.preprocess.atmospheric_correction as atmos_mod
import vhrharmonize.preprocess.fetch_external_data as fetch_mod
import vhrharmonize.preprocess.radiometric_normalization as rad_mod
from vhrharmonize.preprocess.alignment import align_image_pair
from vhrharmonize.preprocess.atmospheric_correction import (
    atmospheric_correction,
    build_flaash_kwargs_from_standardized_metadata,
    build_py6s_kwargs_from_standardized_metadata,
    convert_flaash_params_paths_for_windows,
    init_envi_engine,
    parallel_flaash,
    run_flaash,
    run_flaash_wrapper,
    run_py6s,
    validate_flaash_params,
    wsl_path_to_windows_for_envi,
)
from vhrharmonize.preprocess.cloudmasking import (
    apply_binary_cloud_mask_to_image,
    cloudmask_raster,
    create_cloud_mask_with_omnicloudmask,
)
from vhrharmonize.preprocess.fetch_external_data import (
    download_opentopography_dem_for_bbox,
    fetch_modis_water_vapor_for_bbox,
    fetch_power_atmosphere_for_bbox,
    init_ee_client,
)
from vhrharmonize.preprocess.helpers import log
from vhrharmonize.preprocess.radiometric_normalization import radiometric_normalization


def _metadata() -> SimpleNamespace:
    return SimpleNamespace(
        sun_zenith=17.0,
        sun_azimuth=75.0,
        view_zenith=22.0,
        line_of_sight_zenith=158.0,
        line_of_sight_azimuth=80.0,
        dn_to_radiance_factors=(0.5, 1.0),
        dn_to_radiance_gains=(2.0, 3.0),
        dn_to_radiance_offsets=(-1.0, -2.0),
        resolve_scene_datetime=lambda: __import__("datetime").datetime(2017, 7, 5),
    )


def test_log(capsys) -> None:
    log("hello", enabled=True, step="step", scene_basename="scene")
    assert "[step_scene] hello" in capsys.readouterr().out


def test_align_image_pair(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setitem(sys.modules, "coregix", SimpleNamespace(align_image_pair=lambda **kwargs: SimpleNamespace(output_image_path=kwargs["output_image_path"])))
    output = align_image_pair("moving.tif", "fixed.tif", str(tmp_path / "out.tif"))
    assert output.output_image_path.endswith("out.tif")


def test_py6s_helpers_and_run(monkeypatch, tmp_path: Path) -> None:
    kwargs = build_py6s_kwargs_from_standardized_metadata(_metadata(), ground_elevation_km=1.0, atmosphere_profile="user", aerosol_profile="maritime", aot550=0.2, visibility_km=None, water_vapor=2.5, ozone=0.3, sixs_executable=None, output_scale_factor=10000.0, output_dtype="int16", use_imd_radiance_calibration=True, use_worldview_gain_offset_adjustment=True)
    assert kwargs["dn_to_radiance_factors"] == [1.0, 3.0]
    monkeypatch.setattr(atmos_mod.Py6SCorrector, "run", lambda self, input_raster, output_raster, **kw: output_raster)
    result = run_py6s("in.tif", str(tmp_path / "out.tif"), _metadata(), ground_elevation_km=1.0, atmosphere_profile="user", aerosol_profile="maritime", aot550=0.2, visibility_km=None, water_vapor=2.5, ozone=0.3, sixs_executable=None, output_scale_factor=1.0, output_dtype="int16", use_imd_radiance_calibration=False, use_worldview_gain_offset_adjustment=False)
    assert result.output_raster.endswith("out.tif")


def test_flaash_helpers_and_dispatch(monkeypatch, tmp_path: Path) -> None:
    fake_geospatial = ModuleType("vhrharmonize.io.geospatial")
    fake_geospatial.get_image_percentile_value = lambda *args, **kwargs: 100.0
    monkeypatch.setitem(sys.modules, "vhrharmonize.io.geospatial", fake_geospatial)
    kwargs = build_flaash_kwargs_from_standardized_metadata("in.tif", "dem.tif", "footprint.gpkg", _metadata(), "out.tif", dem_ground_percentile=50.0, modtran_atm="Mid-Latitude Summer", modtran_aer="Rural", use_aerosol="AUTO", default_visibility=10.0)
    assert kwargs["GROUND_ELEVATION"] == 0.1
    assert validate_flaash_params({"MODTRAN_ATM": "x"})["MODTRAN_ATM"] == "x"
    assert wsl_path_to_windows_for_envi("/mnt/c/test/file") == "C:\\test\\file"
    assert convert_flaash_params_paths_for_windows({"INPUT_RASTER": {"url": "/mnt/c/x", "factory": "URLRaster"}, "OUTPUT_RASTER_URI": "/mnt/c/y"})["OUTPUT_RASTER_URI"] == "C:\\y"
    fake_envipyengine = ModuleType("envipyengine")
    fake_envipyengine.Engine = lambda *_: SimpleNamespace(tasks=lambda: None)
    fake_envipyengine_config = ModuleType("envipyengine.config")
    fake_envipyengine_config.set = lambda *args, **kwargs: None
    fake_envipyengine.config = fake_envipyengine_config
    monkeypatch.setitem(sys.modules, "envipyengine", fake_envipyengine)
    monkeypatch.setitem(sys.modules, "envipyengine.config", fake_envipyengine_config)
    assert init_envi_engine("engine").tasks() is None
    monkeypatch.setattr(atmos_mod, "_execute_flaash_task", lambda *args, **kwargs: None)
    assert run_flaash_wrapper(({"OUTPUT_RASTER_URI": "out.tif"}, "params.txt", object())) == "out.tif"
    monkeypatch.setattr(atmos_mod, "_execute_flaash_task", lambda *args, **kwargs: None)
    result = run_flaash("in.tif", str(tmp_path / "out.tif"), params={"INPUT_RASTER": {"url": "in.tif", "factory": "URLRaster"}, "OUTPUT_RASTER_URI": str(tmp_path / "out.tif")}, envi_engine=object())
    assert result.output_raster.endswith("out.tif")
    monkeypatch.setattr(atmos_mod, "run_flaash_wrapper", lambda task: task[0]["OUTPUT_RASTER_URI"])
    class _Future:
        def __init__(self, value): self._value = value
        def result(self): return self._value
    class _Executor:
        def __init__(self, *args, **kwargs): pass
        def __enter__(self): return self
        def __exit__(self, *args): return False
        def submit(self, fn, task): return _Future(fn(task))
    monkeypatch.setattr(atmos_mod, "ProcessPoolExecutor", _Executor)
    monkeypatch.setattr(atmos_mod, "as_completed", lambda futures: futures)
    monkeypatch.setattr(atmos_mod, "tqdm", lambda it, **kwargs: it)
    assert parallel_flaash([({"OUTPUT_RASTER_URI": "a.tif"}, "p.txt")], object()) == ["a.tif"]
    monkeypatch.setattr(atmos_mod.Py6SCorrector, "run", lambda self, input_raster, output_raster, **kw: output_raster)
    assert atmospheric_correction("in.tif", "out.tif", method="py6s", solar_zenith=1, solar_azimuth=1, view_zenith=1, view_azimuth=1, day=1, month=1) == "out.tif"


def test_cloudmask_functions(monkeypatch, tmp_path: Path, make_test_raster) -> None:
    input_path = make_test_raster(tmp_path / "input.tif", count=3, data=np.arange(48, dtype=np.uint16).reshape(3, 4, 4))
    mask_path = tmp_path / "mask.tif"
    output_path = tmp_path / "masked.tif"
    monkeypatch.setitem(sys.modules, "omnicloudmask", SimpleNamespace(predict_from_array=lambda arr, **kwargs: np.array([[0, 1], [1, 0]]) if arr.shape[-1] == 2 else np.array([[0, 1, 0, 0], [0, 0, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0]])))
    create_cloud_mask_with_omnicloudmask(str(input_path), str(mask_path), 1, 2, 3, inference_resolution_m=None)
    apply_binary_cloud_mask_to_image(str(input_path), str(mask_path), str(output_path), output_nodata_value=999)
    result = cloudmask_raster(str(input_path), str(tmp_path / "masked2.tif"), str(tmp_path / "mask2.tif"), red_band_index=1, green_band_index=2, nir_band_index=3, inference_resolution_m=None)
    assert result.mask_pixel_count > 0


def test_fetch_functions(monkeypatch, tmp_path: Path) -> None:
    class _Response:
        headers = {"content-type": "image/tiff"}
        content = b"abc"
        def raise_for_status(self): return None
        def json(self): return {}
        text = ""
    monkeypatch.setattr(fetch_mod.requests, "get", lambda *args, **kwargs: _Response())
    monkeypatch.setenv("OPENTOPOGRAPHY_API_KEY", "key")
    assert download_opentopography_dem_for_bbox(min_lon=0, min_lat=0, max_lon=1, max_lat=1, output_tif_path=str(tmp_path / "dem.tif")).endswith("dem.tif")
    monkeypatch.setattr(fetch_mod, "_fetch_power_daily_point", lambda *args, **kwargs: {"aot550": 0.2, "water_vapor": 2.5, "ozone_cm_atm": 0.3})
    power = fetch_power_atmosphere_for_bbox(day_utc=__import__("datetime").date(2020, 1, 1), min_lon=0, min_lat=0, max_lon=1, max_lat=1)
    assert power.aot550 == 0.2
    fake_ee = SimpleNamespace(
        Authenticate=lambda: None,
        Initialize=lambda **kwargs: None,
    )
    monkeypatch.setitem(sys.modules, "ee", fake_ee)
    assert init_ee_client(env_file=None) is fake_ee
    monkeypatch.setattr(fetch_mod, "_fetch_collection_value", lambda *args, **kwargs: {"collection": "terra", "band_found": True, "raw_value": 2000.0, "abs_time_diff_hours": 1.0, "image_time_utc": "2020-01-01T00:00:00+00:00"})
    monkeypatch.setattr(fetch_mod, "_fetch_first_available_band_value", lambda *args, **kwargs: {"collection": "terra", "band_found": True, "raw_value": 100.0, "abs_time_diff_hours": 1.0, "image_time_utc": "2020-01-01T00:00:00+00:00"})
    modis = fetch_modis_water_vapor_for_bbox(scene_datetime_utc=__import__("datetime").datetime(2020, 1, 1, tzinfo=__import__("datetime").timezone.utc), min_lon=0, min_lat=0, max_lon=1, max_lat=1, ee=object())
    assert modis.status == "ok"


def test_pansharpen_and_radiometric(monkeypatch, tmp_path: Path) -> None:
    class _PanSharpen:
        def __init__(self, pan, mul): self.pan, self.mul = pan, mul
        def process(self, output_image_path, write_mask=False, overwrite=True):
            with rasterio.open(self.mul) as src:
                profile = src.profile.copy()
                data = src.read()
            with rasterio.open(output_image_path, "w", **profile) as dst:
                dst.write(data)
    fake_orthority = ModuleType("orthority")
    fake_orthority.PanSharpen = _PanSharpen
    monkeypatch.setitem(sys.modules, "orthority", fake_orthority)
    pansharpen_mod = importlib.import_module("vhrharmonize.preprocess.pansharpening")
    mul = tmp_path / "mul.tif"
    pan = tmp_path / "pan.tif"
    data = np.ones((1, 4, 4), dtype=np.int16)
    with rasterio.open(mul, "w", driver="GTiff", width=4, height=4, count=1, dtype="int16", transform=rasterio.transform.from_origin(0, 4, 1, 1), crs="EPSG:4326", nodata=-32768) as dst:
        dst.write(data)
    with rasterio.open(pan, "w", driver="GTiff", width=4, height=4, count=1, dtype="int16", transform=rasterio.transform.from_origin(0, 4, 1, 1), crs="EPSG:4326", nodata=-32768) as dst:
        dst.write(data)
    out = tmp_path / "ps.tif"
    pansharpen_mod.pansharpen_image(str(mul), str(pan), str(out), change_nodata_value=0)
    assert out.exists()
    monkeypatch.setattr(rad_mod, "_load_spectralmatch_pipeline", lambda: (lambda **kwargs: kwargs["shared_output_image_path"]))
    result = radiometric_normalization([str(out)], str(tmp_path / "norm.tif"), shared_cache="auto")
    assert result.endswith("norm.tif")


def test_worldview_named_radiometric_grouping(monkeypatch, tmp_path: Path) -> None:
    worldview = importlib.import_module("vhrharmonize.cli.worldview")
    root_output = str(tmp_path / "root.tif")
    child_output = str(tmp_path / "child.tif")
    spec = worldview._normalize_group_by_basename_spec({
        root_output: [
            "auto:*893242343*",
            {child_output: ["file:/external/ref.tif", "auto:*222222222*"]},
            "auto:*123456789*",
        ]
    })
    assert spec[root_output][1][child_output][0] == "file:/external/ref.tif"
    assert worldview._resolve_radiometric_group_output_path(
        output_name="$temp/from_temp.tif",
        temp_root=str(tmp_path / "tmp"),
    ) == str(tmp_path / "tmp" / "from_temp.tif")
    for bad_spec in ({"x.tif": "*missing-prefix*"}, {"": "auto:*"}, {"x.tif": []}):
        try:
            worldview._normalize_group_by_basename_spec(bad_spec)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid group spec was accepted")

    calls = []
    monkeypatch.setattr(
        worldview,
        "radiometric_normalization",
        lambda **kwargs: calls.append(kwargs) or kwargs["shared_output_image_path"],
    )
    args = SimpleNamespace(
        radiometric_normalization_kwargs_json=None,
        run_from_existing=False,
        run_from_existing_check_validity=False,
        log_to_console=False,
        keep_temp_dir=False,
        dtype="int16",
        radiometric_normalization_method="spectralmatch",
        calculate_overviews_radiometric_normalization=False,
    )
    available = [
        "/scene/auto_893242343.tif",
        "/scene/auto_222222222.tif",
        "/scene/auto_123456789.tif",
    ]
    result = worldview._run_named_radiometric_groups(
        spec,
        available_paths=available,
        args=args,
        temp_root=str(tmp_path / "tmp"),
    )
    assert result == root_output
    assert [call["shared_output_image_path"] for call in calls] == [
        child_output,
        root_output,
    ]
    assert calls[0]["shared_input_images"] == ["/external/ref.tif", "/scene/auto_222222222.tif"]
    assert calls[1]["shared_input_images"] == [
        "/scene/auto_893242343.tif",
        child_output,
        "/scene/auto_123456789.tif",
    ]
