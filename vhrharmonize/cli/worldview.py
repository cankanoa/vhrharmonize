#!/usr/bin/env python3
"""CLI wrapper for the WorldView preprocessing workflow."""

from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import uuid
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional

import geopandas as gpd

from vhrharmonize.io.geospatial import get_image_percentile_value, shp_to_gpkg
from vhrharmonize.logging_utils import log
from vhrharmonize.preprocess.atmospheric_correction import run_flaash, run_py6s
from vhrharmonize.preprocess.cloudmasking import cloudmask_raster
from vhrharmonize.preprocess.fetch_external_data import (
    fetch_modis_water_vapor_for_bbox,
    fetch_power_atmosphere_for_bbox,
)
from vhrharmonize.preprocess.orthorectification import (
    gcp_refined_rpc_orthorectification,
    resolve_output_resolution_for_crs,
)
from vhrharmonize.preprocess.pansharpening import pansharpen_image
from vhrharmonize.preprocess.radiometric_normalization import radiometric_normalization
from vhrharmonize.providers.worldview.files import (
    WorldViewImage,
    WorldViewScene,
    load_worldview_scenes_from_tif_files,
)
from vhrharmonize.pipelines.alignment import align_image_pair
from vhrharmonize.workflow_utils import (
    build_output_path_from_input,
    plan_step_outputs,
    resolve_output_dir,
    resolve_relative_to_input,
)

RASTER_STEP_ORDER = [
    "raw",
    "atmospheric_correction",
    "orthorectification",
    "pansharpen",
    "cloud_mask",
    "alignment",
    "radiometric_normalization",
]


@dataclass
class SceneWorkflowState:
    """Workflow state for a single discovered WorldView scene."""

    scene: WorldViewScene
    step_dirs: Dict[str, str]
    scene_bbox_wgs84: tuple[float, float, float, float]
    current_files: List[str]
    current_step: str = "raw"
    pan_ortho_path: Optional[str] = None
    fetch_atmosphere_result: Optional[Dict] = None
    py6s_effective_params: Optional[Dict] = None
    py6s_auto_atmos_estimate: Optional[Dict] = None
    alignment_result: Optional[object] = None
    cloud_mask_pixel_count: Optional[int] = None
    cloud_mask_path: Optional[str] = None


def _require_scene_image(scene: WorldViewScene, role: str) -> WorldViewImage:
    image = scene.get_image(role)
    if image is None:
        raise ValueError(f"WorldView scene is missing required {role} image: {scene.scene_id}_{scene.catalog_id}")
    return image


def _get_worldview_scene_step_path(scene: WorldViewScene, role: str, step_name: str) -> str:
    image = _require_scene_image(scene, role)
    if step_name == "raw":
        return image.tif_file
    step_path = image.step_file_paths.get(step_name)
    if not step_path:
        raise ValueError(
            f"WorldView scene is missing stored {role} output for step {step_name}: "
            f"{scene.scene_id}_{scene.catalog_id}"
        )
    return step_path


def _set_worldview_scene_step_path(scene: WorldViewScene, role: str, step_name: str, output_path: str) -> None:
    image = _require_scene_image(scene, role)
    image.step_file_paths[step_name] = output_path


def _parse_filter_basenames(raw_values: Optional[List[str]]) -> List[str]:
    if not raw_values:
        return []
    parsed: List[str] = []
    for raw in raw_values:
        for value in raw.split(","):
            value = value.strip()
            if value:
                parsed.append(value)
    return parsed


def _scene_bbox_wgs84_from_shp(shp_path: str) -> tuple[float, float, float, float]:
    gdf = gpd.read_file(shp_path)
    if gdf.empty:
        raise ValueError(f"Empty scene footprint shapefile: {shp_path}")
    if gdf.crs is None:
        raise ValueError(f"Scene footprint shapefile has no CRS: {shp_path}")
    gdf = gdf.to_crs(epsg=4326)
    minx, miny, maxx, maxy = gdf.total_bounds
    return float(minx), float(miny), float(maxx), float(maxy)


def _copy_atomic(src_path: str, dst_path: str, *, log_to_console: bool = False) -> None:
    tmp_output_path = f"{dst_path}.tmp-{uuid.uuid4().hex}"
    try:
        shutil.copy2(src_path, tmp_output_path)
    except PermissionError:
        shutil.copyfile(src_path, tmp_output_path)
        log("Copied file without metadata preservation", enabled=log_to_console, step="workflow")
    os.replace(tmp_output_path, dst_path)


def _parse_int_csv(raw_values: str) -> List[int]:
    parsed = []
    for value in raw_values.split(","):
        value = value.strip()
        if value:
            parsed.append(int(value))
    return parsed


def _parse_json_dict(raw_json: Optional[object]) -> Dict:
    if raw_json is None or raw_json == "":
        return {}
    if isinstance(raw_json, dict):
        return raw_json
    if not isinstance(raw_json, str):
        raise ValueError("Expected a JSON string or dictionary mapping.")
    parsed = json.loads(raw_json)
    if not isinstance(parsed, dict):
        raise ValueError("Expected a JSON object.")
    return parsed


def _collect_prefixed_kwargs(
    namespace: argparse.Namespace,
    prefix: str,
    *,
    transform_key=None,
) -> Dict:
    collected = {}
    for key, value in vars(namespace).items():
        if not key.startswith(prefix) or value is None:
            continue
        stripped_key = key[len(prefix):]
        if not stripped_key:
            continue
        if transform_key is not None:
            stripped_key = transform_key(stripped_key)
        collected[stripped_key] = value
    return collected


def _build_radiometric_kwargs(args: argparse.Namespace) -> Dict:
    radiometric_kwargs = _parse_json_dict(args.radiometric_normalization_kwargs_json)
    radiometric_kwargs.update(_collect_prefixed_kwargs(args, "match_"))
    radiometric_kwargs.setdefault("shared_debug_logs", args.log_to_console)
    return radiometric_kwargs


def _coerce_unknown_arg_value(raw_value: str):
    lowered = raw_value.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    if lowered == "null":
        return None
    try:
        if "." in raw_value:
            return float(raw_value)
        return int(raw_value)
    except ValueError:
        return raw_value


def _apply_unknown_prefixed_args(args: argparse.Namespace, unknown_args: List[str]) -> None:
    idx = 0
    while idx < len(unknown_args):
        token = unknown_args[idx]
        if not token.startswith("--"):
            raise SystemExit(f"Unsupported extra argument syntax: {token}")

        key_token = token[2:]
        inline_value = None
        if "=" in key_token:
            key_token, inline_value = key_token.split("=", 1)
        normalized_key = key_token.replace("-", "_")
        if not (normalized_key.startswith("match_") or normalized_key.startswith("flaash_param_")):
            raise SystemExit(f"Unrecognized argument: --{key_token}")

        if inline_value is not None:
            value = _coerce_unknown_arg_value(inline_value)
        elif idx + 1 < len(unknown_args) and not unknown_args[idx + 1].startswith("--"):
            idx += 1
            value = _coerce_unknown_arg_value(unknown_args[idx])
        else:
            value = True

        setattr(args, normalized_key, value)
        idx += 1


def _load_yaml_config(config_yaml_path: str) -> Dict:
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError(
            "PyYAML is required for --config-yaml. Install it with `pip install pyyaml`."
        ) from exc

    with open(config_yaml_path, "r", encoding="utf-8") as f:
        loaded = yaml.safe_load(f) or {}
    if not isinstance(loaded, dict):
        raise ValueError("Config YAML root must be a mapping/dictionary.")

    def _flatten_mapping(mapping: Dict, out: Dict) -> None:
        for key, value in mapping.items():
            normalized_key = str(key).replace("-", "_")
            if isinstance(value, dict):
                _flatten_mapping(value, out)
                continue
            if normalized_key in out:
                raise ValueError(
                    f"Duplicate config key after flattening nested sections: {normalized_key}"
                )
            out[normalized_key] = value

    normalized = {}
    _flatten_mapping(loaded, normalized)
    return normalized


def _normalize_config_defaults(config_defaults: Dict) -> Dict:
    normalized = dict(config_defaults)
    for list_key in ("input_file_glob", "input_dir", "filter_basename"):
        if list_key in normalized and isinstance(normalized[list_key], str):
            normalized[list_key] = [normalized[list_key]]
    if "input_dir" in normalized and "input_file_glob" not in normalized:
        normalized["input_file_glob"] = normalized.pop("input_dir")
    return normalized


def _resolve_fetch_atmosphere_source(args: argparse.Namespace) -> str:
    if args.fetch_atmosphere_source != "auto":
        return args.fetch_atmosphere_source
    if args.atmospheric_method == "flaash":
        return "modis_gee"
    return "nasa_power"


def _write_json(path: str, payload: Dict) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def _read_json(path: str) -> Dict:
    with open(path, "r", encoding="utf-8") as handle:
        loaded = json.load(handle)
    if not isinstance(loaded, dict):
        raise ValueError(f"Expected JSON object at {path}")
    return loaded


def _collect_input_tif_files(input_file_globs: List[str]) -> List[str]:
    tif_files: List[str] = []
    for pattern in input_file_globs:
        matches = glob.glob(pattern, recursive=True)
        tif_files.extend(
            path
            for path in matches
            if os.path.isfile(path) and os.path.splitext(path)[1].lower() == ".tif"
        )
    return sorted({os.path.abspath(path) for path in tif_files})


def _resolve_scene_step_dirs(args: argparse.Namespace, scene: WorldViewScene) -> Dict[str, str]:
    input_folder = scene.root_folder_path
    scene_key = scene.primary_basename or f"{scene.scene_id}_{scene.catalog_id}"
    return {
        "temp_root": resolve_relative_to_input(args.temp_dir, input_folder),
        "scene_work": resolve_output_dir(None, input_folder=input_folder, temp_dir=args.temp_dir, step_name=os.path.join("_scene", scene_key)),
        "fetch_atmosphere": resolve_output_dir(args.fetch_atmosphere_output_dir, input_folder=input_folder, temp_dir=args.temp_dir, step_name="fetch_atmosphere"),
        "atmospheric_correction": resolve_output_dir(args.atmospheric_correction_output_dir, input_folder=input_folder, temp_dir=args.temp_dir, step_name="atmospheric_correction"),
        "orthorectification": resolve_output_dir(args.orthorectification_output_dir, input_folder=input_folder, temp_dir=args.temp_dir, step_name="orthorectification"),
        "pansharpen": resolve_output_dir(args.pansharpen_output_dir, input_folder=input_folder, temp_dir=args.temp_dir, step_name="pansharpen"),
        "cloud_mask": resolve_output_dir(args.cloud_mask_output_dir, input_folder=input_folder, temp_dir=args.temp_dir, step_name="cloud_mask"),
        "cloud_mask_mask": resolve_output_dir(args.cloud_mask_mask_output_dir, input_folder=input_folder, temp_dir=args.temp_dir, step_name="cloud_mask_mask"),
        "alignment": resolve_output_dir(args.alignment_output_dir, input_folder=input_folder, temp_dir=args.temp_dir, step_name="alignment"),
        "radiometric_normalization": resolve_output_dir(args.radiometric_normalization_output_dir, input_folder=input_folder, temp_dir=args.temp_dir, step_name="radiometric_normalization"),
        "final": resolve_output_dir(args.output_dir, input_folder=input_folder, temp_dir=args.temp_dir, step_name="final"),
    }


def _get_atmospheric_extension(args: argparse.Namespace) -> str:
    return ".dat" if args.atmospheric_method == "flaash" else ".tif"


def _get_expected_mul_output_path(state: SceneWorkflowState, args: argparse.Namespace, step_name: str) -> str:
    mul_image = _require_scene_image(state.scene, "mul")
    if step_name == "raw":
        return mul_image.tif_file
    if step_name == "atmospheric_correction":
        return build_output_path_from_input(
            mul_image.tif_file,
            state.step_dirs["atmospheric_correction"],
            suffix=args.atmospheric_correction_output_suffix,
            extension=_get_atmospheric_extension(args),
        )
    if step_name == "orthorectification":
        return build_output_path_from_input(
            mul_image.tif_file,
            state.step_dirs["orthorectification"],
            suffix=args.orthorectification_output_suffix,
        )
    if step_name == "pansharpen":
        return build_output_path_from_input(
            mul_image.tif_file,
            state.step_dirs["pansharpen"],
            suffix=args.pansharpen_output_suffix,
        )
    if step_name == "cloud_mask":
        return build_output_path_from_input(
            mul_image.tif_file,
            state.step_dirs["cloud_mask"],
            suffix=args.cloud_mask_output_suffix,
        )
    if step_name == "alignment":
        return build_output_path_from_input(
            mul_image.tif_file,
            state.step_dirs["alignment"],
            suffix=args.alignment_output_suffix,
        )
    if step_name == "radiometric_normalization":
        return build_output_path_from_input(
            mul_image.tif_file,
            state.step_dirs["radiometric_normalization"],
            suffix=args.radiometric_normalization_output_suffix,
        )
    raise ValueError(f"Unsupported step name: {step_name}")


def _get_expected_pan_ortho_path(state: SceneWorkflowState, args: argparse.Namespace) -> Optional[str]:
    pan_image = state.scene.pan_image
    if pan_image is None:
        return None
    return build_output_path_from_input(
        pan_image.tif_file,
        state.step_dirs["orthorectification"],
        suffix=args.orthorectification_pan_output_suffix,
    )


def _set_scene_current_step(state: SceneWorkflowState, args: argparse.Namespace, step_name: str) -> None:
    current_path = _get_expected_mul_output_path(state, args, step_name)
    if not os.path.isfile(current_path):
        raise ValueError(f"Expected existing output for last_run_step={step_name}: {current_path}")
    _set_worldview_scene_step_path(state.scene, "mul", step_name, current_path)
    state.current_files = [_get_worldview_scene_step_path(state.scene, "mul", step_name)]
    state.current_step = step_name
    state.scene.step_outputs[step_name] = [state.current_files[0]]
    if RASTER_STEP_ORDER.index(step_name) >= RASTER_STEP_ORDER.index("orthorectification"):
        state.pan_ortho_path = args.existing_pan_ortho_input or _get_expected_pan_ortho_path(state, args)
        if args.run_pansharpen and state.pan_ortho_path and not os.path.isfile(state.pan_ortho_path):
            raise ValueError(f"Expected existing panchromatic ortho output: {state.pan_ortho_path}")
        if state.pan_ortho_path:
            _set_worldview_scene_step_path(state.scene, "pan", "orthorectification", state.pan_ortho_path)


def _initialize_scene_state(scene: WorldViewScene, args: argparse.Namespace) -> SceneWorkflowState:
    mul_image = _require_scene_image(scene, "mul")
    pan_image = _require_scene_image(scene, "pan")
    if mul_image.standardized_metadata is None:
        raise ValueError(f"WorldView scene is missing multispectral metadata: {scene.scene_id}_{scene.catalog_id}")
    if pan_image.standardized_metadata is None:
        raise ValueError(f"WorldView scene is missing panchromatic metadata: {scene.scene_id}_{scene.catalog_id}")
    if mul_image.shp_file is None:
        raise ValueError(f"WorldView scene is missing multispectral shapefile: {scene.primary_basename}")

    state = SceneWorkflowState(
        scene=scene,
        step_dirs=_resolve_scene_step_dirs(args, scene),
        scene_bbox_wgs84=_scene_bbox_wgs84_from_shp(mul_image.shp_file),
        current_files=[_get_worldview_scene_step_path(scene, "mul", "raw")],
    )
    if args.last_run_step != "raw":
        _set_scene_current_step(state, args, args.last_run_step)
    return state


def _register_step_outputs(
    state: SceneWorkflowState,
    step_name: str,
    output_paths: List[str],
    *,
    image_role: str = "mul",
) -> List[str]:
    state.scene.step_outputs[step_name] = list(output_paths)
    if output_paths:
        _set_worldview_scene_step_path(state.scene, image_role, step_name, output_paths[0])
    return output_paths


def _run_cloud_mask_command(
    command_template: str,
    input_image_path: str,
    output_image_path: str,
    scene_root_path: str,
    image_basename: str,
    *,
    log_to_console: bool = False,
) -> None:
    command = command_template.format(
        input=input_image_path,
        output=output_image_path,
        scene_root=scene_root_path,
        image_basename=image_basename,
    )
    log("Running external cloud mask command", enabled=log_to_console, step="cloud_mask")
    subprocess.run(command, shell=True, check=True)


def _run_fetch_atmosphere_step(state: SceneWorkflowState, args: argparse.Namespace) -> None:
    if not args.run_fetch_atmosphere:
        return
    mul_image = state.scene.mul_image
    if mul_image is None or mul_image.standardized_metadata is None:
        return
    plan = plan_step_outputs(
        [mul_image.tif_file],
        output_dir=state.step_dirs["fetch_atmosphere"],
        suffix=args.fetch_atmosphere_output_suffix,
        extension=".json",
        skip_existing=args.skip_existing,
    )
    if not plan.pending_output_paths:
        state.fetch_atmosphere_result = _read_json(plan.output_paths[0])
        _register_step_outputs(state, "fetch_atmosphere", plan.output_paths)
        return

    fetch_source = _resolve_fetch_atmosphere_source(args)
    if fetch_source == "nasa_power":
        estimate = fetch_power_atmosphere_for_bbox(
            day_utc=mul_image.standardized_metadata.resolve_scene_datetime().date(),
            min_lon=state.scene_bbox_wgs84[0],
            min_lat=state.scene_bbox_wgs84[1],
            max_lon=state.scene_bbox_wgs84[2],
            max_lat=state.scene_bbox_wgs84[3],
            grid_size=args.fetch_atmosphere_grid_size,
            search_days=args.fetch_atmosphere_search_days,
            timeout_s=args.fetch_atmosphere_timeout_s,
            endpoint=args.fetch_atmosphere_power_endpoint,
            log_to_console=args.log_to_console,
        )
        result = {
            "source": estimate.source,
            "date_used": estimate.date_used,
            "sample_count": estimate.sample_count,
            "aot550": estimate.aot550,
            "water_vapor": estimate.water_vapor,
            "ozone_cm_atm": estimate.ozone_cm_atm,
        }
    elif fetch_source == "modis_gee":
        estimate = fetch_modis_water_vapor_for_bbox(
            scene_datetime_utc=mul_image.standardized_metadata.resolve_scene_datetime(),
            min_lon=state.scene_bbox_wgs84[0],
            min_lat=state.scene_bbox_wgs84[1],
            max_lon=state.scene_bbox_wgs84[2],
            max_lat=state.scene_bbox_wgs84[3],
            ee_project=args.fetch_atmosphere_ee_project,
            authenticate=args.fetch_atmosphere_authenticate,
            env_file=args.fetch_atmosphere_env_file,
            hours_window=args.fetch_atmosphere_hours_window,
            log_to_console=args.log_to_console,
        )
        result = estimate.to_dict()
    else:
        raise ValueError(f"Unsupported fetch atmosphere source: {fetch_source}")

    _write_json(plan.pending_output_paths[0], result)
    state.fetch_atmosphere_result = result
    _register_step_outputs(state, "fetch_atmosphere", plan.output_paths)


def _run_atmospheric_correction_step(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    if not args.run_atmospheric_correction:
        return state.current_files
    mul_image = state.scene.mul_image
    if mul_image is None or mul_image.standardized_metadata is None or mul_image.shp_file is None:
        raise ValueError("WorldView scene is missing multispectral inputs for atmospheric correction.")

    plan = plan_step_outputs(
        state.current_files,
        output_dir=state.step_dirs["atmospheric_correction"],
        suffix=args.atmospheric_correction_output_suffix,
        extension=_get_atmospheric_extension(args),
        skip_existing=args.skip_existing,
    )
    if not plan.pending_output_paths:
        state.current_files = _register_step_outputs(state, "atmospheric_correction", plan.output_paths)
        state.current_step = "atmospheric_correction"
        return state.current_files

    input_raster = plan.pending_input_paths[0]
    output_raster = plan.pending_output_paths[0]

    if args.atmospheric_method == "flaash":
        if args.skip_flaash:
            if not args.existing_flaash_input:
                raise ValueError("--existing-flaash-input is required when --skip-flaash is set.")
            state.current_files = _register_step_outputs(state, "atmospheric_correction", [args.existing_flaash_input])
            state.current_step = "atmospheric_correction"
            return state.current_files

        mul_gpkg_path = os.path.join(state.step_dirs["scene_work"], f"{mul_image.basename}.gpkg")
        shp_to_gpkg(mul_image.shp_file, mul_gpkg_path, args.footprint_epsg)
        custom_flaash_params = _collect_prefixed_kwargs(
            args,
            "flaash_param_",
            transform_key=lambda value: value.upper(),
        )
        modtran_atm = args.flaash_modtran_atm or "Mid-Latitude Summer"
        modtran_aer = args.flaash_modtran_aer or "Maritime"
        use_aerosol = args.flaash_use_aerosol or "Disabled"
        default_visibility = args.flaash_default_visibility
        if state.fetch_atmosphere_result:
            modtran_atm = state.fetch_atmosphere_result.get("modtran_atm") or modtran_atm
            if state.fetch_atmosphere_result.get("default_visibility") is not None:
                default_visibility = float(state.fetch_atmosphere_result["default_visibility"])
            if state.fetch_atmosphere_result.get("water_vapor_preset") is not None:
                custom_flaash_params["WATER_VAPOR_PRESET"] = float(
                    state.fetch_atmosphere_result["water_vapor_preset"]
                )

        run_flaash(
            input_raster=input_raster,
            output_raster=output_raster,
            metadata=mul_image.standardized_metadata,
            dem_file_path=args.dem_file_path,
            footprint_vector_path=mul_gpkg_path,
            envi_engine_path=args.envi_engine_path,
            convert_paths_for_windows=True,
            output_params_path=f"{output_raster}.params.txt",
            dem_ground_percentile=args.flaash_dem_ground_percentile,
            modtran_atm=modtran_atm,
            modtran_aer=modtran_aer,
            use_aerosol=use_aerosol,
            default_visibility=default_visibility,
            custom_params=custom_flaash_params or None,
            log_to_console=args.log_to_console,
        )
    elif args.atmospheric_method == "py6s":
        mul_gpkg_path = os.path.join(state.step_dirs["scene_work"], f"{mul_image.basename}.gpkg")
        shp_to_gpkg(mul_image.shp_file, mul_gpkg_path, args.footprint_epsg)
        ground_elevation_m = get_image_percentile_value(
            args.dem_file_path,
            percentile=args.flaash_dem_ground_percentile,
            mask=mul_gpkg_path,
        )
        py6s_result = run_py6s(
            input_raster=input_raster,
            output_raster=output_raster,
            metadata=mul_image.standardized_metadata,
            ground_elevation_km=ground_elevation_m / 1000.0,
            atmosphere_profile=args.py6s_atmosphere_profile,
            aerosol_profile=args.py6s_aerosol_profile,
            aot550=float(state.fetch_atmosphere_result["aot550"]) if state.fetch_atmosphere_result and state.fetch_atmosphere_result.get("aot550") is not None else args.py6s_aot550,
            visibility_km=args.py6s_visibility,
            water_vapor=float(state.fetch_atmosphere_result["water_vapor"]) if state.fetch_atmosphere_result and state.fetch_atmosphere_result.get("water_vapor") is not None else args.py6s_water_vapor,
            ozone=float(state.fetch_atmosphere_result["ozone_cm_atm"]) if state.fetch_atmosphere_result and state.fetch_atmosphere_result.get("ozone_cm_atm") is not None else args.py6s_ozone,
            sixs_executable=args.py6s_executable,
            output_scale_factor=args.py6s_output_scale_factor,
            output_dtype=args.py6s_output_dtype,
            use_imd_radiance_calibration=args.py6s_use_imd_radiance_calibration,
            use_worldview_gain_offset_adjustment=args.py6s_use_worldview_gain_offset_adjustment,
            auto_atmos_source="none",
            log_to_console=args.log_to_console,
        )
        state.py6s_effective_params = py6s_result.effective_params
        state.py6s_auto_atmos_estimate = py6s_result.auto_atmos_estimate
    else:
        state.current_files = [input_raster]
        state.current_step = "raw"
        return state.current_files

    state.current_files = _register_step_outputs(state, "atmospheric_correction", plan.output_paths)
    state.current_step = "atmospheric_correction"
    return state.current_files


def _run_orthorectification_step(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    if not args.run_orthorectification:
        return state.current_files
    mul_image = state.scene.mul_image
    pan_image = state.scene.pan_image
    if mul_image is None or pan_image is None:
        raise ValueError("WorldView scene is missing multispectral or panchromatic image.")

    if args.existing_mul_ortho_input:
        state.current_files = _register_step_outputs(state, "orthorectification", [args.existing_mul_ortho_input])
    else:
        plan = plan_step_outputs(
            state.current_files,
            output_dir=state.step_dirs["orthorectification"],
            suffix=args.orthorectification_output_suffix,
            skip_existing=args.skip_existing,
        )
        for input_path, output_path in zip(plan.pending_input_paths, plan.pending_output_paths):
            gcp_refined_rpc_orthorectification(
                input_path,
                output_path,
                args.dem_file_path,
                args.epsg,
                output_nodata_value=args.nodata_value,
                dtype=args.dtype,
                output_resolution=resolve_output_resolution_for_crs(
                    args.epsg,
                    mul_image.standardized_metadata.product_resolution,
                ),
                log_to_console=args.log_to_console,
            )
        state.current_files = _register_step_outputs(state, "orthorectification", plan.output_paths)

    if args.run_pansharpen:
        if args.existing_pan_ortho_input:
            state.pan_ortho_path = args.existing_pan_ortho_input
        else:
            pan_plan = plan_step_outputs(
                [pan_image.tif_file],
                output_dir=state.step_dirs["orthorectification"],
                suffix=args.orthorectification_pan_output_suffix,
                skip_existing=args.skip_existing,
            )
            for input_path, output_path in zip(pan_plan.pending_input_paths, pan_plan.pending_output_paths):
                gcp_refined_rpc_orthorectification(
                    input_path,
                    output_path,
                    args.dem_file_path,
                    args.epsg,
                    output_nodata_value=args.nodata_value,
                    dtype=args.dtype,
                    output_resolution=resolve_output_resolution_for_crs(
                        args.epsg,
                        pan_image.standardized_metadata.product_resolution,
                    ),
                    log_to_console=args.log_to_console,
                )
            state.pan_ortho_path = pan_plan.output_paths[0]
            _register_step_outputs(state, "orthorectification_pan", pan_plan.output_paths, image_role="pan")

    state.current_step = "orthorectification"
    return state.current_files


def _run_pansharpen_step(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    if not args.run_pansharpen:
        return state.current_files
    if state.pan_ortho_path is None:
        raise ValueError("Panchromatic orthorectified path is required for pansharpening.")
    plan = plan_step_outputs(
        state.current_files,
        output_dir=state.step_dirs["pansharpen"],
        suffix=args.pansharpen_output_suffix,
        skip_existing=args.skip_existing,
    )
    for input_path, output_path in zip(plan.pending_input_paths, plan.pending_output_paths):
        pansharpen_image(
            input_path,
            state.pan_ortho_path,
            output_path,
            change_nodata_value=args.nodata_value,
            log_to_console=args.log_to_console,
        )
    state.current_files = _register_step_outputs(state, "pansharpen", plan.output_paths)
    state.current_step = "pansharpen"
    return state.current_files


def _run_cloud_mask_step(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    if not args.run_cloud_mask:
        return state.current_files
    mul_image = state.scene.mul_image
    if mul_image is None:
        return state.current_files

    output_plan = plan_step_outputs(
        state.current_files,
        output_dir=state.step_dirs["cloud_mask"],
        suffix=args.cloud_mask_output_suffix,
        skip_existing=args.skip_existing,
    )
    mask_plan = plan_step_outputs(
        state.current_files,
        output_dir=state.step_dirs["cloud_mask_mask"],
        suffix=args.cloud_mask_mask_suffix,
        skip_existing=args.skip_existing,
    )

    if args.cloud_mask_command:
        for input_path, output_path in zip(output_plan.pending_input_paths, output_plan.pending_output_paths):
            _run_cloud_mask_command(
                args.cloud_mask_command,
                input_path,
                output_path,
                state.scene.root_folder_path,
                mul_image.basename,
                log_to_console=args.log_to_console,
            )
    else:
        cloud_classes = _parse_int_csv(args.cloud_mask_classes)
        omnicloud_kwargs = _parse_json_dict(args.cloud_mask_omnicloud_kwargs_json)
        pending_pairs = zip(
            output_plan.pending_input_paths,
            output_plan.pending_output_paths,
            mask_plan.pending_output_paths,
        )
        for input_path, output_path, mask_output_path in pending_pairs:
            cloudmask_result = cloudmask_raster(
                input_image_path=input_path,
                output_raster_path=output_path,
                output_mask_path=mask_output_path,
                red_band_index=args.cloud_mask_red_band_index,
                green_band_index=args.cloud_mask_green_band_index,
                nir_band_index=args.cloud_mask_nir_band_index,
                cloud_classes=cloud_classes,
                buffer_pixels=args.cloud_buffer_pixels,
                omnicloud_kwargs=omnicloud_kwargs,
                inference_resolution_m=args.cloud_mask_inference_resolution_m,
                output_nodata_value=args.nodata_value,
                allow_mask_reprojection=True,
                log_to_console=args.log_to_console,
            )
            state.cloud_mask_pixel_count = cloudmask_result.mask_pixel_count
            state.cloud_mask_path = cloudmask_result.output_mask_path

    state.current_files = _register_step_outputs(state, "cloud_mask", output_plan.output_paths)
    if mask_plan.output_paths:
        _register_step_outputs(state, "cloud_mask_mask", mask_plan.output_paths)
        state.cloud_mask_path = mask_plan.output_paths[0]
    state.current_step = "cloud_mask"
    return state.current_files


def _run_alignment_step(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    if not args.run_alignment:
        return state.current_files
    plan = plan_step_outputs(
        state.current_files,
        output_dir=state.step_dirs["alignment"],
        suffix=args.alignment_output_suffix,
        skip_existing=args.skip_existing,
    )
    resolved_alignment_temp_dir = (
        resolve_relative_to_input(args.alignment_temp_dir, state.scene.root_folder_path)
        if args.alignment_temp_dir
        else state.step_dirs["temp_root"]
    )
    for input_path, output_path in zip(plan.pending_input_paths, plan.pending_output_paths):
        state.alignment_result = align_image_pair(
            moving_image_path=input_path,
            fixed_image_path=args.alignment_fixed_image,
            output_image_path=output_path,
            band_index=0,
            moving_band_index=args.alignment_moving_band_index,
            fixed_band_index=args.alignment_fixed_band_index,
            tiling=not args.alignment_no_tiling,
            tile_size=args.alignment_tile_size,
            tile_buffer=args.alignment_tile_buffer,
            parameter_map=args.alignment_parameter_map,
            moving_nodata=args.alignment_moving_nodata,
            fixed_nodata=args.alignment_fixed_nodata,
            output_nodata=args.alignment_output_nodata,
            min_valid_fraction=args.alignment_min_valid_fraction,
            temp_dir=resolved_alignment_temp_dir,
            keep_temp_dir=args.alignment_keep_temp_dir,
            log_to_console=args.log_to_console or args.alignment_log_to_console,
            clip_fixed_to_moving=args.alignment_clip_fixed_to_moving,
            output_on_moving_grid=args.alignment_output_on_moving_grid,
            enforce_mutual_valid_mask=args.alignment_enforce_mutual_valid_mask,
            registration_mode=args.alignment_registration_mode,
        )
    state.current_files = _register_step_outputs(state, "alignment", plan.output_paths)
    state.current_step = "alignment"
    return state.current_files


def _run_radiometric_normalization_step(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    if not args.run_radiometric_normalization:
        return state.current_files
    plan = plan_step_outputs(
        state.current_files,
        output_dir=state.step_dirs["radiometric_normalization"],
        suffix=args.radiometric_normalization_output_suffix,
        skip_existing=args.skip_existing,
    )
    radiometric_kwargs = _build_radiometric_kwargs(args)
    radiometric_kwargs.setdefault("shared_temp_dir", state.step_dirs["temp_root"])
    radiometric_kwargs.setdefault("shared_custom_nodata_value", args.nodata_value)
    radiometric_kwargs.setdefault("shared_output_dtype", args.dtype)
    for input_path, output_path in zip(plan.pending_input_paths, plan.pending_output_paths):
        radiometric_normalization(
            shared_input_images=[input_path],
            shared_output_image_path=output_path,
            method=args.radiometric_normalization_method,
            log_to_console=args.log_to_console,
            **radiometric_kwargs,
        )
    state.current_files = _register_step_outputs(state, "radiometric_normalization", plan.output_paths)
    state.current_step = "radiometric_normalization"
    return state.current_files


def _final_output_paths(state: SceneWorkflowState, args: argparse.Namespace) -> tuple[str, str]:
    mul_image = state.scene.mul_image
    if mul_image is None:
        raise ValueError("WorldView scene is missing multispectral image.")
    final_image_path = build_output_path_from_input(
        mul_image.tif_file,
        state.step_dirs["final"],
        suffix=args.output_suffix,
    )
    final_metadata_path = os.path.join(
        state.step_dirs["final"],
        f"{mul_image.basename}{args.output_suffix}_metadata.json",
    )
    return final_image_path, final_metadata_path


def _scene_final_outputs_complete(state: SceneWorkflowState, args: argparse.Namespace) -> bool:
    final_image_path, final_metadata_path = _final_output_paths(state, args)
    return os.path.isfile(final_image_path) and os.path.isfile(final_metadata_path)


def _write_scene_report(state: SceneWorkflowState, args: argparse.Namespace, *, scene_started_utc: str) -> None:
    mul_image = state.scene.mul_image
    pan_image = state.scene.pan_image
    if mul_image is None or pan_image is None:
        raise ValueError("WorldView scene is missing required images for metadata reporting.")
    final_scene_path, scene_metadata_path = _final_output_paths(state, args)
    payload = {
        "scene": {
            "scene_id": state.scene.scene_id,
            "catalog_id": state.scene.catalog_id,
            "scene_root": state.scene.root_folder_path,
            "mul_photo_basename": mul_image.basename,
            "pan_photo_basename": pan_image.basename,
            "started_utc": scene_started_utc,
            "completed_utc": datetime.utcnow().isoformat() + "Z",
        },
        "inputs": {
            "mul_imd_file": mul_image.imd_file,
            "mul_tif_file": mul_image.tif_file,
            "mul_shp_file": mul_image.shp_file,
            "pan_imd_file": pan_image.imd_file,
            "pan_tif_file": pan_image.tif_file,
            "pan_shp_file": pan_image.shp_file,
            "dem_file_path": args.dem_file_path,
        },
        "standardized_metadata": {
            "mul": mul_image.standardized_metadata.to_dict() if mul_image.standardized_metadata else None,
            "pan": pan_image.standardized_metadata.to_dict() if pan_image.standardized_metadata else None,
        },
        "workflow": {
            "last_run_step": args.last_run_step,
            "run_fetch_atmosphere": args.run_fetch_atmosphere,
            "run_atmospheric_correction": args.run_atmospheric_correction,
            "run_orthorectification": args.run_orthorectification,
            "run_pansharpen": args.run_pansharpen,
            "run_cloud_mask": args.run_cloud_mask,
            "run_alignment": args.run_alignment,
            "run_radiometric_normalization": args.run_radiometric_normalization,
            "atmospheric_method": args.atmospheric_method,
            "epsg": args.epsg,
            "nodata_value": args.nodata_value,
            "dtype": args.dtype,
            "temp_dir": state.step_dirs["temp_root"],
        },
        "fetch_atmosphere": state.fetch_atmosphere_result,
        "py6s": {
            "effective": state.py6s_effective_params,
            "auto_atmos_estimate": state.py6s_auto_atmos_estimate,
        },
        "cloud_mask": {
            "mask_output_path": state.cloud_mask_path,
            "mask_pixel_count": state.cloud_mask_pixel_count,
        },
        "outputs": {
            "step_dirs": state.step_dirs,
            "step_outputs": state.scene.step_outputs,
            "final_scene_path": final_scene_path,
            "scene_metadata_path": scene_metadata_path,
        },
        "alignment": (
            {
                "fixed_image": args.alignment_fixed_image,
                "result": (
                    {
                        "output_image_path": state.alignment_result.output_image_path,
                        "total_tiles": state.alignment_result.total_tiles,
                        "successful_tiles": state.alignment_result.successful_tiles,
                        "skipped_tiles": state.alignment_result.skipped_tiles,
                        "temp_dir": state.alignment_result.temp_dir,
                    }
                    if state.alignment_result
                    else None
                ),
            }
        ),
    }
    _write_json(scene_metadata_path, payload)
    state.scene.metadata_report_path = scene_metadata_path


def run_workflow(args: argparse.Namespace) -> int:
    """Run full-scene WorldView preprocessing using input tif globs."""
    filter_basenames = _parse_filter_basenames(args.filter_basename)
    tif_files = _collect_input_tif_files(args.input_file_glob)
    if not tif_files:
        raise ValueError("No tif files matched --input-file-glob.")

    scenes = load_worldview_scenes_from_tif_files(tif_files, filter_basenames=filter_basenames)
    for scene in scenes:
        state = _initialize_scene_state(scene, args)
        log(
            f"Processing {scene.primary_basename or f'{scene.scene_id}_{scene.catalog_id}'}",
            enabled=args.log_to_console,
            step="workflow",
        )

        if args.skip_existing and _scene_final_outputs_complete(state, args):
            log("Skipping existing final output", enabled=args.log_to_console, step="workflow")
            continue

        scene_started_utc = datetime.utcnow().isoformat() + "Z"
        _run_fetch_atmosphere_step(state, args)

        if RASTER_STEP_ORDER.index(state.current_step) < RASTER_STEP_ORDER.index("atmospheric_correction"):
            _run_atmospheric_correction_step(state, args)
        if RASTER_STEP_ORDER.index(state.current_step) < RASTER_STEP_ORDER.index("orthorectification"):
            _run_orthorectification_step(state, args)
        if RASTER_STEP_ORDER.index(state.current_step) < RASTER_STEP_ORDER.index("pansharpen"):
            _run_pansharpen_step(state, args)
        if RASTER_STEP_ORDER.index(state.current_step) < RASTER_STEP_ORDER.index("cloud_mask"):
            _run_cloud_mask_step(state, args)
        if RASTER_STEP_ORDER.index(state.current_step) < RASTER_STEP_ORDER.index("alignment"):
            _run_alignment_step(state, args)
        if RASTER_STEP_ORDER.index(state.current_step) < RASTER_STEP_ORDER.index("radiometric_normalization"):
            _run_radiometric_normalization_step(state, args)

        final_image_path, _ = _final_output_paths(state, args)
        final_plan = plan_step_outputs(
            state.current_files,
            output_dir=state.step_dirs["final"],
            suffix=args.output_suffix,
            skip_existing=args.skip_existing,
        )
        for input_path, output_path in zip(final_plan.pending_input_paths, final_plan.pending_output_paths):
            _copy_atomic(input_path, output_path, log_to_console=args.log_to_console)
        state.current_files = final_plan.output_paths
        _register_step_outputs(state, "final", final_plan.output_paths)
        log(f"Wrote final output {final_image_path}", enabled=args.log_to_console, step="workflow")
        _write_scene_report(state, args, scene_started_utc=scene_started_utc)
        log("Scene complete", enabled=args.log_to_console, step="workflow")

    log("All processing complete", enabled=args.log_to_console, step="workflow")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the WorldView preprocessing CLI."""
    parser = argparse.ArgumentParser(description="Run WorldView preprocessing on discovered tif scenes.")
    parser.add_argument("--config-yaml", help="Optional YAML config file.")
    parser.add_argument(
        "--input-file-glob",
        action="append",
        help="Recursive glob used to find input tif files, for example '/data/**/*.tif'.",
    )
    parser.add_argument("--input-dir", dest="input_file_glob", action="append", help=argparse.SUPPRESS)
    parser.add_argument("--dem-file-path", help="DEM GeoTIFF path in WGS84 ellipsoidal height.")
    parser.add_argument("--envi-engine-path", help="Path to ENVI taskengine executable.")
    parser.add_argument("--atmospheric-method", choices=["flaash", "py6s", "none"], default="py6s")
    parser.add_argument("--epsg", type=int, default=4326)
    parser.add_argument("--nodata-value", type=float, default=-9999)
    parser.add_argument("--dtype", default="int16")
    parser.add_argument("--log-to-console", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument(
        "--last-run-step",
        choices=RASTER_STEP_ORDER,
        default="raw",
        help="Most recent raster step already completed for the current scene outputs.",
    )
    parser.add_argument("--flaash-dem-ground-percentile", type=float, default=50.0)
    parser.add_argument("--flaash-modtran-atm")
    parser.add_argument("--flaash-modtran-aer")
    parser.add_argument("--flaash-use-aerosol")
    parser.add_argument("--flaash-default-visibility", type=float)
    parser.add_argument("--py6s-atmosphere-profile", default="midlatitude_summer")
    parser.add_argument("--py6s-aerosol-profile", default="maritime")
    parser.add_argument("--py6s-aot550", type=float, default=0.2)
    parser.add_argument("--py6s-visibility", type=float)
    parser.add_argument("--py6s-water-vapor", type=float, default=2.5)
    parser.add_argument("--py6s-ozone", type=float, default=0.3)
    parser.add_argument("--py6s-output-scale-factor", type=float, default=10000.0)
    parser.add_argument("--py6s-output-dtype", default="int16")
    parser.add_argument("--py6s-executable")
    parser.add_argument("--py6s-use-imd-radiance-calibration", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--py6s-use-worldview-gain-offset-adjustment", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--py6s-auto-atmos-source", choices=["none", "nasa_power"], default="none")
    parser.add_argument("--py6s-auto-atmos-grid-size", type=int, default=3)
    parser.add_argument("--py6s-auto-atmos-search-days", type=int, default=1)
    parser.add_argument("--py6s-auto-atmos-timeout-s", type=float, default=30.0)
    parser.add_argument("--py6s-auto-atmos-power-endpoint", default="https://power.larc.nasa.gov/api/temporal/daily/point")
    parser.add_argument("--footprint-epsg", type=int, default=4326)
    parser.add_argument("--filter-basename", action="append")
    parser.add_argument("--output-dir")
    parser.add_argument("--output-suffix", default="_final")
    parser.add_argument("--fetch-atmosphere-output-dir")
    parser.add_argument("--fetch-atmosphere-output-suffix", default="_atmosphere")
    parser.add_argument("--atmospheric-correction-output-dir")
    parser.add_argument("--atmospheric-correction-output-suffix", default="_atmospheric")
    parser.add_argument("--radiometric-normalization-output-dir")
    parser.add_argument("--radiometric-normalization-output-suffix", default="_normalized")
    parser.add_argument("--orthorectification-output-dir")
    parser.add_argument("--orthorectification-output-suffix", default="_ortho")
    parser.add_argument("--orthorectification-pan-output-suffix", default="_pan_ortho")
    parser.add_argument("--pansharpen-output-dir")
    parser.add_argument("--pansharpen-output-suffix", default="_pansharpen")
    parser.add_argument("--cloud-mask-output-dir")
    parser.add_argument("--cloud-mask-mask-output-dir")
    parser.add_argument("--alignment-output-dir")
    parser.add_argument("--skip-existing", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--temp-dir", default="/tmp")
    parser.add_argument("--scratch-dir", dest="temp_dir", help=argparse.SUPPRESS)
    parser.add_argument("--run-fetch-atmosphere", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--run-atmospheric-correction", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--run-orthorectification", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--run-pansharpen", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--run-cloud-mask", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--run-alignment", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--run-radiometric-normalization", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--skip-flaash", action="store_true")
    parser.add_argument("--existing-flaash-input")
    parser.add_argument("--existing-mul-ortho-input")
    parser.add_argument("--existing-pan-ortho-input")
    parser.add_argument("--fetch-atmosphere-source", choices=["auto", "nasa_power", "modis_gee"], default="auto")
    parser.add_argument("--fetch-atmosphere-grid-size", type=int, default=3)
    parser.add_argument("--fetch-atmosphere-search-days", type=int, default=1)
    parser.add_argument("--fetch-atmosphere-timeout-s", type=float, default=30.0)
    parser.add_argument("--fetch-atmosphere-power-endpoint", default="https://power.larc.nasa.gov/api/temporal/daily/point")
    parser.add_argument("--fetch-atmosphere-ee-project")
    parser.add_argument("--fetch-atmosphere-authenticate", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--fetch-atmosphere-env-file")
    parser.add_argument("--fetch-atmosphere-hours-window", type=int, default=24)
    parser.add_argument("--radiometric-normalization-method", default="spectralmatch")
    parser.add_argument("--radiometric-normalization-kwargs-json")
    parser.add_argument("--cloud-mask-command")
    parser.add_argument("--cloud-mask-method", choices=["omnicloudmask"], default="omnicloudmask")
    parser.add_argument("--cloud-mask-red-band-index", type=int, default=5)
    parser.add_argument("--cloud-mask-green-band-index", type=int, default=3)
    parser.add_argument("--cloud-mask-nir-band-index", type=int, default=7)
    parser.add_argument("--cloud-mask-classes", default="1,2,3")
    parser.add_argument("--cloud-buffer-pixels", type=int, default=10)
    parser.add_argument("--cloud-mask-inference-resolution-m", type=float, default=10.0)
    parser.add_argument("--cloud-mask-omnicloud-kwargs-json")
    parser.add_argument("--cloud-mask-output-suffix", default="_cloudmasked")
    parser.add_argument("--cloud-mask-mask-suffix", default="_cloudmask")
    parser.add_argument("--alignment-fixed-image")
    parser.add_argument("--alignment-output-suffix", default="_aligned")
    parser.add_argument("--alignment-moving-band-index", type=int, default=0)
    parser.add_argument("--alignment-fixed-band-index", type=int, default=0)
    parser.add_argument("--alignment-registration-mode", choices=["default", "structural_wv3_lidar"], default="structural_wv3_lidar")
    parser.add_argument("--alignment-parameter-map", default="rigid")
    parser.add_argument("--alignment-no-tiling", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--alignment-tile-size", type=int, default=1000)
    parser.add_argument("--alignment-tile-buffer", type=int, default=100)
    parser.add_argument("--alignment-clip-fixed-to-moving", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--alignment-enforce-mutual-valid-mask", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--alignment-output-on-moving-grid", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--alignment-moving-nodata", type=float, default=-9999)
    parser.add_argument("--alignment-fixed-nodata", type=float, default=-9999)
    parser.add_argument("--alignment-output-nodata", type=float)
    parser.add_argument("--alignment-min-valid-fraction", type=float, default=0.002)
    parser.add_argument("--alignment-temp-dir")
    parser.add_argument("--alignment-keep-temp-dir", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--alignment-log-to-console", action=argparse.BooleanOptionalAction, default=False)
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """Parse CLI/config arguments, validate them, and execute the workflow."""
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument("--config-yaml")
    config_args, _ = config_parser.parse_known_args(argv)

    config_defaults: Dict = {}
    if config_args.config_yaml:
        config_defaults = _normalize_config_defaults(_load_yaml_config(config_args.config_yaml))

    parser = build_parser()
    if config_defaults:
        parser.set_defaults(**config_defaults)
    args, unknown_args = parser.parse_known_args(argv)
    _apply_unknown_prefixed_args(args, unknown_args)

    if not args.input_file_glob:
        parser.error("--input-file-glob is required (via CLI or --config-yaml).")
    if not args.dem_file_path and not args.existing_mul_ortho_input:
        parser.error("--dem-file-path is required unless --existing-mul-ortho-input is provided.")
    if args.run_pansharpen and not args.run_orthorectification and args.last_run_step not in {"orthorectification", "pansharpen", "cloud_mask", "alignment", "radiometric_normalization"} and not args.existing_mul_ortho_input:
        parser.error("--run-pansharpen requires orthorectified inputs, --run-orthorectification, or a later --last-run-step.")
    if args.run_pansharpen and args.existing_mul_ortho_input and not args.existing_pan_ortho_input:
        parser.error("--existing-pan-ortho-input is required when using --existing-mul-ortho-input with pansharpen.")
    if args.run_alignment and not args.alignment_fixed_image:
        parser.error("--alignment-fixed-image is required when --run-alignment is enabled.")
    if args.run_alignment and args.alignment_fixed_image and not os.path.isfile(args.alignment_fixed_image):
        parser.error(f"--alignment-fixed-image does not exist: {args.alignment_fixed_image}")
    if args.atmospheric_method == "flaash" and args.run_atmospheric_correction and not args.skip_flaash and not args.envi_engine_path:
        parser.error("--envi-engine-path is required when running FLAASH.")
    if args.run_cloud_mask and args.cloud_mask_method and args.cloud_mask_command:
        parser.error("Use either --cloud-mask-method or --cloud-mask-command, not both.")
    if args.last_run_step == "alignment" and not args.run_alignment:
        parser.error("--last-run-step alignment requires --run-alignment.")
    if args.last_run_step == "radiometric_normalization" and not args.run_radiometric_normalization:
        parser.error("--last-run-step radiometric_normalization requires --run-radiometric-normalization.")
    if args.cloud_mask_inference_resolution_m <= 0:
        parser.error("--cloud-mask-inference-resolution-m must be > 0.")
    if args.fetch_atmosphere_grid_size < 1:
        parser.error("--fetch-atmosphere-grid-size must be >= 1.")
    if args.fetch_atmosphere_search_days < 0:
        parser.error("--fetch-atmosphere-search-days must be >= 0.")
    if args.fetch_atmosphere_timeout_s <= 0:
        parser.error("--fetch-atmosphere-timeout-s must be > 0.")
    if args.alignment_tile_size <= 0:
        parser.error("--alignment-tile-size must be > 0.")
    if args.alignment_tile_buffer < 0:
        parser.error("--alignment-tile-buffer must be >= 0.")
    try:
        _parse_json_dict(args.cloud_mask_omnicloud_kwargs_json)
    except Exception as exc:
        parser.error(f"Invalid --cloud-mask-omnicloud-kwargs-json: {exc}")
    try:
        _parse_json_dict(args.radiometric_normalization_kwargs_json)
    except Exception as exc:
        parser.error(f"Invalid --radiometric-normalization-kwargs-json: {exc}")
    return run_workflow(args)


if __name__ == "__main__":
    sys.exit(main())
