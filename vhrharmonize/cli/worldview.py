#!/usr/bin/env python3
"""CLI wrapper for the WorldView preprocessing workflow."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import fnmatch
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional

import geopandas as gpd
from wcmatch import glob

from vhrharmonize.cli.cli_helpers import load_yaml_config
from vhrharmonize.io.geospatial import calculate_raster_overviews, get_image_percentile_value, shp_to_gpkg
from vhrharmonize.io.workflow_utils import (
    plan_step_outputs,
    resolve_output_dir,
    resolve_temp_dir,
    resolve_relative_to_input,
)
from vhrharmonize.preprocess.atmospheric_correction import run_flaash, run_py6s
from vhrharmonize.preprocess.alignment import align_image_pair
from vhrharmonize.preprocess.cloudmasking import cloudmask_raster
from vhrharmonize.preprocess.fetch_external_data import (
    DEFAULT_OPENTOPOGRAPHY_DEMTYPE,
    DEFAULT_OPENTOPOGRAPHY_GLOBALDEM_ENDPOINT,
    download_opentopography_dem_for_bbox,
    fetch_modis_water_vapor_for_bbox,
    fetch_power_atmosphere_for_bbox,
)
from vhrharmonize.preprocess.helpers import log
from vhrharmonize.preprocess.orthorectification import (
    gcp_refined_rpc_orthorectification,
    resolve_output_resolution_for_crs,
)
from vhrharmonize.preprocess.pansharpening import pansharpen_image
from vhrharmonize.preprocess.radiometric_normalization import radiometric_normalization
from vhrharmonize.providers.worldview import (
    WorldViewImage,
    WorldViewScene,
    load_worldview_scenes_from_tif_files,
)

RASTER_STEP_ORDER = [
    "raw",
    "atmospheric_correction",
    "orthorectification",
    "pansharpen",
    "cloud_mask",
    "alignment",
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
    dem_file_path: Optional[str] = None
    fetch_atmosphere_result: Optional[Dict] = None
    py6s_effective_params: Optional[Dict] = None
    py6s_auto_atmos_estimate: Optional[Dict] = None
    alignment_result: Optional[object] = None
    cloud_mask_pixel_count: Optional[int] = None
    cloud_mask_path: Optional[str] = None


def _require_scene_image(scene: WorldViewScene, role: str) -> WorldViewImage:
    """Return a required scene image.
    Args:
        scene: Scene to query.
        role: Requested image role.
    Returns:
        Matching scene image.
    """
    image = scene.get_image(role)
    if image is None:
        raise ValueError(f"WorldView scene is missing required {role} image: {scene.scene_id}_{scene.catalog_id}")
    return image


def _get_worldview_scene_step_path(scene: WorldViewScene, role: str, step_name: str) -> str:
    """Return a stored scene step path.
    Args:
        scene: Scene to query.
        role: Requested image role.
        step_name: Workflow step name.
    Returns:
        Stored step output path.
    """
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
    """Store a scene step path.
    Args:
        scene: Scene to update.
        role: Image role to update.
        step_name: Workflow step name.
        output_path: Output path to store.
    Returns:
        None.
    """
    image = _require_scene_image(scene, role)
    image.step_file_paths[step_name] = output_path


def _parse_filter_basenames(raw_values: Optional[List[str]]) -> List[str]:
    """Normalize optional basename filters.
    Args:
        raw_values: Repeated or comma-delimited basename filters.
    Returns:
        Flattened basename filter list.
    """
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
    """Load a scene footprint bounding box in WGS84.
    Args:
        shp_path: Scene footprint shapefile path.
    Returns:
        Bounding box as min lon, min lat, max lon, max lat.
    """
    gdf = gpd.read_file(shp_path)
    if gdf.empty:
        raise ValueError(f"Empty scene footprint shapefile: {shp_path}")
    if gdf.crs is None:
        raise ValueError(f"Scene footprint shapefile has no CRS: {shp_path}")
    gdf = gdf.to_crs(epsg=4326)
    minx, miny, maxx, maxy = gdf.total_bounds
    return float(minx), float(miny), float(maxx), float(maxy)


def _parse_int_csv(raw_values: str) -> List[int]:
    """Parse comma-delimited integers.
    Args:
        raw_values: Comma-delimited integer string.
    Returns:
        Parsed integer values.
    """
    parsed = []
    for value in raw_values.split(","):
        value = value.strip()
        if value:
            parsed.append(int(value))
    return parsed


def _parse_json_dict(raw_json: Optional[object]) -> Dict:
    """Parse an optional JSON object value.
    Args:
        raw_json: Optional JSON string or mapping.
    Returns:
        Parsed dictionary value.
    """
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
    transform_key: object = None,
) -> Dict:
    """Collect namespace values with a shared prefix.
    Args:
        namespace: Parsed CLI namespace.
        prefix: Prefix to strip from matching keys.
        transform_key: Optional callable used to transform stripped keys.
    Returns:
        Collected keyword arguments.
    """
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
    """Build radiometric normalization keyword arguments.
    Args:
        args: Parsed CLI arguments.
    Returns:
        Radiometric normalization keyword arguments.
    """
    radiometric_kwargs = _parse_json_dict(args.radiometric_normalization_kwargs_json)
    radiometric_kwargs.update(_collect_prefixed_kwargs(args, "match_"))
    return radiometric_kwargs


def _coerce_unknown_arg_value(raw_value: str) -> object:
    """Coerce an unknown CLI value into a Python scalar.
    Args:
        raw_value: Raw CLI token value.
    Returns:
        Parsed Python value.
    """
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
    """Apply passthrough CLI arguments to the namespace.
    Args:
        args: Parsed CLI namespace to mutate.
        unknown_args: Unknown CLI tokens to interpret.
    Returns:
        None.
    """
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


def _load_worldview_yaml_config(config_yaml_path: str) -> Dict:
    """Load and flatten a WorldView workflow config.
    Args:
        config_yaml_path: YAML config file path.
    Returns:
        Flattened config dictionary.
    """
    loaded = load_yaml_config(config_yaml_path)

    def _flatten_mapping(mapping: Dict, out: Dict) -> None:
        """Flatten nested config mappings.
        Args:
            mapping: Mapping to flatten.
            out: Output mapping to populate.
        Returns:
            None.
        """
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
    """Normalize YAML-derived workflow defaults.
    Args:
        config_defaults: Raw config defaults mapping.
    Returns:
        Normalized config defaults mapping.
    """
    normalized = dict(config_defaults)
    for list_key in ("input_file_glob", "input_dir", "filter_basename"):
        if list_key in normalized and isinstance(normalized[list_key], str):
            normalized[list_key] = [normalized[list_key]]
    if "input_dir" in normalized and "input_file_glob" not in normalized:
        normalized["input_file_glob"] = normalized.pop("input_dir")
    return normalized


def _resolve_fetch_atmosphere_source(args: argparse.Namespace) -> str:
    """Resolve the active atmosphere source.
    Args:
        args: Parsed CLI arguments.
    Returns:
        Atmosphere source name.
    """
    if args.fetch_atmosphere_source != "auto":
        return args.fetch_atmosphere_source
    if args.atmospheric_method == "flaash":
        return "modis_gee"
    return "nasa_power"


def _resolve_scene_dem_file_path(state: SceneWorkflowState, args: argparse.Namespace) -> Optional[str]:
    """Resolve the DEM path for a scene.
    Args:
        state: Scene workflow state.
        args: Parsed CLI arguments.
    Returns:
        Resolved DEM path or None.
    """
    if state.dem_file_path:
        return state.dem_file_path
    if args.dem_file_path in (None, ""):
        return None

    mul_image = _require_scene_image(state.scene, "mul")
    mul_folder = os.path.dirname(mul_image.tif_file)
    dem_value = str(args.dem_file_path).strip()
    if dem_value.lower() != "online":
        resolved_dem_path = resolve_relative_to_input(dem_value, mul_folder)
        state.dem_file_path = resolved_dem_path
        log(
            f"Using DEM {os.path.basename(resolved_dem_path)}",
            enabled=args.log_to_console,
            step="dem",
            scene_basename=state.scene.primary_basename,
        )
        return resolved_dem_path

    dem_dir = os.path.join(state.step_dirs["temp_root"], "dem")
    os.makedirs(dem_dir, exist_ok=True)
    dem_output_path = os.path.join(dem_dir, f"{mul_image.basename}_dem.tif")
    if os.path.isfile(dem_output_path):
        log(
            f"Reusing DEM {os.path.basename(dem_output_path)}",
            enabled=args.log_to_console,
            step="dem",
            scene_basename=state.scene.primary_basename,
        )
    else:
        download_opentopography_dem_for_bbox(
            min_lon=state.scene_bbox_wgs84[0],
            min_lat=state.scene_bbox_wgs84[1],
            max_lon=state.scene_bbox_wgs84[2],
            max_lat=state.scene_bbox_wgs84[3],
            output_tif_path=dem_output_path,
            api_key=args.dem_online_api_key,
            demtype=args.dem_online_source,
            endpoint=args.dem_online_api_endpoint,
            timeout_s=args.dem_online_timeout_s,
            log_to_console=args.log_to_console,
            scene_basename=state.scene.primary_basename,
        )
        log(
            f"Downloaded DEM {os.path.basename(dem_output_path)}",
            enabled=args.log_to_console,
            step="dem",
            scene_basename=state.scene.primary_basename,
        )
    state.dem_file_path = dem_output_path
    return dem_output_path


def _write_json(path: str, payload: Dict) -> None:
    """Write a JSON file.
    Args:
        path: Output JSON path.
        payload: JSON-serializable payload.
    Returns:
        None.
    """
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def _read_json(path: str) -> Dict:
    """Read a JSON object file.
    Args:
        path: Input JSON path.
    Returns:
        Parsed JSON dictionary.
    """
    with open(path, "r", encoding="utf-8") as handle:
        loaded = json.load(handle)
    if not isinstance(loaded, dict):
        raise ValueError(f"Expected JSON object at {path}")
    return loaded


def _collect_input_files(input_file_globs: List[str]) -> List[str]:
    """Collect input files from glob patterns.
    Args:
        input_file_globs: Glob patterns to scan.
    Returns:
        Unique absolute file paths.
    """
    matched_files: List[str] = []
    for pattern in input_file_globs:
        matches = glob.glob(pattern, flags=glob.GLOBSTAR | glob.BRACE | glob.EXTGLOB | glob.GLOBTILDE | glob.GLOBSTARLONG | glob.NEGATE)
        matched_files.extend(path for path in matches if os.path.isfile(path))
    return sorted({os.path.abspath(path) for path in matched_files})


def _resolve_concurrent_processing(value: object) -> int:
    """Resolve the scene concurrency setting.
    Args:
        value: Raw concurrency config value.
    Returns:
        Resolved worker count.
    """
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized == "num_cpu":
            return max(1, os.cpu_count() or 1)
        try:
            resolved = int(normalized)
        except ValueError as exc:
            raise ValueError("concurrent_processing must be an integer or 'num_cpu'.") from exc
    else:
        try:
            resolved = int(value)
        except (TypeError, ValueError) as exc:
            raise ValueError("concurrent_processing must be an integer or 'num_cpu'.") from exc
    if resolved < 1:
        raise ValueError("concurrent_processing must be >= 1.")
    return resolved


def _short_path(path: str) -> str:
    """Return a shortened display path.
    Args:
        path: Full file path.
    Returns:
        Path basename.
    """
    return os.path.basename(path)


def _short_paths(paths: List[str]) -> str:
    """Join shortened display paths.
    Args:
        paths: Full file paths.
    Returns:
        Comma-delimited basenames.
    """
    return ", ".join(_short_path(path) for path in paths)


def _log_step_plan(
    step: str,
    *,
    inputs: List[str] | None = None,
    outputs: List[str] | None = None,
    message: str | None = None,
    enabled: bool = False,
    scene_basename: str | None = None,
) -> None:
    """Log a concise step plan message.
    Args:
        step: Step name for the log prefix.
        inputs: Optional input paths.
        outputs: Optional output paths.
        message: Optional human-readable message.
        enabled: Whether logging is enabled.
        scene_basename: Optional scene basename for the log prefix.
    Returns:
        None.
    """
    parts: List[str] = []
    if message:
        parts.append(message)
    if inputs:
        parts.append(f"in={_short_paths(inputs)}")
    if outputs:
        parts.append(f"out={_short_paths(outputs)}")
    log(" | ".join(parts), enabled=enabled, step=step, scene_basename=scene_basename)


def _resolve_step_save_dir(
    save_value: Optional[str],
    *,
    temp_root: str,
    output_root: str,
    relative_base_folder: str,
) -> str:
    """Resolve a configured step save directory.
    Args:
        save_value: Step save target configuration value.
        temp_root: Temp root directory.
        output_root: Output root directory.
        relative_base_folder: Base folder used for relative resolution.
    Returns:
        Resolved step save directory.
    """
    save_mode = (save_value or "$temp").strip()
    if save_mode == "$temp":
        resolved_dir = temp_root
    elif save_mode.startswith("$temp/"):
        resolved_dir = os.path.join(temp_root, save_mode[len("$temp/"):])
    elif save_mode == "$output":
        resolved_dir = output_root
    elif save_mode.startswith("$output/"):
        resolved_dir = os.path.join(output_root, save_mode[len("$output/"):])
    elif save_mode.startswith("./"):
        resolved_dir = resolve_relative_to_input(save_mode, relative_base_folder)
    elif os.path.isabs(save_mode):
        resolved_dir = save_mode
    else:
        resolved_dir = os.path.abspath(save_mode)
    os.makedirs(resolved_dir, exist_ok=True)
    return resolved_dir


def _is_temp_save_value(save_value: Optional[str]) -> bool:
    """Return whether a save target points into temp storage.
    Args:
        save_value: Step save target configuration value.
    Returns:
        True when the save target uses the temp root.
    """
    normalized = (save_value or "$temp").strip()
    return normalized == "$temp" or normalized.startswith("$temp/")


def _get_last_enabled_raster_step(args: argparse.Namespace) -> str:
    """Return the last enabled raster step.
    Args:
        args: Parsed CLI arguments.
    Returns:
        Last enabled raster step name.
    """
    if args.run_alignment:
        return "alignment"
    if args.run_cloud_mask:
        return "cloud_mask"
    if args.run_pansharpen:
        return "pansharpen"
    if args.run_orthorectification:
        return "orthorectification"
    if args.run_atmospheric_correction:
        return "atmospheric_correction"
    return "raw"


def _step_outputs_exist(output_paths: List[str]) -> bool:
    """Return whether all output paths exist.
    Args:
        output_paths: Output file paths to check.
    Returns:
        True when every output path exists.
    """
    return bool(output_paths) and all(os.path.exists(path) for path in output_paths)


def _resolve_scene_step_dirs(args: argparse.Namespace, scene: WorldViewScene) -> Dict[str, str]:
    """Resolve per-scene step directories.
    Args:
        args: Parsed CLI arguments.
        scene: Scene being prepared.
    Returns:
        Mapping of step names to resolved directories.
    """
    mul_image = _require_scene_image(scene, "mul")
    relative_output_base = os.path.dirname(mul_image.tif_file)
    default_output_root = os.path.normpath(os.path.join(relative_output_base, "..", "Processed"))
    resolved_temp_root = resolve_temp_dir(args.temp_dir, input_folder=relative_output_base)
    resolved_output_root = (
        resolve_relative_to_input(args.output_dir, relative_output_base)
        if args.output_dir not in (None, "")
        else default_output_root
    )
    os.makedirs(resolved_output_root, exist_ok=True)
    step_dirs: Dict[str, str] = {
        "temp_root": resolved_temp_root,
        "output_root": resolved_output_root,
        "scene_work": resolve_output_dir(
            None,
            temp_dir=resolved_temp_root,
            step_name="shared",
        ),
    }

    if args.run_fetch_atmosphere:
        step_dirs["fetch_atmosphere"] = _resolve_step_save_dir(
            args.save_fetch_atmosphere,
            temp_root=resolved_temp_root,
            output_root=resolved_output_root,
            relative_base_folder=relative_output_base,
        )
    if args.run_atmospheric_correction:
        step_dirs["atmospheric_correction"] = _resolve_step_save_dir(
            args.save_atmospheric_correction,
            temp_root=resolved_temp_root,
            output_root=resolved_output_root,
            relative_base_folder=relative_output_base,
        )
    if args.run_orthorectification:
        step_dirs["orthorectification"] = _resolve_step_save_dir(
            args.save_orthorectification,
            temp_root=resolved_temp_root,
            output_root=resolved_output_root,
            relative_base_folder=relative_output_base,
        )
    if args.run_pansharpen:
        step_dirs["pansharpen"] = _resolve_step_save_dir(
            args.save_pansharpen,
            temp_root=resolved_temp_root,
            output_root=resolved_output_root,
            relative_base_folder=relative_output_base,
        )
    if args.run_cloud_mask:
        step_dirs["cloud_mask"] = _resolve_step_save_dir(
            args.save_cloud_mask,
            temp_root=resolved_temp_root,
            output_root=resolved_output_root,
            relative_base_folder=relative_output_base,
        )
    if args.run_alignment:
        step_dirs["alignment"] = _resolve_step_save_dir(
            args.save_alignment,
            temp_root=resolved_temp_root,
            output_root=resolved_output_root,
            relative_base_folder=relative_output_base,
        )
    if args.run_radiometric_normalization:
        step_dirs["radiometric_normalization"] = _resolve_step_save_dir(
            args.save_radiometric_normalization,
            temp_root=resolved_temp_root,
            output_root=resolved_output_root,
            relative_base_folder=relative_output_base,
        )
    return step_dirs


def _get_atmospheric_extension(args: argparse.Namespace) -> str:
    """Return the atmospheric output extension.
    Args:
        args: Parsed CLI arguments.
    Returns:
        Atmospheric output file extension.
    """
    return ".dat" if args.atmospheric_method == "flaash" else ".tif"


def _get_expected_scene_step_outputs(state: SceneWorkflowState, args: argparse.Namespace) -> Dict[str, List[str]]:
    """Plan expected outputs for all enabled scene steps.
    Args:
        state: Scene workflow state.
        args: Parsed CLI arguments.
    Returns:
        Mapping of step names to expected output paths.
    """
    mul_image = _require_scene_image(state.scene, "mul")
    expected_outputs: Dict[str, List[str]] = {
        "raw": [mul_image.tif_file],
    }
    current_mul_outputs = [mul_image.tif_file]

    if args.run_fetch_atmosphere:
        expected_outputs["fetch_atmosphere"] = plan_step_outputs(
            [mul_image.tif_file],
            output_dir=state.step_dirs["fetch_atmosphere"],
            suffix=args.fetch_atmosphere_output_suffix,
            extension=".json",
            skip_existing=False,
        ).output_paths

    if args.run_atmospheric_correction:
        current_mul_outputs = plan_step_outputs(
            current_mul_outputs,
            output_dir=state.step_dirs["atmospheric_correction"],
            suffix=args.atmospheric_correction_output_suffix,
            extension=_get_atmospheric_extension(args),
            skip_existing=False,
        ).output_paths
        expected_outputs["atmospheric_correction"] = list(current_mul_outputs)

    if args.run_orthorectification:
        current_mul_outputs = plan_step_outputs(
            current_mul_outputs,
            output_dir=state.step_dirs["orthorectification"],
            suffix=args.orthorectification_output_suffix,
            skip_existing=False,
        ).output_paths
        expected_outputs["orthorectification"] = list(current_mul_outputs)

        if args.run_pansharpen and state.scene.pan_image is not None:
            expected_outputs["orthorectification_pan"] = plan_step_outputs(
                [state.scene.pan_image.tif_file],
                output_dir=state.step_dirs["orthorectification"],
                suffix=args.orthorectification_pan_output_suffix,
                skip_existing=False,
            ).output_paths

    if args.run_pansharpen:
        current_mul_outputs = plan_step_outputs(
            current_mul_outputs,
            output_dir=state.step_dirs["pansharpen"],
            suffix=args.pansharpen_output_suffix,
            skip_existing=False,
        ).output_paths
        expected_outputs["pansharpen"] = list(current_mul_outputs)

    if args.run_cloud_mask:
        expected_outputs["cloud_mask_mask"] = plan_step_outputs(
            current_mul_outputs,
            output_dir=state.step_dirs["cloud_mask"],
            suffix=args.cloud_mask_mask_suffix,
            skip_existing=False,
        ).output_paths
        current_mul_outputs = plan_step_outputs(
            current_mul_outputs,
            output_dir=state.step_dirs["cloud_mask"],
            suffix=args.cloud_mask_output_suffix,
            skip_existing=False,
        ).output_paths
        expected_outputs["cloud_mask"] = list(current_mul_outputs)

    if args.run_alignment:
        current_mul_outputs = plan_step_outputs(
            current_mul_outputs,
            output_dir=state.step_dirs["alignment"],
            suffix=args.alignment_output_suffix,
            skip_existing=False,
        ).output_paths
        expected_outputs["alignment"] = list(current_mul_outputs)

    expected_outputs["final_raster"] = list(current_mul_outputs)
    return expected_outputs


def _initialize_scene_state(scene: WorldViewScene, args: argparse.Namespace) -> SceneWorkflowState:
    """Initialize workflow state for a scene.
    Args:
        scene: Scene to initialize.
        args: Parsed CLI arguments.
    Returns:
        Initialized scene workflow state.
    """
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
    return state


def _register_step_outputs(
    state: SceneWorkflowState,
    step_name: str,
    output_paths: List[str],
    *,
    image_role: str = "mul",
) -> List[str]:
    """Register step outputs on the scene state.
    Args:
        state: Scene workflow state to update.
        step_name: Workflow step name.
        output_paths: Step output paths to register.
        image_role: Image role whose stored step path should be updated.
    Returns:
        Registered output paths.
    """
    state.scene.step_outputs[step_name] = list(output_paths)
    if output_paths:
        _set_worldview_scene_step_path(state.scene, image_role, step_name, output_paths[0])
    return output_paths


def _mark_scene_complete_from_existing_output(state: SceneWorkflowState, args: argparse.Namespace) -> None:
    """Mark a scene as complete from existing outputs.
    Args:
        state: Scene workflow state to update.
        args: Parsed CLI arguments.
    Returns:
        None.
    """
    expected_outputs = _get_expected_scene_step_outputs(state, args)
    last_step = _get_last_enabled_raster_step(args)
    existing_output = expected_outputs[last_step][0]
    state.current_files = [existing_output]
    state.current_step = last_step
    state.scene.step_outputs[last_step] = [existing_output]
    _set_worldview_scene_step_path(state.scene, "mul", last_step, existing_output)


def _scene_skip_required_outputs(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    """Return required outputs for scene-level skipping.
    Args:
        state: Scene workflow state.
        args: Parsed CLI arguments.
    Returns:
        Required non-temp output paths.
    """
    expected_outputs = _get_expected_scene_step_outputs(state, args)
    required_outputs: List[str] = []
    if args.run_fetch_atmosphere and not _is_temp_save_value(args.save_fetch_atmosphere):
        required_outputs.extend(expected_outputs.get("fetch_atmosphere", []))
    if args.run_atmospheric_correction and not _is_temp_save_value(args.save_atmospheric_correction):
        required_outputs.extend(expected_outputs.get("atmospheric_correction", []))
    if args.run_orthorectification and not _is_temp_save_value(args.save_orthorectification):
        required_outputs.extend(expected_outputs.get("orthorectification", []))
        if args.run_pansharpen:
            required_outputs.extend(expected_outputs.get("orthorectification_pan", []))
    if args.run_pansharpen and not _is_temp_save_value(args.save_pansharpen):
        required_outputs.extend(expected_outputs.get("pansharpen", []))
    if args.run_cloud_mask and not _is_temp_save_value(args.save_cloud_mask):
        required_outputs.extend(expected_outputs.get("cloud_mask", []))
        required_outputs.extend(expected_outputs.get("cloud_mask_mask", []))
    if args.run_alignment and not _is_temp_save_value(args.save_alignment):
        required_outputs.extend(expected_outputs.get("alignment", []))
    return required_outputs


def _normalize_group_by_basename_spec(raw_spec: object) -> object:
    """Normalize a radiometric grouping specification.
    Args:
        raw_spec: Raw grouping specification value.
    Returns:
        Normalized grouping specification.
    """
    if raw_spec in (None, "", []):
        return None
    if isinstance(raw_spec, str):
        stripped = raw_spec.strip()
        if stripped.startswith("["):
            return _normalize_group_by_basename_spec(json.loads(stripped))
        return raw_spec
    if not isinstance(raw_spec, list):
        raise ValueError("group_by_basename must be a string or nested list of strings.")
    return [_normalize_group_by_basename_spec(item) for item in raw_spec]


def _match_radiometric_input_patterns(pattern: str, available_paths: List[str]) -> List[str]:
    """Match radiometric input patterns against available paths.
    Args:
        pattern: Glob-like pattern to match.
        available_paths: Available scene output paths.
    Returns:
        Matching scene output paths.
    """
    matches = [
        path
        for path in available_paths
        if fnmatch.fnmatch(os.path.basename(path), pattern) or fnmatch.fnmatch(path, pattern)
    ]
    if not matches:
        raise ValueError(f"No radiometric normalization inputs matched pattern: {pattern}")
    return matches


def _dedupe_paths(paths: List[str]) -> List[str]:
    """Remove duplicate paths while preserving order.
    Args:
        paths: Candidate file paths.
    Returns:
        Deduplicated file paths.
    """
    seen = set()
    deduped: List[str] = []
    for path in paths:
        if path in seen:
            continue
        seen.add(path)
        deduped.append(path)
    return deduped


def _resolve_radiometric_group_output_path(
    args: argparse.Namespace,
    output_dir: str,
    *,
    group_label: str,
    is_root: bool,
) -> str:
    """Resolve a radiometric group output path.
    Args:
        args: Parsed CLI arguments.
        output_dir: Base output directory for radiometric products.
        group_label: Group label used in default output naming.
        is_root: Whether this is the root radiometric group.
    Returns:
        Resolved radiometric group output path.
    """
    configured_output = _build_radiometric_kwargs(args).get("shared_output_image_path")
    if configured_output is None and is_root:
        configured_output = args.radiometric_normalization_output
    if configured_output:
        if "$" in configured_output:
            configured_output = configured_output.replace("$", group_label)
        if os.path.isabs(configured_output):
            return configured_output
        return os.path.normpath(os.path.join(output_dir, configured_output))
    return os.path.join(output_dir, f"{group_label}.tif")


def _run_radiometric_group(
    group_spec: object,
    *,
    available_paths: List[str],
    args: argparse.Namespace,
    output_dir: str,
    temp_root: str,
    label_parts: List[int],
) -> str | List[str]:
    """Run or resolve a radiometric group.
    Args:
        group_spec: Group specification string or nested list.
        available_paths: Available scene output paths.
        args: Parsed CLI arguments.
        output_dir: Radiometric output directory.
        temp_root: Temp root used for SpectralMatch temp files.
        label_parts: Hierarchical label parts for nested groups.
    Returns:
        Group output path or matched input paths.
    """
    if isinstance(group_spec, str):
        return _match_radiometric_input_patterns(group_spec, available_paths)
    if not isinstance(group_spec, list) or not group_spec:
        raise ValueError("Each radiometric normalization group must be a non-empty string or list.")

    child_inputs: List[str] = []
    for index, item in enumerate(group_spec):
        child_label_parts = [*label_parts, index]
        child_result = _run_radiometric_group(
            item,
            available_paths=available_paths,
            args=args,
            output_dir=output_dir,
            temp_root=temp_root,
            label_parts=child_label_parts,
        )
        if isinstance(child_result, list):
            child_inputs.extend(child_result)
        else:
            child_inputs.append(child_result)

    child_inputs = _dedupe_paths(child_inputs)
    group_label = "radiometric_root" if not label_parts else "radiometric_group_" + "_".join(str(part) for part in label_parts)
    output_path = _resolve_radiometric_group_output_path(
        args,
        output_dir,
        group_label=group_label,
        is_root=not label_parts,
    )
    if args.run_from_existing and os.path.exists(output_path):
        _log_step_plan(
            "radiometric",
            outputs=[output_path],
            message=f"Skipping group {group_label} because output exists",
            enabled=args.log_to_console,
        )
        return output_path

    radiometric_kwargs = _build_radiometric_kwargs(args)
    radiometric_kwargs.setdefault("shared_input_images", child_inputs)
    radiometric_kwargs.setdefault("shared_output_image_path", output_path)
    radiometric_kwargs.setdefault("shared_temp_dir", os.path.join(temp_root, "spectralmatch"))
    radiometric_kwargs.setdefault("delete_temp_dir", args.keep_temp_dir is not True)
    radiometric_kwargs.setdefault("shared_debug_logs", args.log_to_console)
    radiometric_kwargs.setdefault("shared_output_dtype", args.dtype)
    _log_step_plan(
        "radiometric",
        inputs=child_inputs,
        outputs=[output_path],
        message=f"Running SpectralMatch group {group_label}",
        enabled=args.log_to_console,
    )
    radiometric_normalization(
        method=args.radiometric_normalization_method,
        log_to_console=args.log_to_console,
        **radiometric_kwargs,
    )
    if args.calculate_overviews_radiometric_normalization:
        log("Calculating overviews for step radiometric_normalization", enabled=args.log_to_console, step="overviews")
        calculate_raster_overviews(output_path, args.overview_scales)
    return output_path


def _run_radiometric_normalization_workflow(
    scene_output_paths: List[str],
    *,
    args: argparse.Namespace,
    reference_state: SceneWorkflowState,
) -> Optional[str]:
    """Run grouped radiometric normalization.
    Args:
        scene_output_paths: Per-scene output raster paths.
        args: Parsed CLI arguments.
        reference_state: Reference scene state used for directory resolution.
    Returns:
        Final radiometric output path or None.
    """
    if not args.run_radiometric_normalization:
        return None
    available_paths = _dedupe_paths([str(path) for path in scene_output_paths if str(path)])
    if not available_paths:
        raise ValueError("No scene outputs were available for radiometric normalization.")

    group_spec = _normalize_group_by_basename_spec(args.group_by_basename)
    if group_spec is None:
        group_spec = available_paths
    elif isinstance(group_spec, str):
        group_spec = [group_spec]

    return _run_radiometric_group(
        group_spec,
        available_paths=available_paths,
        args=args,
        output_dir=reference_state.step_dirs["radiometric_normalization"],
        temp_root=reference_state.step_dirs["temp_root"],
        label_parts=[],
    )


def _run_cloud_mask_command(
    command_template: str,
    input_image_path: str,
    output_image_path: str,
    scene_root_path: str,
    image_basename: str,
    *,
    log_to_console: bool = False,
    scene_basename: str | None = None,
) -> None:
    """Run an external cloud mask command template.
    Args:
        command_template: Shell command template to execute.
        input_image_path: Input raster path.
        output_image_path: Output raster path.
        scene_root_path: Scene root directory.
        image_basename: Scene image basename.
        log_to_console: Whether to emit console logs.
        scene_basename: Optional scene basename for log prefixes.
    Returns:
        None.
    """
    command = command_template.format(
        input=input_image_path,
        output=output_image_path,
        scene_root=scene_root_path,
        image_basename=image_basename,
    )
    log(
        "Running external cloud mask command",
        enabled=log_to_console,
        step="cloud_mask",
        scene_basename=scene_basename,
    )
    subprocess.run(command, shell=True, check=True)


def _run_fetch_atmosphere_step(state: SceneWorkflowState, args: argparse.Namespace) -> None:
    """Run the fetch-atmosphere step.
    Args:
        state: Scene workflow state.
        args: Parsed CLI arguments.
    Returns:
        None.
    """
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
        skip_existing=False,
    )
    if args.run_from_existing and _step_outputs_exist(plan.output_paths):
        _log_step_plan(
            "fetch_atmosphere",
            outputs=plan.output_paths,
            message="Skipping because output exists",
            enabled=args.log_to_console,
            scene_basename=state.scene.primary_basename,
        )
        state.fetch_atmosphere_result = _read_json(plan.output_paths[0])
        _register_step_outputs(state, "fetch_atmosphere", plan.output_paths)
        return

    fetch_source = _resolve_fetch_atmosphere_source(args)
    _log_step_plan(
        "fetch_atmosphere",
        inputs=[mul_image.tif_file],
        outputs=plan.output_paths,
        message=f"Fetching atmosphere via {fetch_source}",
        enabled=args.log_to_console,
        scene_basename=state.scene.primary_basename,
    )
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
            scene_basename=state.scene.primary_basename,
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
            scene_basename=state.scene.primary_basename,
        )
        result = estimate.to_dict()
    else:
        raise ValueError(f"Unsupported fetch atmosphere source: {fetch_source}")

    _write_json(plan.pending_output_paths[0], result)
    _log_step_plan(
        "fetch_atmosphere",
        outputs=plan.output_paths,
        message="Wrote atmosphere metadata",
        enabled=args.log_to_console,
        scene_basename=state.scene.primary_basename,
    )
    state.fetch_atmosphere_result = result
    _register_step_outputs(state, "fetch_atmosphere", plan.output_paths)


def _run_atmospheric_correction_step(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    """Run the atmospheric correction step.
    Args:
        state: Scene workflow state.
        args: Parsed CLI arguments.
    Returns:
        Current raster paths after atmospheric correction.
    """
    if not args.run_atmospheric_correction:
        return state.current_files
    mul_image = state.scene.mul_image
    if mul_image is None or mul_image.standardized_metadata is None or mul_image.shp_file is None:
        raise ValueError("WorldView scene is missing multispectral inputs for atmospheric correction.")
    resolved_dem_file_path = _resolve_scene_dem_file_path(state, args)
    if resolved_dem_file_path is None:
        raise ValueError("DEM file path is required for atmospheric correction.")

    plan = plan_step_outputs(
        state.current_files,
        output_dir=state.step_dirs["atmospheric_correction"],
        suffix=args.atmospheric_correction_output_suffix,
        extension=_get_atmospheric_extension(args),
        skip_existing=False,
    )
    if args.run_from_existing and _step_outputs_exist(plan.output_paths):
        _log_step_plan(
            "atmospheric_correction",
            outputs=plan.output_paths,
            message="Skipping because output exists",
            enabled=args.log_to_console,
            scene_basename=state.scene.primary_basename,
        )
        state.current_files = _register_step_outputs(state, "atmospheric_correction", plan.output_paths)
        state.current_step = "atmospheric_correction"
        return state.current_files

    input_raster = plan.pending_input_paths[0]
    output_raster = plan.pending_output_paths[0]
    _log_step_plan(
        "atmospheric_correction",
        inputs=[input_raster],
        outputs=[output_raster],
        message=f"Running {args.atmospheric_method}",
        enabled=args.log_to_console,
        scene_basename=state.scene.primary_basename,
    )

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
            dem_file_path=resolved_dem_file_path,
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
            resolved_dem_file_path,
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
            scene_basename=state.scene.primary_basename,
        )
        state.py6s_effective_params = py6s_result.effective_params
        state.py6s_auto_atmos_estimate = py6s_result.auto_atmos_estimate
    else:
        state.current_files = [input_raster]
        state.current_step = "raw"
        return state.current_files

    state.current_files = _register_step_outputs(state, "atmospheric_correction", plan.output_paths)
    if args.calculate_overviews_atmospheric_correction:
        log(
            "Calculating overviews for step atmospheric_correction",
            enabled=args.log_to_console,
            step="overviews",
            scene_basename=state.scene.primary_basename,
        )
        for output_path in plan.output_paths:
            calculate_raster_overviews(output_path, args.overview_scales)
    state.current_step = "atmospheric_correction"
    return state.current_files


def _run_orthorectification_step(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    """Run the orthorectification step.
    Args:
        state: Scene workflow state.
        args: Parsed CLI arguments.
    Returns:
        Current raster paths after orthorectification.
    """
    if not args.run_orthorectification:
        return state.current_files
    mul_image = state.scene.mul_image
    pan_image = state.scene.pan_image
    if mul_image is None or pan_image is None:
        raise ValueError("WorldView scene is missing multispectral or panchromatic image.")
    resolved_dem_file_path = _resolve_scene_dem_file_path(state, args)
    if resolved_dem_file_path is None:
        raise ValueError("DEM file path is required for orthorectification.")

    if args.existing_mul_ortho_input:
        state.current_files = _register_step_outputs(state, "orthorectification", [args.existing_mul_ortho_input])
    else:
        plan = plan_step_outputs(
            state.current_files,
            output_dir=state.step_dirs["orthorectification"],
            suffix=args.orthorectification_output_suffix,
            skip_existing=False,
        )
        if args.run_from_existing and _step_outputs_exist(plan.output_paths):
            _log_step_plan(
                "orthorectification",
                outputs=plan.output_paths,
                message="Skipping because output exists",
                enabled=args.log_to_console,
                scene_basename=state.scene.primary_basename,
            )
        else:
            _log_step_plan(
                "orthorectification",
                inputs=plan.pending_input_paths,
                outputs=plan.pending_output_paths,
                message=f"Projecting to EPSG:{args.epsg}",
                enabled=args.log_to_console,
                scene_basename=state.scene.primary_basename,
            )
            for input_path, output_path in zip(plan.pending_input_paths, plan.pending_output_paths):
                gcp_refined_rpc_orthorectification(
                    input_path,
                    output_path,
                    resolved_dem_file_path,
                    args.epsg,
                    gcp_geojson_file_path=args.orthorectification_rpc_refinement_geojson,
                    output_nodata_value=args.nodata_value,
                    dtype=args.dtype,
                    output_resolution=resolve_output_resolution_for_crs(
                        args.epsg,
                        mul_image.standardized_metadata.product_resolution,
                    ),
                    log_to_console=args.log_to_console,
                    scene_basename=state.scene.primary_basename,
                )
        state.current_files = _register_step_outputs(state, "orthorectification", plan.output_paths)
        if args.calculate_overviews_orthorectification:
            log(
                "Calculating overviews for step orthorectification",
                enabled=args.log_to_console,
                step="overviews",
                scene_basename=state.scene.primary_basename,
            )
            for output_path in plan.output_paths:
                calculate_raster_overviews(output_path, args.overview_scales)

    if args.run_pansharpen:
        if args.existing_pan_ortho_input:
            state.pan_ortho_path = args.existing_pan_ortho_input
        else:
            pan_plan = plan_step_outputs(
                [pan_image.tif_file],
                output_dir=state.step_dirs["orthorectification"],
                suffix=args.orthorectification_pan_output_suffix,
                skip_existing=False,
            )
            if args.run_from_existing and _step_outputs_exist(pan_plan.output_paths):
                _log_step_plan(
                    "orthorectification_pan",
                    outputs=pan_plan.output_paths,
                    message="Skipping because output exists",
                    enabled=args.log_to_console,
                    scene_basename=state.scene.primary_basename,
                )
            else:
                _log_step_plan(
                    "orthorectification_pan",
                    inputs=pan_plan.pending_input_paths,
                    outputs=pan_plan.pending_output_paths,
                    message=f"Projecting to EPSG:{args.epsg}",
                    enabled=args.log_to_console,
                    scene_basename=state.scene.primary_basename,
                )
                for input_path, output_path in zip(pan_plan.pending_input_paths, pan_plan.pending_output_paths):
                    gcp_refined_rpc_orthorectification(
                        input_path,
                        output_path,
                        resolved_dem_file_path,
                        args.epsg,
                        gcp_geojson_file_path=args.orthorectification_rpc_refinement_geojson,
                        output_nodata_value=args.nodata_value,
                        dtype=args.dtype,
                        output_resolution=resolve_output_resolution_for_crs(
                            args.epsg,
                            pan_image.standardized_metadata.product_resolution,
                        ),
                        log_to_console=args.log_to_console,
                        scene_basename=state.scene.primary_basename,
                    )
            state.pan_ortho_path = pan_plan.output_paths[0]
            _register_step_outputs(state, "orthorectification_pan", pan_plan.output_paths, image_role="pan")
            if args.calculate_overviews_orthorectification:
                log(
                    "Calculating overviews for step orthorectification_pan",
                    enabled=args.log_to_console,
                    step="overviews",
                    scene_basename=state.scene.primary_basename,
                )
                for output_path in pan_plan.output_paths:
                    calculate_raster_overviews(output_path, args.overview_scales)

    state.current_step = "orthorectification"
    return state.current_files


def _run_pansharpen_step(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    """Run the pansharpen step.
    Args:
        state: Scene workflow state.
        args: Parsed CLI arguments.
    Returns:
        Current raster paths after pansharpening.
    """
    if not args.run_pansharpen:
        return state.current_files
    if state.pan_ortho_path is None:
        raise ValueError("Panchromatic orthorectified path is required for pansharpening.")
    plan = plan_step_outputs(
        state.current_files,
        output_dir=state.step_dirs["pansharpen"],
        suffix=args.pansharpen_output_suffix,
        skip_existing=False,
    )
    if args.run_from_existing and _step_outputs_exist(plan.output_paths):
        _log_step_plan(
            "pansharpen",
            outputs=plan.output_paths,
            message="Skipping because output exists",
            enabled=args.log_to_console,
            scene_basename=state.scene.primary_basename,
        )
    else:
        _log_step_plan(
            "pansharpen",
            inputs=plan.pending_input_paths + ([state.pan_ortho_path] if state.pan_ortho_path else []),
            outputs=plan.pending_output_paths,
            message="Running pansharpen",
            enabled=args.log_to_console,
            scene_basename=state.scene.primary_basename,
        )
        for input_path, output_path in zip(plan.pending_input_paths, plan.pending_output_paths):
            pansharpen_image(
                input_path,
                state.pan_ortho_path,
                output_path,
                change_nodata_value=args.nodata_value,
                log_to_console=args.log_to_console,
                scene_basename=state.scene.primary_basename,
            )
    state.current_files = _register_step_outputs(state, "pansharpen", plan.output_paths)
    if args.calculate_overviews_pansharpen:
        log(
            "Calculating overviews for step pansharpen",
            enabled=args.log_to_console,
            step="overviews",
            scene_basename=state.scene.primary_basename,
        )
        for output_path in plan.output_paths:
            calculate_raster_overviews(output_path, args.overview_scales)
    state.current_step = "pansharpen"
    return state.current_files


def _run_cloud_mask_step(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    """Run the cloud mask step.
    Args:
        state: Scene workflow state.
        args: Parsed CLI arguments.
    Returns:
        Current raster paths after cloud masking.
    """
    if not args.run_cloud_mask:
        return state.current_files
    mul_image = state.scene.mul_image
    if mul_image is None:
        return state.current_files

    output_plan = plan_step_outputs(
        state.current_files,
        output_dir=state.step_dirs["cloud_mask"],
        suffix=args.cloud_mask_output_suffix,
        skip_existing=False,
    )
    mask_plan = plan_step_outputs(
        state.current_files,
        output_dir=state.step_dirs["cloud_mask"],
        suffix=args.cloud_mask_mask_suffix,
        skip_existing=False,
    )

    if args.run_from_existing and _step_outputs_exist(output_plan.output_paths) and _step_outputs_exist(mask_plan.output_paths):
        _log_step_plan(
            "cloud_mask",
            outputs=output_plan.output_paths + mask_plan.output_paths,
            message="Skipping because output exists",
            enabled=args.log_to_console,
            scene_basename=state.scene.primary_basename,
        )
    elif args.cloud_mask_command:
        _log_step_plan(
            "cloud_mask",
            inputs=output_plan.pending_input_paths,
            outputs=output_plan.pending_output_paths,
            message="Running external cloud mask command",
            enabled=args.log_to_console,
            scene_basename=state.scene.primary_basename,
        )
        for input_path, output_path in zip(output_plan.pending_input_paths, output_plan.pending_output_paths):
            _run_cloud_mask_command(
                args.cloud_mask_command,
                input_path,
                output_path,
                state.scene.root_folder_path,
                mul_image.basename,
                log_to_console=args.log_to_console,
                scene_basename=state.scene.primary_basename,
            )
    else:
        cloud_classes = _parse_int_csv(args.cloud_mask_classes)
        omnicloud_kwargs = _parse_json_dict(args.cloud_mask_omnicloud_kwargs_json)
        _log_step_plan(
            "cloud_mask",
            inputs=output_plan.pending_input_paths,
            outputs=output_plan.pending_output_paths + mask_plan.pending_output_paths,
            message=f"Running OmniCloudMask classes={cloud_classes} buffer={args.cloud_buffer_pixels}",
            enabled=args.log_to_console,
            scene_basename=state.scene.primary_basename,
        )
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
                scene_basename=state.scene.primary_basename,
            )
            state.cloud_mask_pixel_count = cloudmask_result.mask_pixel_count
            state.cloud_mask_path = cloudmask_result.output_mask_path

    state.current_files = _register_step_outputs(state, "cloud_mask", output_plan.output_paths)
    if args.calculate_overviews_cloud_mask:
        log(
            "Calculating overviews for step cloud_mask",
            enabled=args.log_to_console,
            step="overviews",
            scene_basename=state.scene.primary_basename,
        )
        for output_path in output_plan.output_paths:
            calculate_raster_overviews(output_path, args.overview_scales)
    if mask_plan.output_paths:
        _register_step_outputs(state, "cloud_mask_mask", mask_plan.output_paths)
        state.cloud_mask_path = mask_plan.output_paths[0]
    state.current_step = "cloud_mask"
    return state.current_files


def _run_alignment_step(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    """Run the alignment step.
    Args:
        state: Scene workflow state.
        args: Parsed CLI arguments.
    Returns:
        Current raster paths after alignment.
    """
    if not args.run_alignment:
        return state.current_files
    plan = plan_step_outputs(
        state.current_files,
        output_dir=state.step_dirs["alignment"],
        suffix=args.alignment_output_suffix,
        skip_existing=False,
    )
    if args.run_from_existing and _step_outputs_exist(plan.output_paths):
        _log_step_plan(
            "alignment",
            outputs=plan.output_paths,
            message="Skipping because output exists",
            enabled=args.log_to_console,
            scene_basename=state.scene.primary_basename,
        )
    else:
        _log_step_plan(
            "alignment",
            inputs=plan.pending_input_paths + [args.alignment_fixed_image],
            outputs=plan.pending_output_paths,
            message=f"Running coregistration split_factor={args.alignment_split_factor}",
            enabled=args.log_to_console,
            scene_basename=state.scene.primary_basename,
        )
        for input_path, output_path in zip(plan.pending_input_paths, plan.pending_output_paths):
            state.alignment_result = align_image_pair(
                moving_image_path=input_path,
                fixed_image_path=args.alignment_fixed_image,
                output_image_path=output_path,
                band_index=args.alignment_band_index,
                moving_band_index=args.alignment_moving_band_index,
                fixed_band_index=args.alignment_fixed_band_index,
                moving_nodata=args.alignment_moving_nodata,
                fixed_nodata=args.alignment_fixed_nodata,
                output_nodata=args.alignment_output_nodata if args.alignment_output_nodata is not None else args.nodata_value,
                min_valid_fraction=args.alignment_min_valid_fraction,
                temp_dir=state.step_dirs["temp_root"],
                keep_temp_dir=args.keep_temp_dir,
                split_factor=args.alignment_split_factor,
                clip_fixed_to_moving=args.alignment_clip_fixed_to_moving,
                output_on_moving_grid=args.alignment_output_on_moving_grid,
                trim_edge_invalid=args.alignment_trim_edge_invalid,
                edge_trim_depth=args.alignment_edge_trim_depth,
                edge_trim_detection_band_index=args.alignment_edge_trim_detection_band_index,
                edge_trim_invalid_below=args.alignment_edge_trim_invalid_below,
                edge_trim_invalid_above=args.alignment_edge_trim_invalid_above,
                enforce_mutual_valid_mask=args.alignment_enforce_mutual_valid_mask,
                use_edge_proxies=args.alignment_use_edge_proxies,
                solve_resolution=args.alignment_solve_resolution,
                log_to_console=args.log_to_console,
                scene_basename=state.scene.primary_basename,
            )
    state.current_files = _register_step_outputs(state, "alignment", plan.output_paths)
    if args.calculate_overviews_alignment:
        log(
            "Calculating overviews for step alignment",
            enabled=args.log_to_console,
            step="overviews",
            scene_basename=state.scene.primary_basename,
        )
        for output_path in plan.output_paths:
            calculate_raster_overviews(output_path, args.overview_scales)
    state.current_step = "alignment"
    return state.current_files


def _final_output_paths(state: SceneWorkflowState, args: argparse.Namespace) -> tuple[str, str]:
    """Resolve final scene output paths.
    Args:
        state: Scene workflow state.
        args: Parsed CLI arguments.
    Returns:
        Final raster path and final metadata report path.
    """
    final_image_path = _get_expected_scene_step_outputs(state, args)["final_raster"][0]
    final_base = os.path.splitext(final_image_path)[0]
    final_metadata_path = f"{final_base}_metadata.json"
    return final_image_path, final_metadata_path


def _scene_final_outputs_complete(state: SceneWorkflowState, args: argparse.Namespace) -> bool:
    """Return whether final scene outputs are complete.
    Args:
        state: Scene workflow state.
        args: Parsed CLI arguments.
    Returns:
        True when the scene has all required outputs.
    """
    expected_outputs = _get_expected_scene_step_outputs(state, args)
    required_outputs = _scene_skip_required_outputs(state, args)
    if required_outputs and not _step_outputs_exist(required_outputs):
        return False
    if args.run_radiometric_normalization:
        return os.path.isfile(expected_outputs["final_raster"][0])
    return True if required_outputs else os.path.isfile(expected_outputs["final_raster"][0])


def _scene_saved_output_paths(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    """Collect saved scene output paths.
    Args:
        state: Scene workflow state.
        args: Parsed CLI arguments.
    Returns:
        Saved non-temp output paths.
    """
    saved_paths: List[str] = []

    def _extend(step_name: str) -> None:
        """Extend saved paths from a step.
        Args:
            step_name: Step name whose outputs should be added.
        Returns:
            None.
        """
        saved_paths.extend(state.scene.step_outputs.get(step_name, []))

    if args.run_fetch_atmosphere and not _is_temp_save_value(args.save_fetch_atmosphere):
        _extend("fetch_atmosphere")
    if args.run_atmospheric_correction and not _is_temp_save_value(args.save_atmospheric_correction):
        _extend("atmospheric_correction")
    if args.run_orthorectification and not _is_temp_save_value(args.save_orthorectification):
        _extend("orthorectification")
        if args.run_pansharpen:
            _extend("orthorectification_pan")
    if args.run_pansharpen and not _is_temp_save_value(args.save_pansharpen):
        _extend("pansharpen")
    if args.run_cloud_mask and not _is_temp_save_value(args.save_cloud_mask):
        _extend("cloud_mask")
        _extend("cloud_mask_mask")
    if args.run_alignment and not _is_temp_save_value(args.save_alignment):
        _extend("alignment")
    return _dedupe_paths([path for path in saved_paths if path])


def _scene_cloud_cover_percent(scene: WorldViewScene) -> Optional[float]:
    """Return the scene cloud cover percent.
    Args:
        scene: Scene to inspect.
    Returns:
        Cloud cover percent or None.
    """
    image = scene.mul_image or scene.pan_image
    if image is None or image.standardized_metadata is None:
        return None
    return image.standardized_metadata.cloud_cover


def _delete_files(paths: List[str]) -> None:
    """Delete files from a path list.
    Args:
        paths: File paths to delete if they exist.
    Returns:
        None.
    """
    for path in _dedupe_paths([str(path) for path in paths if str(path)]):
        if os.path.isfile(path):
            os.remove(path)


def _scene_temp_cleanup_paths(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    """Collect temp cleanup paths for a scene.
    Args:
        state: Scene workflow state.
        args: Parsed CLI arguments.
    Returns:
        Temp-backed file paths safe to delete after scene completion.
    """
    temp_paths: List[str] = []

    def _extend(step_name: str) -> None:
        """Extend temp paths from a step.
        Args:
            step_name: Step name whose outputs should be added.
        Returns:
            None.
        """
        temp_paths.extend(state.scene.step_outputs.get(step_name, []))

    if args.run_fetch_atmosphere and _is_temp_save_value(args.save_fetch_atmosphere):
        _extend("fetch_atmosphere")
    if args.run_atmospheric_correction and _is_temp_save_value(args.save_atmospheric_correction):
        _extend("atmospheric_correction")
        if args.atmospheric_method == "flaash":
            temp_paths.extend([f"{path}.params.txt" for path in state.scene.step_outputs.get("atmospheric_correction", [])])
    if args.run_orthorectification and _is_temp_save_value(args.save_orthorectification):
        _extend("orthorectification")
        if args.run_pansharpen:
            _extend("orthorectification_pan")
    if args.run_pansharpen and _is_temp_save_value(args.save_pansharpen):
        _extend("pansharpen")
    if args.run_cloud_mask and _is_temp_save_value(args.save_cloud_mask):
        _extend("cloud_mask")
        _extend("cloud_mask_mask")
    if args.run_alignment and _is_temp_save_value(args.save_alignment):
        _extend("alignment")

    mul_image = state.scene.mul_image
    if mul_image is not None:
        temp_paths.append(os.path.join(state.step_dirs["scene_work"], f"{mul_image.basename}.gpkg"))

    return _dedupe_paths(temp_paths)


def _write_scene_report(state: SceneWorkflowState, args: argparse.Namespace, *, scene_started_utc: str) -> None:
    """Write the final scene report.
    Args:
        state: Scene workflow state.
        args: Parsed CLI arguments.
        scene_started_utc: Scene start timestamp in UTC.
    Returns:
        None.
    """
    mul_image = state.scene.mul_image
    pan_image = state.scene.pan_image
    if mul_image is None or pan_image is None:
        raise ValueError("WorldView scene is missing required images for metadata reporting.")
    resolved_dem_file_path = _resolve_scene_dem_file_path(state, args)
    final_scene_path, scene_metadata_path = _final_output_paths(state, args)
    saved_output_paths = _scene_saved_output_paths(state, args)
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
            "dem_file_path": resolved_dem_file_path,
            "dem_file_path_requested": args.dem_file_path,
        },
        "standardized_metadata": {
            "mul": mul_image.standardized_metadata.to_dict() if mul_image.standardized_metadata else None,
            "pan": pan_image.standardized_metadata.to_dict() if pan_image.standardized_metadata else None,
        },
        "workflow": {
            "run_from_existing": args.run_from_existing,
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
            "saved_output_paths": saved_output_paths,
            "final_scene_path": final_scene_path,
            "scene_metadata_path": scene_metadata_path,
        },
        "alignment": (
            {
                "fixed_image": args.alignment_fixed_image,
                "result": (
                    {
                        "output_image_path": state.alignment_result.output_image_path,
                    }
                    if state.alignment_result
                    else None
                ),
            }
        ),
    }
    _write_json(scene_metadata_path, payload)
    state.scene.metadata_report_path = scene_metadata_path


def _write_cloud_cover_skip_report(
    state: SceneWorkflowState,
    args: argparse.Namespace,
    *,
    cloud_cover: float,
) -> None:
    """Write a cloud-cover skip report.
    Args:
        state: Scene workflow state.
        args: Parsed CLI arguments.
        cloud_cover: Scene cloud cover percent.
    Returns:
        None.
    """
    _, scene_metadata_path = _final_output_paths(state, args)
    message = (
        f"Max cloud cover of {cloud_cover:.2f}% does not meet "
        f"max_cloud_cover_to_process of {args.max_cloud_cover_to_process:.2f}%"
    )
    _write_json(scene_metadata_path, {"message": message})
    state.scene.metadata_report_path = scene_metadata_path


def _process_scene(scene: WorldViewScene, args: argparse.Namespace) -> SceneWorkflowState:
    """Process a single WorldView scene.
    Args:
        scene: Scene to process.
        args: Parsed CLI arguments.
    Returns:
        Final scene workflow state.
    """
    state = _initialize_scene_state(scene, args)
    cloud_cover = _scene_cloud_cover_percent(scene)
    if args.max_cloud_cover_to_process is not None and cloud_cover is not None and cloud_cover > args.max_cloud_cover_to_process:
        _write_cloud_cover_skip_report(state, args, cloud_cover=cloud_cover)
        log(
            f"Skipping because cloud cover {cloud_cover:.2f}% exceeds max {args.max_cloud_cover_to_process:.2f}%",
            enabled=args.log_to_console,
            step="workflow",
            scene_basename=scene.primary_basename,
        )
        return state
    log(
        f"Processing {scene.primary_basename or f'{scene.scene_id}_{scene.catalog_id}'}",
        enabled=args.log_to_console,
        step="workflow",
        scene_basename=scene.primary_basename,
    )

    if args.skip_existing and _scene_final_outputs_complete(state, args):
        _log_step_plan(
            "workflow",
            outputs=_scene_skip_required_outputs(state, args),
            message="Skipping whole scene because desired outputs exist",
            enabled=args.log_to_console,
            scene_basename=scene.primary_basename,
        )
        _mark_scene_complete_from_existing_output(state, args)
        return state

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

    saved_output_paths = _scene_saved_output_paths(state, args)
    if saved_output_paths:
        log(
            "Wrote scene outputs: " + ", ".join(saved_output_paths),
            enabled=args.log_to_console,
            step="workflow",
            scene_basename=scene.primary_basename,
        )
    _write_scene_report(state, args, scene_started_utc=scene_started_utc)
    if not args.keep_temp_dir:
        temp_cleanup_paths = _scene_temp_cleanup_paths(state, args)
        if temp_cleanup_paths:
            log(
                f"Deleting {len(temp_cleanup_paths)} temp files",
                enabled=args.log_to_console,
                step="workflow",
                scene_basename=scene.primary_basename,
            )
            _delete_files(temp_cleanup_paths)
    log("Scene complete", enabled=args.log_to_console, step="workflow", scene_basename=scene.primary_basename)
    return state


def _run_workflow(args: argparse.Namespace) -> int:
    """Run the full WorldView preprocessing workflow.
    Args:
        args: Parsed CLI arguments.
    Returns:
        Process exit code.
    """
    filter_basenames = _parse_filter_basenames(args.filter_basename)
    input_files = _collect_input_files(args.input_file_glob)
    if not input_files:
        raise ValueError("No files matched --input-file-glob.")

    scenes = load_worldview_scenes_from_tif_files(input_files, filter_basenames=filter_basenames)
    log(
        f"Discovered {len(input_files)} input files across {len(scenes)} scenes",
        enabled=args.log_to_console,
        step="workflow",
    )
    processed_states: List[SceneWorkflowState] = []
    worker_count = _resolve_concurrent_processing(args.concurrent_processing)
    if worker_count > 1 and len(scenes) > 1:
        log(
            f"Running per-scene processing with {min(worker_count, len(scenes))} processes",
            enabled=args.log_to_console,
            step="workflow",
        )
        ordered_results: List[SceneWorkflowState | None] = [None] * len(scenes)
        with ProcessPoolExecutor(max_workers=min(worker_count, len(scenes))) as executor:
            future_to_index = {
                executor.submit(_process_scene, scene, args): index
                for index, scene in enumerate(scenes)
            }
            for future in as_completed(future_to_index):
                index = future_to_index[future]
                ordered_results[index] = future.result()
        processed_states = [state for state in ordered_results if state is not None]
    else:
        for scene in scenes:
            processed_states.append(_process_scene(scene, args))

    if processed_states and args.run_radiometric_normalization:
        log(
            f"Preparing grouped radiometric normalization for {len([state for state in processed_states if state.current_files])} scene outputs",
            enabled=args.log_to_console,
            step="workflow",
        )
        radiometric_output = _run_radiometric_normalization_workflow(
            [state.current_files[0] for state in processed_states if state.current_files],
            args=args,
            reference_state=processed_states[0],
        )
        if radiometric_output:
            log(
                f"Wrote radiometric normalization output {radiometric_output}",
                enabled=args.log_to_console,
                step="workflow",
            )

    log("All processing complete", enabled=args.log_to_console, step="workflow")
    return 0


def _build_parser() -> argparse.ArgumentParser:
    """Build the WorldView CLI parser.
    Args:
        None.
    Returns:
        Configured argument parser.
    """
    parser = argparse.ArgumentParser(description="Run WorldView preprocessing on discovered tif scenes.")
    parser.add_argument("--config-yaml", help="Optional YAML config file.")
    parser.add_argument(
        "--input-file-glob",
        action="append",
        help="Glob used to find input files, for example '/data/**/*.TIF'.",
    )
    parser.add_argument("--input-dir", dest="input_file_glob", action="append", help=argparse.SUPPRESS)
    parser.add_argument(
        "--dem-file-path",
        default="online",
        help="DEM GeoTIFF path in WGS84 ellipsoidal height, or 'online' to download SRTM GL1 ellipsoidal to the temp dir.",
    )
    parser.add_argument("--dem-online-api-key")
    parser.add_argument("--dem-online-source", default=DEFAULT_OPENTOPOGRAPHY_DEMTYPE)
    parser.add_argument("--dem-online-api-endpoint", default=DEFAULT_OPENTOPOGRAPHY_GLOBALDEM_ENDPOINT)
    parser.add_argument("--dem-online-timeout-s", type=float, default=120.0)
    parser.add_argument("--envi-engine-path", help="Path to ENVI taskengine executable.")
    parser.add_argument("--atmospheric-method", choices=["flaash", "py6s", "none"], default="py6s")
    parser.add_argument("--epsg", type=int, default=4326)
    parser.add_argument("--nodata-value", type=float, default=-9999)
    parser.add_argument("--dtype", default="int16")
    parser.add_argument("--log-to-console", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--run-from-existing", action=argparse.BooleanOptionalAction, default=True)
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
    parser.add_argument("--py6s-auto-atmos-source", choices=["none", "nasa_power"], default="nasa_power")
    parser.add_argument("--py6s-auto-atmos-grid-size", type=int, default=3)
    parser.add_argument("--py6s-auto-atmos-search-days", type=int, default=1)
    parser.add_argument("--py6s-auto-atmos-timeout-s", type=float, default=30.0)
    parser.add_argument("--py6s-auto-atmos-power-endpoint", default="https://power.larc.nasa.gov/api/temporal/daily/point")
    parser.add_argument("--footprint-epsg", type=int, default=4326)
    parser.add_argument("--filter-basename", action="append")
    parser.add_argument("--max-cloud-cover-to-process", type=float)
    parser.add_argument("--output-dir")
    parser.add_argument("--fetch-atmosphere-output-suffix", default="_atmosphere")
    parser.add_argument("--atmospheric-correction-output-suffix", default="_atmospheric")
    parser.add_argument("--radiometric-normalization-output")
    parser.add_argument("--orthorectification-output-suffix", default="_ortho")
    parser.add_argument("--orthorectification-pan-output-suffix", default="_pan_ortho")
    parser.add_argument("--orthorectification-rpc-refinement-geojson")
    parser.add_argument("--pansharpen-output-suffix", default="_pansharpen")
    parser.add_argument("--skip-existing", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--concurrent-processing", default=1)
    parser.add_argument("--overview-scales", nargs="+")
    parser.add_argument("--temp-dir")
    parser.add_argument("--keep-temp-dir", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--scratch-dir", dest="temp_dir", help=argparse.SUPPRESS)
    parser.add_argument("--run-fetch-atmosphere", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--save-fetch-atmosphere", default="$temp")
    parser.add_argument("--run-atmospheric-correction", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--save-atmospheric-correction", default="$temp")
    parser.add_argument("--calculate-overviews-atmospheric-correction", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--run-orthorectification", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--save-orthorectification", default="$temp")
    parser.add_argument("--calculate-overviews-orthorectification", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--run-pansharpen", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--save-pansharpen", default="$temp")
    parser.add_argument("--calculate-overviews-pansharpen", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--run-cloud-mask", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--save-cloud-mask", default="$output")
    parser.add_argument("--calculate-overviews-cloud-mask", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--run-alignment", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--save-alignment", default="$temp")
    parser.add_argument("--calculate-overviews-alignment", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--run-radiometric-normalization", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--save-radiometric-normalization", default="$temp")
    parser.add_argument("--calculate-overviews-radiometric-normalization", action=argparse.BooleanOptionalAction, default=False)
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
    parser.add_argument("--group-by-basename")
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
    parser.add_argument("--alignment-band-index", type=int, default=0)
    parser.add_argument("--alignment-moving-band-index", type=int)
    parser.add_argument("--alignment-fixed-band-index", type=int)
    parser.add_argument("--alignment-moving-nodata", type=float)
    parser.add_argument("--alignment-fixed-nodata", type=float)
    parser.add_argument("--alignment-output-nodata", type=float)
    parser.add_argument("--alignment-min-valid-fraction", type=float, default=0.01)
    parser.add_argument("--alignment-split-factor", type=int, default=2)
    parser.add_argument("--alignment-clip-fixed-to-moving", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--alignment-output-on-moving-grid", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--alignment-trim-edge-invalid", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--alignment-edge-trim-depth", type=int, default=8)
    parser.add_argument("--alignment-edge-trim-detection-band-index", type=int, default=0)
    parser.add_argument("--alignment-edge-trim-invalid-below", type=float)
    parser.add_argument("--alignment-edge-trim-invalid-above", type=float)
    parser.add_argument("--alignment-enforce-mutual-valid-mask", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--alignment-use-edge-proxies", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--alignment-solve-resolution", type=float)
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """Parse CLI/config arguments, validate them, and execute the workflow."""
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument("--config-yaml")
    config_args, _ = config_parser.parse_known_args(argv)

    config_defaults: Dict = {}
    if config_args.config_yaml:
        config_defaults = _normalize_config_defaults(_load_worldview_yaml_config(config_args.config_yaml))

    parser = _build_parser()
    if config_defaults:
        parser.set_defaults(**config_defaults)
    args, unknown_args = parser.parse_known_args(argv)
    _apply_unknown_prefixed_args(args, unknown_args)

    if not args.input_file_glob:
        parser.error("--input-file-glob is required (via CLI or --config-yaml).")
    if not args.dem_file_path and not args.existing_mul_ortho_input:
        parser.error("--dem-file-path is required unless --existing-mul-ortho-input is provided.")
    if args.run_pansharpen and not args.run_orthorectification and not args.existing_mul_ortho_input:
        parser.error("--run-pansharpen requires orthorectified inputs, --run-orthorectification, or --existing-mul-ortho-input.")
    if args.run_pansharpen and args.existing_mul_ortho_input and not args.existing_pan_ortho_input:
        parser.error("--existing-pan-ortho-input is required when using --existing-mul-ortho-input with pansharpen.")
    if (
        args.orthorectification_rpc_refinement_geojson
        and not os.path.isfile(args.orthorectification_rpc_refinement_geojson)
    ):
        parser.error(
            "--orthorectification-rpc-refinement-geojson does not exist: "
            f"{args.orthorectification_rpc_refinement_geojson}"
        )
    if args.run_alignment and not args.alignment_fixed_image:
        parser.error("--alignment-fixed-image is required when --run-alignment is enabled.")
    if args.run_alignment and args.alignment_fixed_image and not os.path.isfile(args.alignment_fixed_image):
        parser.error(f"--alignment-fixed-image does not exist: {args.alignment_fixed_image}")
    if args.atmospheric_method == "flaash" and args.run_atmospheric_correction and not args.skip_flaash and not args.envi_engine_path:
        parser.error("--envi-engine-path is required when running FLAASH.")
    if args.run_cloud_mask and args.cloud_mask_method and args.cloud_mask_command:
        parser.error("Use either --cloud-mask-method or --cloud-mask-command, not both.")
    if args.cloud_mask_inference_resolution_m <= 0:
        parser.error("--cloud-mask-inference-resolution-m must be > 0.")
    if args.fetch_atmosphere_grid_size < 1:
        parser.error("--fetch-atmosphere-grid-size must be >= 1.")
    if args.fetch_atmosphere_search_days < 0:
        parser.error("--fetch-atmosphere-search-days must be >= 0.")
    if args.fetch_atmosphere_timeout_s <= 0:
        parser.error("--fetch-atmosphere-timeout-s must be > 0.")
    if args.dem_online_timeout_s <= 0:
        parser.error("--dem-online-timeout-s must be > 0.")
    for arg_name in (
        "save_fetch_atmosphere",
        "save_atmospheric_correction",
        "save_orthorectification",
        "save_pansharpen",
        "save_cloud_mask",
        "save_alignment",
        "save_radiometric_normalization",
    ):
        save_value = getattr(args, arg_name)
        if save_value in (None, ""):
            continue
        normalized = str(save_value).strip()
        if normalized in ("$temp", "$output"):
            continue
        if normalized.startswith(("$temp/", "$output/", "./")):
            continue
        if os.path.isabs(normalized):
            continue
        if normalized.startswith("../"):
            continue
        if normalized in ("temp", "output") or normalized.startswith(("temp/", "output/")):
            parser.error(
                f"--{arg_name.replace('_', '-')} no longer supports legacy 'temp'/'output' prefixes; "
                "use '$temp' or '$output' instead."
            )
    if (
        args.max_cloud_cover_to_process is not None
        and (args.max_cloud_cover_to_process < 0 or args.max_cloud_cover_to_process > 100)
    ):
        parser.error("--max-cloud-cover-to-process must be in [0, 100].")
    try:
        args.concurrent_processing = _resolve_concurrent_processing(args.concurrent_processing)
    except ValueError as exc:
        parser.error(str(exc))
    if args.overview_scales is not None:
        if isinstance(args.overview_scales, str):
            args.overview_scales = [int(value.strip()) for value in args.overview_scales.split(",") if value.strip()]
        else:
            args.overview_scales = [int(value) for value in args.overview_scales]
    if (
        args.calculate_overviews_atmospheric_correction
        or args.calculate_overviews_orthorectification
        or args.calculate_overviews_pansharpen
        or args.calculate_overviews_cloud_mask
        or args.calculate_overviews_alignment
        or args.calculate_overviews_radiometric_normalization
    ) and not args.overview_scales:
        parser.error("--overview-scales is required when any calculate-overviews-* option is enabled.")
    if args.alignment_band_index < 0:
        parser.error("--alignment-band-index must be >= 0.")
    if args.alignment_moving_band_index is not None and args.alignment_moving_band_index < 0:
        parser.error("--alignment-moving-band-index must be >= 0.")
    if args.alignment_fixed_band_index is not None and args.alignment_fixed_band_index < 0:
        parser.error("--alignment-fixed-band-index must be >= 0.")
    if args.alignment_min_valid_fraction <= 0 or args.alignment_min_valid_fraction > 1:
        parser.error("--alignment-min-valid-fraction must be in (0, 1].")
    if args.alignment_split_factor < 0:
        parser.error("--alignment-split-factor must be >= 0.")
    if args.alignment_edge_trim_depth <= 0:
        parser.error("--alignment-edge-trim-depth must be > 0.")
    if args.alignment_edge_trim_detection_band_index < 0:
        parser.error("--alignment-edge-trim-detection-band-index must be >= 0.")
    if args.alignment_solve_resolution is not None and args.alignment_solve_resolution <= 0:
        parser.error("--alignment-solve-resolution must be > 0.")
    try:
        _parse_json_dict(args.cloud_mask_omnicloud_kwargs_json)
    except Exception as exc:
        parser.error(f"Invalid --cloud-mask-omnicloud-kwargs-json: {exc}")
    try:
        _parse_json_dict(args.radiometric_normalization_kwargs_json)
    except Exception as exc:
        parser.error(f"Invalid --radiometric-normalization-kwargs-json: {exc}")
    try:
        _normalize_group_by_basename_spec(args.group_by_basename)
    except Exception as exc:
        parser.error(f"Invalid --group-by-basename: {exc}")
    return _run_workflow(args)


if __name__ == "__main__":
    sys.exit(main())
