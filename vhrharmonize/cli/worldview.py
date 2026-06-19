#!/usr/bin/env python3
"""CLI wrapper for the WorldView preprocessing workflow."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
import os
import re
import shutil
import struct
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from typing import Any, Dict, List, Mapping, Optional

import geopandas as gpd
from osgeo import gdal
from wcmatch import fnmatch as wc_fnmatch
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
from vhrharmonize.preprocess.seamline_metadata import write_seamline_metadata_gpkg
from vhrharmonize.providers.worldview import (
    WorldViewImage,
    WorldViewScene,
    load_worldview_scenes_from_tif_files,
)

RASTER_STEP_ORDER = [
    "file_source",
    "atmospheric_correction",
    "orthorectification",
    "pansharpen",
    "cloud_mask",
    "alignment",
]

WCMATCH_INPUT_FLAGS = (
    glob.GLOBSTAR
    | glob.BRACE
    | glob.EXTGLOB
    | glob.GLOBTILDE
    | glob.GLOBSTARLONG
    | glob.NEGATE
)
WCMATCH_GROUP_FLAGS = (
    getattr(wc_fnmatch, "BRACE", 0)
    | getattr(wc_fnmatch, "EXTMATCH", 0)
    | getattr(wc_fnmatch, "EXTGLOB", 0)
    | getattr(wc_fnmatch, "NEGATE", 0)
    | getattr(wc_fnmatch, "GLOBSTAR", 0)
)

@dataclass
class SceneWorkflowState:
    """Workflow state for a single discovered WorldView scene."""

    scene: WorldViewScene
    step_dirs: Dict[str, str]
    scene_bbox_wgs84: tuple[float, float, float, float]
    current_files: List[str]
    current_step: str = "file_source"
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
    if step_name in {"file_source", "raw"}:
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


def _radiometric_steps_include(args: argparse.Namespace, step_name: str) -> bool:
    """Return whether the SpectralMatch steps list includes a step."""
    steps = getattr(args, "match_steps", None)
    if steps is None:
        radiometric_kwargs = _parse_json_dict(getattr(args, "radiometric_normalization_kwargs_json", None))
        steps = radiometric_kwargs.get("steps")
    if steps is None:
        return False
    if isinstance(steps, str):
        return steps == step_name
    if isinstance(steps, (list, tuple)):
        return step_name in steps
    return False


def _explicit_cli_arg_present(argv: List[str], arg_name: str) -> bool:
    """Return whether a long CLI option was explicitly provided."""
    option = f"--{arg_name.replace('_', '-')}"
    return any(token == option or token.startswith(f"{option}=") for token in argv)


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
        matches = glob.glob(pattern, flags=WCMATCH_INPUT_FLAGS)
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
        if not re.fullmatch(r"[-+]?\d+", normalized):
            raise ValueError("concurrent_processing must be an integer or 'num_cpu'.")
        resolved = int(normalized)
    else:
        if not isinstance(value, int):
            raise ValueError("concurrent_processing must be an integer or 'num_cpu'.")
        resolved = int(value)
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


def _classify_save_target(save_value: Optional[str], *, default: str) -> tuple[str, str]:
    """Return a normalized save target and target kind."""
    normalized = (save_value if save_value not in (None, "") else default).strip()
    if normalized == "$temp":
        return normalized, "temp_root"
    if normalized.startswith("$temp/"):
        return normalized, "temp_child"
    if normalized == "$output":
        return normalized, "output_root"
    if normalized.startswith("$output/"):
        return normalized, "output_child"
    if normalized.startswith("./"):
        return normalized, "input_relative"
    if os.path.isabs(normalized):
        return normalized, "absolute"
    return normalized, "cwd_relative"


def _format_allowed_save_targets(allowed_modes: set[str]) -> str:
    """Return a human-readable list of accepted save target forms."""
    labels = [
        ("temp_root", "$temp"),
        ("temp_child", "$temp/..."),
        ("output_root", "$output"),
        ("output_child", "$output/..."),
        ("input_relative", "./relative/to/input"),
        ("absolute", "/custom/absolute"),
        ("cwd_relative", "relative/to/pwd"),
    ]
    return " | ".join(label for mode, label in labels if mode in allowed_modes)


def _resolve_save_target(
    save_value: Optional[str],
    *,
    default: str,
    temp_root: str,
    output_root: str,
    relative_base_folder: str,
    accepted_modes: set[str],
) -> str:
    """Resolve a save target after checking that its target kind is accepted."""
    save_mode, save_kind = _classify_save_target(save_value, default=default)
    if save_kind not in accepted_modes:
        raise ValueError(
            f"Unsupported save target '{save_mode}'. Accepted forms: {_format_allowed_save_targets(accepted_modes)}"
        )
    if save_kind == "temp_root":
        return temp_root
    if save_kind == "temp_child":
        return os.path.join(temp_root, save_mode[len("$temp/"):])
    if save_kind == "output_root":
        return output_root
    if save_kind == "output_child":
        return os.path.join(output_root, save_mode[len("$output/"):])
    if save_kind == "input_relative":
        return resolve_relative_to_input(save_mode, relative_base_folder)
    if save_kind == "absolute":
        return save_mode
    return os.path.abspath(save_mode)


def _resolve_step_save_dir(
    save_value: Optional[str],
    *,
    temp_root: str,
    output_root: str,
    relative_base_folder: str,
    accepted_modes: set[str] | None = None,
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
    resolved_dir = _resolve_save_target(
        save_value,
        default="$temp",
        temp_root=temp_root,
        output_root=output_root,
        relative_base_folder=relative_base_folder,
        accepted_modes=accepted_modes
        or {"temp_root", "temp_child", "output_root", "output_child", "input_relative", "absolute", "cwd_relative"},
    )
    os.makedirs(resolved_dir, exist_ok=True)
    return resolved_dir


def _resolve_single_output_save_path(
    save_value: Optional[str],
    *,
    default: str,
    temp_root: str,
    output_root: str,
    relative_base_folder: str,
    accepted_modes: set[str] | None = None,
) -> str:
    """Resolve a configured single-output save path."""
    resolved_path = _resolve_save_target(
        save_value,
        default=default,
        temp_root=temp_root,
        output_root=output_root,
        relative_base_folder=relative_base_folder,
        accepted_modes=accepted_modes or {"temp_child", "absolute", "cwd_relative"},
    )
    os.makedirs(os.path.dirname(resolved_path) or ".", exist_ok=True)
    return resolved_path


def _validate_save_target_value(
    save_value: Optional[str],
    *,
    arg_name: str,
    default: str,
    accepted_modes: set[str],
) -> None:
    """Validate a save target before workflow execution."""
    normalized, save_kind = _classify_save_target(save_value, default=default)
    if save_kind not in accepted_modes:
        raise ValueError(
            f"--{arg_name.replace('_', '-')} does not support '{normalized}'. "
            f"Accepted forms: {_format_allowed_save_targets(accepted_modes)}"
        )


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
    return "file_source"


def _enabled_raster_steps(args: argparse.Namespace) -> List[str]:
    """Return enabled raster processing steps in workflow order."""
    enabled_steps: List[str] = []
    if args.run_file_source:
        enabled_steps.append("file_source")
    if args.run_atmospheric_correction:
        enabled_steps.append("atmospheric_correction")
    if args.run_orthorectification:
        enabled_steps.append("orthorectification")
    if args.run_pansharpen:
        enabled_steps.append("pansharpen")
    if args.run_cloud_mask:
        enabled_steps.append("cloud_mask")
    if args.run_alignment:
        enabled_steps.append("alignment")
    return enabled_steps


def _step_outputs_exist(output_paths: List[str]) -> bool:
    """Return whether all output paths exist.
    Args:
        output_paths: Output file paths to check.
    Returns:
        True when every output path exists.
    """
    return bool(output_paths) and all(os.path.exists(path) for path in output_paths)


def _is_gdal_raster_path(path: str) -> bool:
    """Return whether a path should be validated as a GDAL raster."""
    extension = os.path.splitext(path)[1].lower()
    return extension in {".tif", ".tiff", ".dat", ".img", ".vrt"}


def _gdal_raster_is_valid(path: str, *, validity_check_grid_size: int = 0) -> tuple[bool, str | None]:
    """Check that a raster can be opened and read by GDAL.
    Args:
        path: Raster path to validate.
        validity_check_grid_size: Pixel sampling grid size. 0 disables pixel validity sampling.
    Returns:
        Tuple of validity and optional reason.
    """
    if not os.path.exists(path):
        return False, "missing"
    dataset = gdal.OpenEx(path, gdal.OF_RASTER)
    if dataset is None:
        return False, "GDAL open failed"
    try:
        band_count = dataset.RasterCount
        width = dataset.RasterXSize
        height = dataset.RasterYSize
        if band_count < 1 or width < 1 or height < 1:
            return False, f"invalid raster shape bands={band_count} size={width}x{height}"

        read_w = min(width, 256)
        read_h = min(height, 256)
        offsets = [
            (0, 0),
            (max(0, width - read_w), max(0, height - read_h)),
        ]
        for band_index in range(1, band_count + 1):
            band = dataset.GetRasterBand(band_index)
            if band is None:
                return False, f"missing band {band_index}"
            for xoff, yoff in offsets:
                data = band.ReadRaster(xoff, yoff, read_w, read_h)
                if data is None:
                    return False, f"GDAL ReadRaster failed for band {band_index}"
        if validity_check_grid_size > 0 and not _gdal_raster_has_valid_sample(
            dataset,
            width,
            height,
            validity_check_grid_size,
        ):
            return False, "no finite valid pixels found in validity sample"
    except Exception as exc:
        return False, str(exc)
    finally:
        dataset = None
    return True, None


def _gdal_raster_has_valid_sample(dataset: gdal.Dataset, width: int, height: int, grid_size: int) -> bool:
    """Return whether any sampled raster pixel is finite and non-nodata."""
    if grid_size < 1:
        return True
    block_size = 512 if grid_size == 1 else grid_size
    for band_index in range(1, dataset.RasterCount + 1):
        band = dataset.GetRasterBand(band_index)
        if band is None:
            continue
        for yoff in range(0, height, block_size):
            read_h = min(block_size, height - yoff)
            sample_y = 0 if grid_size == 1 else min(grid_size // 2, read_h - 1)
            for xoff in range(0, width, block_size):
                read_w = min(block_size, width - xoff)
                sample_x = 0 if grid_size == 1 else min(grid_size // 2, read_w - 1)
                window_xoff = xoff if grid_size == 1 else xoff + sample_x
                window_yoff = yoff if grid_size == 1 else yoff + sample_y
                window_w = read_w if grid_size == 1 else 1
                window_h = read_h if grid_size == 1 else 1
                if _gdal_band_window_has_valid_pixel(
                    band,
                    window_xoff,
                    window_yoff,
                    window_w,
                    window_h,
                ):
                    return True
    return False


def _gdal_band_window_has_valid_pixel(
    band: gdal.Band,
    xoff: int,
    yoff: int,
    read_w: int,
    read_h: int,
) -> bool:
    """Return whether a sampled band window has a finite non-nodata pixel."""
    struct_format = _gdal_data_type_struct_format(band.DataType)
    if struct_format is None:
        return True
    data = band.ReadRaster(xoff, yoff, read_w, read_h)
    if data is None:
        return False
    mask = None
    mask_band = band.GetMaskBand()
    if mask_band is not None:
        mask = mask_band.ReadRaster(xoff, yoff, read_w, read_h, buf_type=gdal.GDT_Byte)

    nodata = band.GetNoDataValue()
    item_size = struct.calcsize(struct_format)
    if read_w == 1 and read_h == 1:
        return _sample_pixel_is_valid(
            struct.unpack_from(struct_format, data, 0)[0],
            nodata,
            mask[0] if mask else 255,
        )

    for pixel_index in range(0, read_w * read_h):
        mask_value = mask[pixel_index] if mask else 255
        value = struct.unpack_from(struct_format, data, pixel_index * item_size)[0]
        if _sample_pixel_is_valid(value, nodata, mask_value):
            return True
    return False


def _gdal_data_type_struct_format(data_type: int) -> str | None:
    """Return a struct format for scalar GDAL data types."""
    formats = {
        gdal.GDT_Byte: "=B",
        gdal.GDT_UInt16: "=H",
        gdal.GDT_Int16: "=h",
        gdal.GDT_UInt32: "=I",
        gdal.GDT_Int32: "=i",
        gdal.GDT_Float32: "=f",
        gdal.GDT_Float64: "=d",
    }
    if hasattr(gdal, "GDT_Int8"):
        formats[getattr(gdal, "GDT_Int8")] = "=b"
    if hasattr(gdal, "GDT_UInt64"):
        formats[getattr(gdal, "GDT_UInt64")] = "=Q"
    if hasattr(gdal, "GDT_Int64"):
        formats[getattr(gdal, "GDT_Int64")] = "=q"
    return formats.get(data_type)


def _sample_pixel_is_valid(value: object, nodata: float | int | None, mask_value: int = 255) -> bool:
    """Return whether a scalar sample is finite and not nodata."""
    if mask_value == 0:
        return False
    try:
        sample = float(value)
    except (TypeError, ValueError):
        return False
    if not math.isfinite(sample):
        return False
    if nodata is None:
        return True
    if isinstance(nodata, float) and math.isnan(nodata):
        return True
    return sample != nodata


def _existing_outputs_are_reusable(
    output_paths: List[str],
    *,
    check_validity: bool,
    validity_check_grid_size: int,
    log_to_console: bool,
    step: str,
    scene_basename: str | None = None,
) -> bool:
    """Return whether existing outputs can be reused.
    Args:
        output_paths: Output file paths expected for reuse.
        check_validity: Whether GDAL-readable raster outputs must be valid.
        log_to_console: Whether to log validation failures.
        step: Step name used for concise logging.
        scene_basename: Optional scene basename for the log prefix.
    Returns:
        True when all outputs exist and enabled validity checks pass.
    """
    if not _step_outputs_exist(output_paths):
        return False
    if not check_validity:
        return True
    for output_path in output_paths:
        if not _is_gdal_raster_path(output_path):
            continue
        is_valid, reason = _gdal_raster_is_valid(
            output_path,
            validity_check_grid_size=validity_check_grid_size,
        )
        if not is_valid:
            _log_step_plan(
                step,
                outputs=[output_path],
                message=f"Existing output invalid; rerunning ({reason})",
                enabled=log_to_console,
                scene_basename=scene_basename,
            )
            return False
    return True


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

    if args.run_file_source:
        step_dirs["file_source"] = _resolve_step_save_dir(
            args.save_file_source,
            temp_root=resolved_temp_root,
            output_root=resolved_output_root,
            relative_base_folder=relative_output_base,
        )
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
    if args.run_seamline_metadata:
        step_dirs["seamline_metadata"] = _resolve_single_output_save_path(
            args.save_seamline_metadata,
            default="$temp/seamline_metadata.gpkg",
            temp_root=resolved_temp_root,
            output_root=resolved_output_root,
            relative_base_folder=relative_output_base,
        )
    if args.run_radiometric_normalization:
        step_dirs["radiometric_normalization"] = _resolve_single_output_save_path(
            args.save_radiometric_normalization,
            default="$temp/radiometric_root.tif",
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
    pan_image = _require_scene_image(state.scene, "pan")
    if args.run_file_source:
        file_source_map: Dict[str, str] = {}
        file_source_map.update(_image_source_file_map(mul_image, state.step_dirs["file_source"]))
        file_source_map.update(_image_source_file_map(pan_image, state.step_dirs["file_source"]))
        file_source_outputs = list(file_source_map.values())
        current_mul_outputs = [file_source_map.get(os.path.abspath(mul_image.tif_file), mul_image.tif_file)]
    else:
        file_source_outputs = [mul_image.tif_file]
        current_mul_outputs = [mul_image.tif_file]

    expected_outputs: Dict[str, List[str]] = {
        "file_source": file_source_outputs,
        "raw": [mul_image.tif_file],
    }

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


def _scene_step_output_keys(step_name: str) -> List[str]:
    """Return expected-output keys that belong to one workflow step."""
    if step_name == "file_source":
        return ["file_source"]
    if step_name == "orthorectification":
        return ["orthorectification", "orthorectification_pan"]
    if step_name == "cloud_mask":
        return ["cloud_mask", "cloud_mask_mask"]
    return [step_name]


def _scene_step_expected_outputs(
    expected_outputs: Dict[str, List[str]],
    step_name: str,
) -> List[str]:
    """Return expected output paths for a workflow step."""
    outputs: List[str] = []
    for key in _scene_step_output_keys(step_name):
        outputs.extend(expected_outputs.get(key, []))
    return outputs


def _get_scene_upload_source_step(
    state: SceneWorkflowState,
    args: argparse.Namespace,
) -> str:
    """Return the latest enabled processed raster step, or file_source."""
    expected_outputs = _get_expected_scene_step_outputs(state, args)
    for step_name in reversed(_enabled_raster_steps(args)):
        if step_name == "file_source":
            continue
        step_outputs = _scene_step_expected_outputs(expected_outputs, step_name)
        if _existing_outputs_are_reusable(
            step_outputs,
            check_validity=args.run_from_existing_check_validity,
            validity_check_grid_size=args.validity_check_grid_size,
            log_to_console=False,
            step=step_name,
            scene_basename=state.scene.primary_basename,
        ):
            return step_name
    return "file_source"


def _get_scene_upload_source_files(
    state: SceneWorkflowState,
    args: argparse.Namespace,
) -> List[str]:
    """Return files from the latest enabled processed raster step, or file_source."""
    source_step = _get_scene_upload_source_step(state, args)
    expected_outputs = _get_expected_scene_step_outputs(state, args)
    return _scene_step_expected_outputs(expected_outputs, source_step)


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
        current_files=[_get_worldview_scene_step_path(scene, "mul", "file_source")],
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
    if args.run_file_source and not _is_temp_save_value(args.save_file_source):
        required_outputs.extend(expected_outputs.get("file_source", []))
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
    if raw_spec in (None, ""):
        return None
    if isinstance(raw_spec, str):
        stripped = raw_spec.strip()
        if stripped.startswith("{"):
            return _normalize_group_by_basename_spec(json.loads(stripped))
        raise ValueError("group_by_basename must be a JSON object such as {'name.tif': ['auto:*123*', 'file:/tmp/ref.tif']}.")
    return _validate_radiometric_group_spec(raw_spec)


def _validate_radiometric_group_spec(group_spec: object) -> Dict[str, object]:
    """Validate the named radiometric grouping JSON shape."""
    if not isinstance(group_spec, dict) or not group_spec:
        raise ValueError("group_by_basename must be a non-empty object with output filename keys.")
    normalized: Dict[str, object] = {}
    for output_name, value in group_spec.items():
        if not isinstance(output_name, str) or not output_name:
            raise ValueError("Each radiometric group key must be a non-empty output filename string.")
        normalized[output_name] = _validate_radiometric_group_value(value)
    return normalized


def _validate_radiometric_group_value(value: object) -> object:
    """Validate a radiometric group value."""
    if isinstance(value, str):
        if not (value.startswith("auto:") or value.startswith("file:")):
            raise ValueError("Radiometric group strings must start with 'auto:' or 'file:'.")
        if value in {"auto:", "file:"}:
            raise ValueError("Radiometric group strings must include a pattern or path after the prefix.")
        return value
    if isinstance(value, list):
        if not value:
            raise ValueError("Radiometric group lists must not be empty.")
        return [_validate_radiometric_group_item(item) for item in value]
    raise ValueError("Radiometric group values must be strings or lists.")


def _validate_radiometric_group_item(item: object) -> object:
    """Validate an item inside a radiometric group list."""
    if isinstance(item, str):
        return _validate_radiometric_group_value(item)
    if isinstance(item, dict):
        return _validate_radiometric_group_spec(item)
    raise ValueError("Radiometric group list items must be prefixed strings or nested group objects.")


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
        if wc_fnmatch.fnmatch(os.path.basename(path), pattern, flags=WCMATCH_GROUP_FLAGS)
    ]
    if not matches:
        raise ValueError(f"No radiometric normalization inputs matched pattern: {pattern}")
    return matches


def _resolve_radiometric_input_token(token: str, available_paths: List[str]) -> List[str]:
    """Resolve a prefixed radiometric group token into input paths."""
    source, value = token.split(":", 1)
    if source == "auto":
        return _match_radiometric_input_patterns(value, available_paths)
    if source == "file":
        return [value]
    raise ValueError("Radiometric group strings must start with 'auto:' or 'file:'.")


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
    *,
    output_name: str,
    temp_root: str,
    output_root: str,
) -> str:
    """Resolve a radiometric group output path.
    Args:
        output_name: Group output filename.
        temp_root: Temp root directory.
        output_root: Output root directory.
    Returns:
        Resolved radiometric group output path.
    """
    del temp_root
    group_output_path = os.path.join(output_root, os.path.basename(output_name))
    os.makedirs(os.path.dirname(group_output_path) or ".", exist_ok=True)
    return group_output_path


def _run_named_radiometric_group(
    output_name: str,
    group_value: object,
    *,
    available_paths: List[str],
    args: argparse.Namespace,
    temp_root: str,
    output_root: str,
) -> str:
    """Run one named radiometric group and return its output path."""
    child_inputs: List[str] = []
    if isinstance(group_value, str):
        child_inputs.extend(_resolve_radiometric_input_token(group_value, available_paths))
    elif isinstance(group_value, list):
        for item in group_value:
            if isinstance(item, str):
                child_inputs.extend(_resolve_radiometric_input_token(item, available_paths))
            elif isinstance(item, dict):
                for child_output_name, child_value in item.items():
                    child_inputs.append(
                        _run_named_radiometric_group(
                            child_output_name,
                            child_value,
                            available_paths=available_paths,
                            args=args,
                            temp_root=temp_root,
                            output_root=output_root,
                        )
                    )
            else:
                raise ValueError("Radiometric group list items must be prefixed strings or nested group objects.")
    else:
        raise ValueError("Radiometric group values must be strings or lists.")

    child_inputs = _dedupe_paths(child_inputs)
    radiometric_kwargs = _build_radiometric_kwargs(args)
    radiometric_kwargs.pop("shared_output_image_path", None)
    group_output_path = _resolve_radiometric_group_output_path(
        output_name=output_name,
        temp_root=temp_root,
        output_root=output_root,
    )
    if args.run_from_existing and _existing_outputs_are_reusable(
        [group_output_path],
        check_validity=args.run_from_existing_check_validity,
        validity_check_grid_size=args.validity_check_grid_size,
        log_to_console=args.log_to_console,
        step="radiometric",
    ):
        _log_step_plan(
            "radiometric",
            outputs=[group_output_path],
            message=f"Skipping group {output_name} because output exists",
            enabled=args.log_to_console,
        )
        return group_output_path

    radiometric_kwargs.setdefault("shared_input_images", child_inputs)
    radiometric_kwargs.setdefault("shared_output_image_path", group_output_path)
    radiometric_kwargs.setdefault("shared_temp_dir", os.path.join(temp_root, "spectralmatch"))
    radiometric_kwargs.setdefault("delete_temp_dir", args.keep_temp_dir is not True)
    radiometric_kwargs.setdefault("shared_debug_logs", args.log_to_console)
    radiometric_kwargs.setdefault("shared_output_dtype", args.dtype)
    _log_step_plan(
        "radiometric",
        inputs=child_inputs,
        outputs=[group_output_path],
        message=f"Running SpectralMatch group {output_name}",
        enabled=args.log_to_console,
    )
    radiometric_normalization(
        method=args.radiometric_normalization_method,
        log_to_console=args.log_to_console,
        **radiometric_kwargs,
    )
    if args.calculate_overviews_radiometric_normalization:
        log("Calculating overviews for step radiometric_normalization", enabled=args.log_to_console, step="overviews")
        calculate_raster_overviews(group_output_path, args.overview_scales)
    return group_output_path


def _run_named_radiometric_groups(
    group_spec: Dict[str, object],
    *,
    available_paths: List[str],
    args: argparse.Namespace,
    temp_root: str,
    output_root: str,
) -> Optional[str]:
    """Run named radiometric groups in order and return the final output path."""
    final_output: Optional[str] = None
    for output_name, group_value in group_spec.items():
        final_output = _run_named_radiometric_group(
            output_name,
            group_value,
            available_paths=available_paths,
            args=args,
            temp_root=temp_root,
            output_root=output_root,
        )
    return final_output


def _run_default_radiometric_normalization(
    available_paths: List[str],
    *,
    args: argparse.Namespace,
    output_path: str,
    temp_root: str,
) -> str:
    """Run one default radiometric normalization over all scene outputs."""
    radiometric_kwargs = _build_radiometric_kwargs(args)
    group_output_path = str(radiometric_kwargs.get("shared_output_image_path") or output_path)
    if args.run_from_existing and _existing_outputs_are_reusable(
        [group_output_path],
        check_validity=args.run_from_existing_check_validity,
        validity_check_grid_size=args.validity_check_grid_size,
        log_to_console=args.log_to_console,
        step="radiometric",
    ):
        _log_step_plan(
            "radiometric",
            outputs=[group_output_path],
            message="Skipping radiometric normalization because output exists",
            enabled=args.log_to_console,
        )
        return group_output_path

    radiometric_kwargs.setdefault("shared_input_images", available_paths)
    radiometric_kwargs.setdefault("shared_output_image_path", group_output_path)
    radiometric_kwargs.setdefault("shared_temp_dir", os.path.join(temp_root, "spectralmatch"))
    radiometric_kwargs.setdefault("delete_temp_dir", args.keep_temp_dir is not True)
    radiometric_kwargs.setdefault("shared_debug_logs", args.log_to_console)
    radiometric_kwargs.setdefault("shared_output_dtype", args.dtype)
    _log_step_plan(
        "radiometric",
        inputs=available_paths,
        outputs=[group_output_path],
        message="Running SpectralMatch radiometric normalization",
        enabled=args.log_to_console,
    )
    radiometric_normalization(
        method=args.radiometric_normalization_method,
        log_to_console=args.log_to_console,
        **radiometric_kwargs,
    )
    if args.calculate_overviews_radiometric_normalization:
        log("Calculating overviews for step radiometric_normalization", enabled=args.log_to_console, step="overviews")
        calculate_raster_overviews(group_output_path, args.overview_scales)
    return group_output_path


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
    if group_spec is not None:
        return _run_named_radiometric_groups(
            group_spec,
            available_paths=available_paths,
            args=args,
            temp_root=reference_state.step_dirs["temp_root"],
            output_root=reference_state.step_dirs["output_root"],
        )

    return _run_default_radiometric_normalization(
        available_paths,
        args=args,
        output_path=reference_state.step_dirs["radiometric_normalization"],
        temp_root=reference_state.step_dirs["temp_root"],
    )


def _run_seamline_metadata_workflow(
    states: List[SceneWorkflowState],
    *,
    args: argparse.Namespace,
    reference_state: SceneWorkflowState,
) -> str | None:
    """Run the aggregate seamline metadata GeoPackage step.
    Args:
        states: Processed scene states.
        args: Parsed CLI arguments.
        reference_state: State used to resolve aggregate output directories.
    Returns:
        Seamline metadata GeoPackage path or None.
    """
    if not args.run_seamline_metadata:
        return None

    output_path = reference_state.step_dirs["seamline_metadata"]
    if args.run_from_existing and _existing_outputs_are_reusable(
        [output_path],
        check_validity=False,
        validity_check_grid_size=args.validity_check_grid_size,
        log_to_console=args.log_to_console,
        step="seamline_metadata",
    ):
        _log_step_plan(
            "seamline_metadata",
            outputs=[output_path],
            message="Skipping because output exists",
            enabled=args.log_to_console,
        )
        return output_path

    _log_step_plan(
        "seamline_metadata",
        inputs=[state.current_files[0] for state in states if state.current_files],
        outputs=[output_path],
        message="Writing footprint metadata GeoPackage",
        enabled=args.log_to_console,
    )
    return write_seamline_metadata_gpkg(
        states,
        output_path,
        layer=args.seamline_metadata_layer,
        image_field_name=args.seamline_metadata_image_field_name,
        footprint_source=args.seamline_metadata_footprint_source,
        calculate_bounds_eight_connected=args.seamline_metadata_calculate_bounds_eight_connected,
        epsg=args.epsg,
    )


def _apply_weighted_seamline_metadata_defaults(args: argparse.Namespace, seamline_metadata_output: str | None) -> None:
    """Default weighted seamline inputs from the generated WorldView metadata GPKG."""
    if not seamline_metadata_output:
        return
    if not _radiometric_steps_include(args, "weighted_seamline"):
        return
    if not getattr(args, "match_weighted_seamline_input_polygons", None):
        setattr(args, "match_weighted_seamline_input_polygons", seamline_metadata_output)
    if not getattr(args, "match_weighted_seamline_input_layer", None):
        setattr(args, "match_weighted_seamline_input_layer", args.seamline_metadata_layer)
    if not getattr(args, "match_weighted_seamline_image_field_name", None):
        setattr(args, "match_weighted_seamline_image_field_name", args.seamline_metadata_image_field_name)


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


def _image_source_file_map(image: Optional[WorldViewImage], output_dir: str) -> Dict[str, str]:
    """Return source bundle files mapped to their staged output paths."""
    if image is None:
        return {}
    source_paths = glob.glob(os.path.join(os.path.dirname(image.tif_file), f"{image.basename}.*"))
    if image.shp_file:
        source_paths.extend(glob.glob(f"{os.path.splitext(image.shp_file)[0]}.*"))
    path_map: Dict[str, str] = {}
    for source_path in _dedupe_paths([path for path in source_paths if os.path.isfile(path)]):
        path_map[os.path.abspath(source_path)] = os.path.join(output_dir, os.path.basename(source_path))
    return path_map


def _set_image_file_paths_from_source_map(image: Optional[WorldViewImage], path_map: Mapping[str, str]) -> None:
    """Update a WorldView image to point at staged source bundle files."""
    if image is None:
        return
    for attr_name in ("tif_file", "imd_file", "shp_file", "til_file"):
        current_path = getattr(image, attr_name)
        if current_path is None:
            continue
        staged_path = path_map.get(os.path.abspath(current_path))
        if staged_path:
            setattr(image, attr_name, staged_path)


def _copy_file_source_bundle(path_map: Mapping[str, str]) -> None:
    """Copy source bundle files to their staged file_source paths."""
    for source_path, output_path in path_map.items():
        if os.path.abspath(source_path) == os.path.abspath(output_path):
            continue
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        shutil.copy2(source_path, output_path)


def _run_file_source_step(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    """Run the file_source staging step."""
    if not args.run_file_source:
        return state.current_files
    output_dir = state.step_dirs["file_source"]
    path_map: Dict[str, str] = {}
    path_map.update(_image_source_file_map(state.scene.mul_image, output_dir))
    path_map.update(_image_source_file_map(state.scene.pan_image, output_dir))
    output_paths = list(path_map.values())
    expected_outputs = _get_expected_scene_step_outputs(state, args)

    if args.run_from_existing and _existing_outputs_are_reusable(
        output_paths,
        check_validity=args.run_from_existing_check_validity,
        validity_check_grid_size=args.validity_check_grid_size,
        log_to_console=args.log_to_console,
        step="file_source",
        scene_basename=state.scene.primary_basename,
    ):
        _log_step_plan(
            "file_source",
            outputs=output_paths,
            message="Skipping because output exists",
            enabled=args.log_to_console,
            scene_basename=state.scene.primary_basename,
        )
    else:
        _log_step_plan(
            "file_source",
            inputs=list(path_map.keys()),
            outputs=output_paths,
            message="Staging source bundle files",
            enabled=args.log_to_console,
            scene_basename=state.scene.primary_basename,
        )
        _copy_file_source_bundle(path_map)

    _set_image_file_paths_from_source_map(state.scene.mul_image, path_map)
    _set_image_file_paths_from_source_map(state.scene.pan_image, path_map)
    state.scene.step_outputs["file_source"] = output_paths
    mul_output = path_map.get(
        os.path.abspath(state.scene.mul_image.tif_file),
        state.scene.mul_image.tif_file,
    ) if state.scene.mul_image is not None else expected_outputs["file_source"][0]
    pan_output = path_map.get(
        os.path.abspath(state.scene.pan_image.tif_file),
        state.scene.pan_image.tif_file,
    ) if state.scene.pan_image is not None else None
    state.current_files = [mul_output]
    _set_worldview_scene_step_path(state.scene, "mul", "file_source", mul_output)
    if pan_output:
        _set_worldview_scene_step_path(state.scene, "pan", "file_source", pan_output)
    state.current_step = "file_source"
    if args.calculate_overviews_file_source:
        log(
            "Calculating overviews for step file_source",
            enabled=args.log_to_console,
            step="overviews",
            scene_basename=state.scene.primary_basename,
        )
        for output_path in [mul_output, pan_output]:
            if not output_path:
                continue
            calculate_raster_overviews(output_path, args.overview_scales)
    return state.current_files


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
    if args.run_from_existing and _existing_outputs_are_reusable(
        plan.output_paths,
        check_validity=args.run_from_existing_check_validity,
        validity_check_grid_size=args.validity_check_grid_size,
        log_to_console=args.log_to_console,
        step="fetch_atmosphere",
        scene_basename=state.scene.primary_basename,
    ):
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
    if args.run_from_existing and _existing_outputs_are_reusable(
        plan.output_paths,
        check_validity=args.run_from_existing_check_validity,
        validity_check_grid_size=args.validity_check_grid_size,
        log_to_console=args.log_to_console,
        step="atmospheric_correction",
        scene_basename=state.scene.primary_basename,
    ):
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
        state.current_step = "file_source"
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
        if args.run_from_existing and _existing_outputs_are_reusable(
            plan.output_paths,
            check_validity=args.run_from_existing_check_validity,
            validity_check_grid_size=args.validity_check_grid_size,
            log_to_console=args.log_to_console,
            step="orthorectification",
            scene_basename=state.scene.primary_basename,
        ):
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
            if args.run_from_existing and _existing_outputs_are_reusable(
                pan_plan.output_paths,
                check_validity=args.run_from_existing_check_validity,
                validity_check_grid_size=args.validity_check_grid_size,
                log_to_console=args.log_to_console,
                step="orthorectification_pan",
                scene_basename=state.scene.primary_basename,
            ):
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
    if args.run_from_existing and _existing_outputs_are_reusable(
        plan.output_paths,
        check_validity=args.run_from_existing_check_validity,
        validity_check_grid_size=args.validity_check_grid_size,
        log_to_console=args.log_to_console,
        step="pansharpen",
        scene_basename=state.scene.primary_basename,
    ):
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

    if args.run_from_existing and _existing_outputs_are_reusable(
        output_plan.output_paths + mask_plan.output_paths,
        check_validity=args.run_from_existing_check_validity,
        validity_check_grid_size=args.validity_check_grid_size,
        log_to_console=args.log_to_console,
        step="cloud_mask",
        scene_basename=state.scene.primary_basename,
    ):
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
    if args.run_from_existing and _existing_outputs_are_reusable(
        plan.output_paths,
        check_validity=args.run_from_existing_check_validity,
        validity_check_grid_size=args.validity_check_grid_size,
        log_to_console=args.log_to_console,
        step="alignment",
        scene_basename=state.scene.primary_basename,
    ):
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
    if required_outputs and not _existing_outputs_are_reusable(
        required_outputs,
        check_validity=args.skip_existing_check_validity,
        validity_check_grid_size=args.validity_check_grid_size,
        log_to_console=args.log_to_console,
        step="workflow",
        scene_basename=state.scene.primary_basename,
    ):
        return False
    if args.run_radiometric_normalization:
        return _existing_outputs_are_reusable(
            expected_outputs["final_raster"],
            check_validity=args.skip_existing_check_validity,
            validity_check_grid_size=args.validity_check_grid_size,
            log_to_console=args.log_to_console,
            step="workflow",
            scene_basename=state.scene.primary_basename,
        )
    if required_outputs:
        return True
    return _existing_outputs_are_reusable(
        expected_outputs["final_raster"],
        check_validity=args.skip_existing_check_validity,
        validity_check_grid_size=args.validity_check_grid_size,
        log_to_console=args.log_to_console,
        step="workflow",
        scene_basename=state.scene.primary_basename,
    )


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

    if args.run_file_source and not _is_temp_save_value(args.save_file_source):
        _extend("file_source")
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

    if args.run_file_source and _is_temp_save_value(args.save_file_source):
        _extend("file_source")
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
            "run_from_existing_check_validity": args.run_from_existing_check_validity,
            "skip_existing_check_validity": args.skip_existing_check_validity,
            "run_file_source": args.run_file_source,
            "run_fetch_atmosphere": args.run_fetch_atmosphere,
            "run_atmospheric_correction": args.run_atmospheric_correction,
            "run_orthorectification": args.run_orthorectification,
            "run_pansharpen": args.run_pansharpen,
            "run_cloud_mask": args.run_cloud_mask,
            "run_alignment": args.run_alignment,
            "run_seamline_metadata": args.run_seamline_metadata,
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
    _run_file_source_step(state, args)
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

    seamline_metadata_output = None
    if processed_states and args.run_seamline_metadata:
        seamline_metadata_output = _run_seamline_metadata_workflow(
            processed_states,
            args=args,
            reference_state=processed_states[0],
        )
        if seamline_metadata_output:
            log(
                f"Wrote seamline metadata {seamline_metadata_output}",
                enabled=args.log_to_console,
                step="seamline_metadata",
            )
        _apply_weighted_seamline_metadata_defaults(args, seamline_metadata_output)

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
    parser.add_argument("--run-from-existing-check-validity", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--validity-check-grid-size", type=int, default=0)
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
    parser.add_argument("--orthorectification-output-suffix", default="_ortho")
    parser.add_argument("--orthorectification-pan-output-suffix", default="_pan_ortho")
    parser.add_argument("--orthorectification-rpc-refinement-geojson")
    parser.add_argument("--pansharpen-output-suffix", default="_pansharpen")
    parser.add_argument("--skip-existing", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--skip-existing-check-validity", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--concurrent-processing", default=1)
    parser.add_argument("--overview-scales", nargs="+")
    parser.add_argument("--temp-dir")
    parser.add_argument("--keep-temp-dir", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--scratch-dir", dest="temp_dir", help=argparse.SUPPRESS)
    parser.add_argument("--run-file-source", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--save-file-source", default="$temp/file_source")
    parser.add_argument("--calculate-overviews-file-source", action=argparse.BooleanOptionalAction, default=False)
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
    parser.add_argument("--run-seamline-metadata", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--save-seamline-metadata", default="$temp/seamline_metadata.gpkg")
    parser.add_argument("--seamline-metadata-layer", default="footprints")
    parser.add_argument("--seamline-metadata-image-field-name", default="image")
    parser.add_argument(
        "--seamline-metadata-footprint-source",
        choices=["package_bounds", "calculate_bounds"],
        default="package_bounds",
    )
    parser.add_argument(
        "--seamline-metadata-calculate-bounds-eight-connected",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--run-radiometric-normalization", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--save-radiometric-normalization", default="$temp/radiometric_root.tif")
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
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument("--config-yaml")
    config_args, _ = config_parser.parse_known_args(raw_argv)

    config_defaults: Dict = {}
    if config_args.config_yaml:
        config_defaults = _normalize_config_defaults(_load_worldview_yaml_config(config_args.config_yaml))

    parser = _build_parser()
    if config_defaults:
        parser.set_defaults(**config_defaults)
    args, unknown_args = parser.parse_known_args(raw_argv)
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
    if args.validity_check_grid_size < 0:
        parser.error("--validity-check-grid-size must be >= 0.")
    per_scene_save_modes = {
        "temp_root",
        "temp_child",
        "output_root",
        "output_child",
        "input_relative",
        "absolute",
        "cwd_relative",
    }
    aggregate_single_output_modes = {"temp_child", "absolute", "cwd_relative"}
    for arg_name in (
        "save_file_source",
        "save_fetch_atmosphere",
        "save_atmospheric_correction",
        "save_orthorectification",
        "save_pansharpen",
        "save_cloud_mask",
        "save_alignment",
    ):
        _validate_save_target_value(
            getattr(args, arg_name),
            arg_name=arg_name,
            default="$temp",
            accepted_modes=per_scene_save_modes,
        )
    for arg_name, default_value in (
        ("save_seamline_metadata", "$temp/seamline_metadata.gpkg"),
        ("save_radiometric_normalization", "$temp/radiometric_root.tif"),
    ):
        _validate_save_target_value(
            getattr(args, arg_name),
            arg_name=arg_name,
            default=default_value,
            accepted_modes=aggregate_single_output_modes,
        )
    if (
        args.max_cloud_cover_to_process is not None
        and (args.max_cloud_cover_to_process < 0 or args.max_cloud_cover_to_process > 100)
    ):
        parser.error("--max-cloud-cover-to-process must be in [0, 100].")
    args.concurrent_processing = _resolve_concurrent_processing(args.concurrent_processing)
    if args.overview_scales is not None:
        if isinstance(args.overview_scales, str):
            args.overview_scales = [int(value.strip()) for value in args.overview_scales.split(",") if value.strip()]
        else:
            args.overview_scales = [int(value) for value in args.overview_scales]
    if (
        args.calculate_overviews_file_source
        or args.calculate_overviews_atmospheric_correction
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
    _parse_json_dict(args.cloud_mask_omnicloud_kwargs_json)
    radiometric_kwargs_json = _parse_json_dict(args.radiometric_normalization_kwargs_json)
    removed_spectralmatch_keys = {
        "matching_order": "steps",
        "match_steps": "steps",
        "seamline_steps": "steps",
        "seamline_method": "steps",
        "global_regression_output_images": "shared_temp_dir/global or shared_output_image_path when it is the final step",
        "local_block_adjustment_output_images": "shared_temp_dir/local or shared_output_image_path when it is the final step",
        "align_method": "steps",
        "align_rasters_output_images": "shared_temp_dir/aligned or shared_output_image_path when it is the final step",
        "voronoi_center_seamline_output_mask": "shared_output_image_path when it is the final step",
        "weighted_seamline_output_mask": "shared_output_image_path when it is the final step",
        "clip_method": "steps",
        "mask_rasters_output_images": "shared_temp_dir/clip or shared_output_image_path when it is the final step",
        "merge_method": "steps",
    }
    for old_key, replacement in removed_spectralmatch_keys.items():
        if getattr(args, f"match_{old_key}", None) is not None or old_key in radiometric_kwargs_json:
            parser.error(f"SpectralMatch no longer accepts {old_key}; use {replacement}.")
    save_radiometric_output_is_explicit = (
        "save_radiometric_normalization" in config_defaults
        or _explicit_cli_arg_present(raw_argv, "save_radiometric_normalization")
    )
    match_shared_output_is_explicit = (
        getattr(args, "match_shared_output_image_path", None) is not None
        or "shared_output_image_path" in radiometric_kwargs_json
    )
    if save_radiometric_output_is_explicit and match_shared_output_is_explicit:
        parser.error(
            "Cannot set both save_radiometric_normalization and "
            "match_shared_output_image_path/shared_output_image_path; use only one radiometric output path option."
        )
    group_by_basename_spec = _normalize_group_by_basename_spec(args.group_by_basename)
    if group_by_basename_spec is not None and save_radiometric_output_is_explicit:
        parser.error("Cannot set save_radiometric_normalization when group_by_basename is set; use group keys as output filenames.")
    if group_by_basename_spec is not None and match_shared_output_is_explicit:
        parser.error("Cannot set match_shared_output_image_path/shared_output_image_path when group_by_basename is set; use group keys as output filenames.")
    return _run_workflow(args)


if __name__ == "__main__":
    sys.exit(main())
