#!/usr/bin/env python3
"""CLI wrapper for the WorldView preprocessing workflow."""

from __future__ import annotations


import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from datetime import datetime
from functools import wraps
from typing import Any, Dict, List, Mapping, Optional

from osgeo import gdal
from tifffile import TiffFile
from wcmatch import fnmatch as wc_fnmatch
from wcmatch import glob

from vhrharmonize.providers.worldview.core import _worldview_image_source_files
from vhrharmonize.providers.standardized import materialize_scene_bounds
from vhrharmonize.cli.cli_helpers import _load_yaml_config
from vhrharmonize.io.geospatial import calculate_raster_overviews, get_image_percentile_value
from vhrharmonize.io.workflow_utils import (
    plan_step_outputs,
    remove_output_files,
    resolve_output_dir,
    resolve_temp_dir,
    resolve_relative_to_input,
)
from vhrharmonize.preprocess.atmospheric_correction import run_flaash, run_py6s
from vhrharmonize.preprocess.alignment import align_image_pair
from vhrharmonize.preprocess.cloudmasking import cloudmask_raster
from vhrharmonize.preprocess.concurrency import (
    _resolve_concurrent_processing, _resolve_concurrent_processing_backend, _make_dask_client,
)
from vhrharmonize.preprocess.fetch_external_data import (
    DEFAULT_OPENTOPOGRAPHY_DEMTYPE,
    DEFAULT_OPENTOPOGRAPHY_GLOBALDEM_ENDPOINT,
    download_opentopography_dem_for_bbox,
    fetch_modis_water_vapor_for_bbox,
    fetch_power_atmosphere_for_bbox,
)
from vhrharmonize.preprocess.helpers import (
    _log, _log_step_start, _log_image_start, _log_image_completed, _processing_step,
)
from vhrharmonize.preprocess.orthorectification import (
    gcp_refined_rpc_orthorectification,
    resolve_output_resolution_for_crs,
)
from vhrharmonize.preprocess.pansharpening import pansharpen_image
from vhrharmonize.preprocess.spectralmatch import DEFAULT_PIPELINE_STEPS, spectralmatch
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

SPECTRALMATCH_OVERVIEW_STEPS = {
    "joint_coregistration": "joint_coregistration_build_overviews",
    "global_regression": "global_regression_build_overviews",
    "local_block_adjustment": "local_block_adjustment_build_overviews",
    "merge": "merge_rasters_build_overviews",
}

WCMATCH_INPUT_FLAGS = (
    glob.GLOBSTAR
    | glob.BRACE
    | glob.EXTGLOB
    | glob.GLOBTILDE
    | glob.GLOBSTARLONG
    | glob.NEGATE
)

INPUT_FILE_STAGES = set(RASTER_STEP_ORDER)
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
    source_files: List[str] = field(default_factory=list)
    cleanup_step_outputs: Dict[str, List[str]] = field(default_factory=dict)


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
        return image.step_file_paths.get("file_source", image.tif_file)
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


def _normalize_input_file_stage_key(stage_key: str) -> str:
    """Normalize an input_file_glob stage key to an internal raster step."""
    if not isinstance(stage_key, str) or not stage_key.strip():
        raise ValueError("input_file_glob stage keys must be non-empty strings.")
    normalized = stage_key.strip().replace("-", "_")
    if normalized not in INPUT_FILE_STAGES:
        allowed = ", ".join(RASTER_STEP_ORDER)
        raise ValueError(f"Unsupported input_file_glob stage '{stage_key}'. Allowed stages: {allowed}")
    return normalized


def _input_file_stage_config_key(stage_name: str) -> str:
    """Return the external input_file_glob key for an internal stage name."""
    return _normalize_input_file_stage_key(stage_name)


def _normalize_input_file_glob_entries(input_file_globs: object) -> List[Dict[str, str]]:
    """Normalize input_file_glob to one-key dictionaries with string paths."""
    raw_entries = input_file_globs if isinstance(input_file_globs, list) else [input_file_globs]
    normalized_entries: List[Dict[str, str]] = []
    for entry in raw_entries:
        if isinstance(entry, str):
            normalized_entries.append({"file_source": entry})
            continue
        if not isinstance(entry, dict) or len(entry) != 1:
            raise ValueError("Each input_file_glob entry must be a one-key dictionary of stage: path.")
        stage_key, pattern = next(iter(entry.items()))
        if not isinstance(stage_key, str) or not isinstance(pattern, str):
            raise ValueError("Each input_file_glob entry must have string stage and path values.")
        normalized_entries.append({_input_file_stage_config_key(stage_key): pattern})
    return normalized_entries


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


def _build_spectralmatch_kwargs(args: argparse.Namespace) -> Dict:
    """Build SpectralMatch keyword arguments.
    Args:
        args: Parsed CLI arguments.
    Returns:
        SpectralMatch keyword arguments.
    """
    spectralmatch_kwargs = _parse_json_dict(args.spectralmatch_kwargs_json)
    match_kwargs = _collect_prefixed_kwargs(args, "match_")
    # Explicit null disables these upstream defaults instead of omitting the option.
    for key in ("shared_window_scales", "global_regression_pif_max_samples", "global_regression_pif_min_samples"):
        if hasattr(args, "match_" + key):
            match_kwargs[key] = getattr(args, "match_" + key)
    spectralmatch_kwargs.update(match_kwargs)
    spectralmatch_kwargs.setdefault("delete_temp_dir", getattr(args, "delete_temp_dir", True))
    overview_scales = getattr(args, "overview_scales", None)
    if overview_scales is not None:
        spectralmatch_kwargs.setdefault("shared_window_scales", tuple(overview_scales))
    for key, value in _spectralmatch_runtime_kwargs(args).items():
        spectralmatch_kwargs.setdefault(key, value)
    if getattr(args, "calculate_overviews_spectralmatch", False):
        conflicts = [
            "match_" + key for key, value in spectralmatch_kwargs.items()
            if key.endswith("_build_overviews") and value is True
        ]
        if conflicts:
            raise ValueError(
                "calculate_overviews_spectralmatch cannot be combined with enabled "
                + ", ".join(conflicts)
            )
        steps = spectralmatch_kwargs.get("steps", DEFAULT_PIPELINE_STEPS)
        overview_option = next(
            (SPECTRALMATCH_OVERVIEW_STEPS[step] for step in reversed(steps)
             if step in SPECTRALMATCH_OVERVIEW_STEPS),
            None,
        )
        if overview_option is None:
            raise ValueError("calculate_overviews_spectralmatch requires an overview-capable step in match_steps.")
        if not spectralmatch_kwargs.get("shared_window_scales"):
            raise ValueError("calculate_overviews_spectralmatch requires overview_scales or match_shared_window_scales.")
        spectralmatch_kwargs[overview_option] = True
    return spectralmatch_kwargs


def _spectralmatch_runtime_kwargs(args: argparse.Namespace) -> Dict:
    """Convert VHRHarmonize concurrency settings to SpectralMatch pipeline kwargs."""
    backend = getattr(args, "concurrent_processing_backend", "process_pool")
    if backend != "dask":
        return {"shared_concurrent_processing_backend": backend, "shared_dask_scheduler": None}

    scheduler_file = getattr(args, "dask_scheduler_file", None)
    scheduler = ("file", scheduler_file) if scheduler_file else (
        "address",
        getattr(args, "dask_scheduler_address", None),
    )
    return {
        "shared_concurrent_processing_backend": backend,
        "shared_dask_scheduler": scheduler,
        "shared_image_threads": None,
    }


def _spectralmatch_steps_include(args: argparse.Namespace, step_name: str) -> bool:
    """Return whether the SpectralMatch steps list includes a step."""
    steps = getattr(args, "match_steps", None)
    if steps is None:
        spectralmatch_kwargs = _parse_json_dict(getattr(args, "spectralmatch_kwargs_json", None))
        steps = spectralmatch_kwargs.get("steps")
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
    loaded = _load_yaml_config(config_yaml_path)

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
    if "keep_temp_dir" in normalized:
        raise ValueError("Unrecognized config key: keep_temp_dir")
    retired = [key for key in normalized if "radiometric_normalization" in key]
    if retired:
        raise ValueError("Use spectralmatch instead of radiometric_normalization in config keys: " + ", ".join(retired))
    if "input_file_glob" in normalized:
        normalized["input_file_glob"] = _normalize_input_file_glob_entries(normalized["input_file_glob"])
    for list_key in ("filter_basename", "match_steps"):
        if list_key in normalized and isinstance(normalized[list_key], str):
            normalized[list_key] = [normalized[list_key]]
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
        _log(
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
        _log(
            f"Reusing DEM {os.path.basename(dem_output_path)}",
            enabled=args.log_to_console,
            step="dem",
            scene_basename=state.scene.primary_basename,
        )
    else:
        scene_bbox = materialize_scene_bounds(mul_image.standardized_metadata.source_metadata).bounds
        download_opentopography_dem_for_bbox(
            min_lon=scene_bbox[0],
            min_lat=scene_bbox[1],
            max_lon=scene_bbox[2],
            max_lat=scene_bbox[3],
            output_tif_path=dem_output_path,
            api_key=args.dem_online_api_key,
            demtype=args.dem_online_source,
            endpoint=args.dem_online_api_endpoint,
            timeout_s=args.dem_online_timeout_s,
            log_to_console=args.log_to_console,
            scene_basename=state.scene.primary_basename,
        )
        _log(
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


def _collect_input_files_by_stage(input_file_globs: object) -> Dict[str, List[str]]:
    """Collect input files from stage-keyed glob patterns."""
    matched_by_stage: Dict[str, List[str]] = {}
    for entry in _normalize_input_file_glob_entries(input_file_globs):
        stage_key, pattern = next(iter(entry.items()))
        stage_name = _normalize_input_file_stage_key(stage_key)
        matches = glob.glob(pattern, flags=WCMATCH_INPUT_FLAGS)
        matched_by_stage.setdefault(stage_name, []).extend(path for path in matches if os.path.isfile(path))
    return {stage: sorted({os.path.abspath(path) for path in paths}) for stage, paths in matched_by_stage.items()}


def _collect_input_files(input_file_globs: object) -> List[str]:
    """Collect unique input files from all stage-keyed glob patterns."""
    matched_by_stage = _collect_input_files_by_stage(input_file_globs)
    return sorted({path for paths in matched_by_stage.values() for path in paths})


def _load_worldview_scenes_from_stage_paths(
    input_files_by_stage: Mapping[str, List[str]],
    *,
    filter_basenames: Optional[List[str]],
) -> List[WorldViewScene]:
    """Load scenes while recording which workflow stage each input satisfies."""
    scenes_by_key: Dict[tuple[str, str], WorldViewScene] = {}
    for stage_name in RASTER_STEP_ORDER:
        stage_files = input_files_by_stage.get(stage_name) or []
        if not stage_files:
            continue
        for stage_scene in load_worldview_scenes_from_tif_files(stage_files, filter_basenames=filter_basenames):
            key = (stage_scene.scene_id, stage_scene.catalog_id)
            existing = scenes_by_key.get(key)
            if existing is None:
                scenes_by_key[key] = stage_scene
                existing = stage_scene
            for image in stage_scene.iter_images():
                target_image = existing.get_image(image.image_role or "")
                if target_image is None:
                    existing.set_image(image)
                    target_image = image
                target_image.step_file_paths[stage_name] = image.tif_file
                existing.step_outputs.setdefault(stage_name, [])
                if image.tif_file not in existing.step_outputs[stage_name]:
                    existing.step_outputs[stage_name].append(image.tif_file)
    return [scenes_by_key[key] for key in sorted(scenes_by_key)]


def _process_scenes_with_process_pool(
    scenes: List[WorldViewScene],
    args: argparse.Namespace,
    worker_count: int,
) -> List[SceneWorkflowState]:
    """Run independent scene preprocessing with ProcessPoolExecutor."""
    ordered_results: List[SceneWorkflowState | None] = [None] * len(scenes)
    with ProcessPoolExecutor(max_workers=min(worker_count, len(scenes))) as executor:
        future_to_index = {
            executor.submit(_process_scene, scene, args): index
            for index, scene in enumerate(scenes)
        }
        for future in as_completed(future_to_index):
            index = future_to_index[future]
            ordered_results[index] = future.result()
    return [state for state in ordered_results if state is not None]


def _process_scenes_with_dask(
    scenes: List[WorldViewScene],
    args: argparse.Namespace,
) -> List[SceneWorkflowState]:
    """Run independent scene preprocessing on an existing Dask cluster."""
    client = _make_dask_client(args)
    try:
        futures = client.map(_process_scene, scenes, [args] * len(scenes))
        return list(client.gather(futures))
    finally:
        client.close()


def _process_scenes(
    scenes: List[WorldViewScene],
    args: argparse.Namespace,
) -> List[SceneWorkflowState]:
    """Process scenes independently before aggregate workflow steps."""
    args.scene_indices = {
        getattr(scene, "primary_basename", str(scene)): index
        for index, scene in enumerate(scenes, start=1)
    }
    args.scene_total = len(scenes)
    if getattr(args, "delete_temp_steps_proactively", False) or getattr(args, "delete_temp_dir", False):
        args._cleanup_source_files = _dedupe_paths([
            path for scene in scenes for image in scene.iter_images()
            for path in _worldview_image_source_files(image)
        ] + [path for scene in scenes for paths in scene.step_outputs.values() for path in paths])
    if getattr(args, "delete_temp_steps_proactively", False):
        # Sweep all discovered scenes before dispatch, including scenes that will
        # be skipped or may never reach a worker if another scene fails.
        args._cleanup_source_keys = _cleanup_file_keys(_files_with_sidecars(args._cleanup_source_files))
        if args.temp_dir:
            for scene in scenes:
                state = _initialize_scene_state(scene, args, validate_inputs=False)
                _cleanup_completed_scene_temp_steps(state, args)
    worker_count = _resolve_concurrent_processing(args.concurrent_processing)
    backend = _resolve_concurrent_processing_backend(args.concurrent_processing_backend)
    if backend == "dask":
        if worker_count != 1:
            raise ValueError("concurrent_processing must be 1 when concurrent_processing_backend is 'dask'.")
        _log(
            f"Running per-scene processing with Dask tasks={len(scenes)}",
            enabled=args.log_to_console,
            step="workflow",
        )
        return _process_scenes_with_dask(scenes, args)
    if worker_count <= 1 or len(scenes) <= 1:
        return [_process_scene(scene, args) for scene in scenes]
    _log(
        f"Running per-scene processing with {min(worker_count, len(scenes))} processes",
        enabled=args.log_to_console,
        step="workflow",
    )
    return _process_scenes_with_process_pool(scenes, args, worker_count)


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
    _log(" | ".join(parts), enabled=enabled, step=step, scene_basename=scene_basename)


def _logged_scene_step(step):
    """Give each scene step a consistent lifecycle, including reused outputs."""
    def _decorate(function):
        @wraps(function)
        def _wrapped(state, args):
            if not args.log_to_console or not getattr(args, f"run_{step}", False):
                return function(state, args)
            if step == "file_source" and state.current_step != "file_source":
                return function(state, args)
            expected_outputs = _get_expected_scene_step_outputs(state, args)
            inputs = list(state.current_files)
            if step == "pansharpen" and state.pan_ortho_path:
                inputs.append(state.pan_ortho_path)
            elif step == "orthorectification" and args.run_pansharpen and state.scene.pan_image:
                inputs.append(state.scene.pan_image.tif_file)
            elif step == "alignment":
                inputs.append(args.alignment_fixed_image)
            scene_basename = state.scene.primary_basename
            with _processing_step(
                step, scene_basename, inputs,
                _scene_step_expected_outputs(expected_outputs, step),
                enabled=args.log_to_console,
                index=getattr(args, "scene_indices", {}).get(scene_basename, 1),
                total=getattr(args, "scene_total", 1),
                announce_step=False,
            ):
                return function(state, args)
        return _wrapped
    return _decorate


def _classify_save_target(save_value: Optional[str], *, default: str) -> tuple[str, str]:
    """Return a normalized save target and target kind."""
    normalized = os.path.expanduser(
        (save_value if save_value not in (None, "") else default).strip()
    )
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


def _is_gdal_raster_path(path: str) -> bool:
    """Return whether a path should be validated as a GDAL raster."""
    extension = os.path.splitext(path)[1].lower()
    return extension in {".tif", ".tiff", ".dat", ".img", ".vrt"}


def _tiff_truncation_reason(path: str) -> str | None:
    """Check strip/tile bounds in all TIFF image directories without decoding pixels."""
    try:
        with TiffFile(path) as tif:
            file_size = tif.filehandle.size
            pending = [tif.pages]
            while pending:
                for page in pending.pop():
                    for offset, byte_count in zip(page.dataoffsets, page.databytecounts):
                        if offset and byte_count and offset + byte_count > file_size:
                            return (
                                f"Truncated TIFF: offset {offset} + byte count {byte_count} "
                                f"exceeds file size {file_size} (IFD at {page.offset})"
                            )
                    if page.pages is not None:
                        pending.append(page.pages)
    except (OSError, ValueError, NotImplementedError, struct.error):
        # Inconclusive TIFF inspection: let the existing GDAL checks decide.
        pass
    return None


def _gdal_raster_is_valid(path: str, *, validity_check_grid_size: int = 0) -> tuple[bool, str | None]:
    """Check raster readability and TIFF strip/tile bounds.
    Args:
        path: Raster path to validate.
        validity_check_grid_size: Pixel sampling grid size. 0 disables pixel validity sampling.
    Returns:
        Tuple of validity and optional reason.
    """
    if not os.path.exists(path):
        return False, "missing"
    dataset = None
    try:
        dataset = gdal.OpenEx(path, gdal.OF_RASTER)
        if dataset is None:
            return False, "GDAL open failed"
        if dataset.GetDriver().ShortName == "GTiff":
            reason = _tiff_truncation_reason(path)
            if reason:
                return False, reason
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


def _json_file_is_valid(path: str) -> tuple[bool, str | None]:
    """Check that a file can be parsed as JSON without validating its contents."""
    try:
        with open(path, "r", encoding="utf-8") as handle:
            json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        return False, str(exc)
    return True, None


def _existing_output_failures(
    output_paths: List[str],
    *,
    check_validity: bool,
    validity_check_grid_size: int,
    log_to_console: bool,
    step: str,
    scene_basename: str | None = None,
) -> Dict[str, str]:
    """Inspect every expected output without changing files, including partial sets."""
    failures: Dict[str, str] = {}
    for output_path in output_paths:
        if not os.path.exists(output_path):
            failures[output_path] = "missing"
            continue
        if not check_validity:
            continue
        if os.path.splitext(output_path)[1].lower() == ".json":
            is_valid, reason = _json_file_is_valid(output_path)
        elif _is_gdal_raster_path(output_path):
            is_valid, reason = _gdal_raster_is_valid(
                output_path,
                validity_check_grid_size=validity_check_grid_size,
            )
        else:
            continue
        if not is_valid:
            failures[output_path] = reason or "invalid"
            _log_step_plan(
                step,
                outputs=[output_path],
                message=f"Existing output invalid ({reason})",
                enabled=log_to_console,
                scene_basename=scene_basename,
            )
    return failures


def _existing_outputs_are_reusable(
    output_paths: List[str],
    *,
    check_validity: bool,
    validity_check_grid_size: int,
    log_to_console: bool,
    step: str,
    scene_basename: str | None = None,
) -> bool:
    """Return whether all outputs exist and pass enabled checks, without deleting files."""
    return bool(output_paths) and not _existing_output_failures(
        output_paths,
        check_validity=check_validity,
        validity_check_grid_size=validity_check_grid_size,
        log_to_console=log_to_console,
        step=step,
        scene_basename=scene_basename,
    )


def _prepare_step_outputs(
    output_paths: List[str],
    *,
    input_paths: List[str],
    args: argparse.Namespace,
    step: str,
    scene_basename: str | None = None,
    remove_invalid: bool = True,
) -> bool:
    """Check reuse and remove invalid outputs before executing a processing step.

    Valid outputs are preserved. With validity checking disabled, files are only
    checked for existence. Inspection/counting code must use the read-only helper.
    Set remove_invalid=False when the caller will bypass generation of these outputs.
    """
    failures = _existing_output_failures(
        output_paths,
        check_validity=args.run_from_existing_check_validity,
        validity_check_grid_size=args.validity_check_grid_size,
        log_to_console=args.log_to_console,
        step=step,
        scene_basename=scene_basename,
    )
    if remove_invalid and args.run_from_existing_check_validity and failures:
        remove_output_files(failures, input_paths=input_paths)
    return args.run_from_existing and bool(output_paths) and not failures


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
    if args.run_spectralmatch:
        step_dirs["spectralmatch"] = _resolve_single_output_save_path(
            args.save_spectralmatch,
            default="$temp/spectralmatch_root.tif",
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
    pan_image = state.scene.get_image("pan")
    if args.run_file_source:
        file_source_map: Dict[str, str] = {}
        file_source_map.update(_image_source_file_map(mul_image, state.step_dirs["file_source"]))
        if pan_image is not None:
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


def _step_will_run_from(start_step: str, step_name: str, args: argparse.Namespace) -> bool:
    """Return whether a raster step remains after the supplied input stage."""
    return (
        RASTER_STEP_ORDER.index(start_step) < RASTER_STEP_ORDER.index(step_name)
        and bool(getattr(args, f"run_{step_name}", False))
    )


def _remaining_scene_input_requirements(start_step: str, args: argparse.Namespace) -> set[str]:
    """Return source-side files required by remaining enabled raster steps."""
    required: set[str] = set()
    if _step_will_run_from(start_step, "atmospheric_correction", args):
        required.add("mul_metadata")
    if _step_will_run_from(start_step, "orthorectification", args):
        required.add("mul_metadata")
        if args.run_pansharpen:
            required.update({"pan_image", "pan_metadata"})
    elif _step_will_run_from(start_step, "pansharpen", args):
        required.add("pan_ortho")
    return required


def _validate_remaining_scene_inputs(
    scene: WorldViewScene,
    *,
    start_step: str,
    args: argparse.Namespace,
) -> None:
    """Validate files required by remaining enabled steps."""
    required = _remaining_scene_input_requirements(start_step, args)
    mul_image = _require_scene_image(scene, "mul")
    pan_image = scene.get_image("pan")
    if "mul_metadata" in required and mul_image.standardized_metadata is None:
        raise ValueError(f"WorldView scene is missing multispectral metadata: {scene.scene_id}_{scene.catalog_id}")
    if "pan_image" in required and pan_image is None:
        raise ValueError(
            "WorldView scene needs a panchromatic image before pansharpening: "
            f"{scene.scene_id}_{scene.catalog_id}"
        )
    if "pan_metadata" in required and pan_image is not None and pan_image.standardized_metadata is None:
        raise ValueError(f"WorldView scene is missing panchromatic metadata: {scene.scene_id}_{scene.catalog_id}")
    if "pan_ortho" in required:
        raise ValueError(
            "Input stage is before pansharpen but no prior orthorectification step will create a panchromatic ortho image."
        )


def _initialize_scene_state(
    scene: WorldViewScene, args: argparse.Namespace, *, validate_inputs: bool = True,
) -> SceneWorkflowState:
    """Initialize workflow state for a scene.
    Args:
        scene: Scene to initialize.
        args: Parsed CLI arguments.
        validate_inputs: Check remaining processing inputs; false for cleanup inspection only.
    Returns:
        Initialized scene workflow state.
    """
    start_step = _get_scene_input_start_step(scene)
    if validate_inputs:
        _validate_remaining_scene_inputs(scene, start_step=start_step, args=args)
    start_path = _get_worldview_scene_step_path(scene, "mul", start_step)

    state = SceneWorkflowState(
        scene=scene,
        step_dirs=_resolve_scene_step_dirs(args, scene),
        current_files=[start_path],
        current_step=start_step,
    )
    if getattr(args, "delete_temp_steps_proactively", False) or getattr(args, "delete_temp_dir", False):
        state.source_files = _dedupe_paths([
            path for image in scene.iter_images() for path in _worldview_image_source_files(image)
        ] + [path for paths in scene.step_outputs.values() for path in paths])
        if args.dem_file_path not in (None, "", "online"):
            state.source_files.append(resolve_relative_to_input(args.dem_file_path, os.path.dirname(scene.mul_image.tif_file)))
    if getattr(args, "delete_temp_steps_proactively", False):
        state.cleanup_step_outputs = _discovered_step_outputs(state, args)
    return state


def _get_scene_input_start_step(scene: WorldViewScene) -> str:
    """Return the latest raster step supplied by input_file_glob for a scene."""
    for step_name in reversed(RASTER_STEP_ORDER):
        if scene.step_outputs.get(step_name):
            return step_name
    return "file_source"


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
    """Normalize a spectralmatch grouping specification.
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
    return _validate_spectralmatch_group_spec(raw_spec)


def _validate_spectralmatch_group_spec(group_spec: object) -> Dict[str, object]:
    """Validate the named spectralmatch grouping JSON shape."""
    if not isinstance(group_spec, dict) or not group_spec:
        raise ValueError("group_by_basename must be a non-empty object with output filename keys.")
    normalized: Dict[str, object] = {}
    for output_name, value in group_spec.items():
        if not isinstance(output_name, str) or not output_name:
            raise ValueError("Each spectralmatch group key must be a non-empty output filename string.")
        normalized[output_name] = _validate_spectralmatch_group_value(value)
    return normalized


def _validate_spectralmatch_group_value(value: object) -> object:
    """Validate a spectralmatch group value."""
    if isinstance(value, str):
        if not (value.startswith("auto:") or value.startswith("file:")):
            raise ValueError("SpectralMatch group strings must start with 'auto:' or 'file:'.")
        if value in {"auto:", "file:"}:
            raise ValueError("SpectralMatch group strings must include a pattern or path after the prefix.")
        return value
    if isinstance(value, list):
        if not value:
            raise ValueError("SpectralMatch group lists must not be empty.")
        return [_validate_spectralmatch_group_item(item) for item in value]
    raise ValueError("SpectralMatch group values must be strings or lists.")


def _validate_spectralmatch_group_item(item: object) -> object:
    """Validate an item inside a spectralmatch group list."""
    if isinstance(item, str):
        return _validate_spectralmatch_group_value(item)
    if isinstance(item, dict):
        return _validate_spectralmatch_group_spec(item)
    raise ValueError("SpectralMatch group list items must be prefixed strings or nested group objects.")


def _match_spectralmatch_input_patterns(pattern: str, available_paths: List[str]) -> List[str]:
    """Match spectralmatch input patterns against available paths.
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
        raise ValueError(f"No SpectralMatch inputs matched pattern: {pattern}")
    return matches


def _resolve_spectralmatch_input_token(token: str, available_paths: List[str]) -> List[str]:
    """Resolve a prefixed spectralmatch group token into input paths."""
    source, value = token.split(":", 1)
    if source == "auto":
        return _match_spectralmatch_input_patterns(value, available_paths)
    if source == "file":
        return [value]
    raise ValueError("SpectralMatch group strings must start with 'auto:' or 'file:'.")


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


def _resolve_spectralmatch_group_output_path(
    *,
    output_name: str,
    temp_root: str,
    output_root: str,
) -> str:
    """Resolve a spectralmatch group output path.
    Args:
        output_name: Group output filename.
        temp_root: Temp root directory.
        output_root: Output root directory.
    Returns:
        Resolved spectralmatch group output path.
    """
    del temp_root
    group_output_path = os.path.join(output_root, os.path.basename(output_name))
    os.makedirs(os.path.dirname(group_output_path) or ".", exist_ok=True)
    return group_output_path


def _run_named_spectralmatch_group(
    output_name: str,
    group_value: object,
    *,
    available_paths: List[str],
    args: argparse.Namespace,
    temp_root: str,
    output_root: str,
) -> str:
    """Run one named spectralmatch group and return its output path."""
    child_inputs: List[str] = []
    if isinstance(group_value, str):
        child_inputs.extend(_resolve_spectralmatch_input_token(group_value, available_paths))
    elif isinstance(group_value, list):
        for item in group_value:
            if isinstance(item, str):
                child_inputs.extend(_resolve_spectralmatch_input_token(item, available_paths))
            elif isinstance(item, dict):
                for child_output_name, child_value in item.items():
                    child_inputs.append(
                        _run_named_spectralmatch_group(
                            child_output_name,
                            child_value,
                            available_paths=available_paths,
                            args=args,
                            temp_root=temp_root,
                            output_root=output_root,
                        )
                    )
            else:
                raise ValueError("SpectralMatch group list items must be prefixed strings or nested group objects.")
    else:
        raise ValueError("SpectralMatch group values must be strings or lists.")

    # Only full-resolution tiles feed parent groups; pyramid subfolders are excluded.
    expanded_inputs = []
    for path in child_inputs:
        if os.path.isdir(path):
            expanded_inputs.extend(sorted(
                str(tile) for tile in Path(path).iterdir()
                if tile.is_file() and tile.suffix.lower() in {".tif", ".tiff"}
            ))
        else:
            expanded_inputs.append(path)
    child_inputs = _dedupe_paths(expanded_inputs)
    spectralmatch_kwargs = _build_spectralmatch_kwargs(args)
    spectralmatch_kwargs.pop("shared_output_image_path", None)
    group_output_path = _resolve_spectralmatch_group_output_path(
        output_name=output_name,
        temp_root=temp_root,
        output_root=output_root,
    )
    spectralmatch_kwargs.setdefault("shared_input_images", child_inputs)
    spectralmatch_kwargs.setdefault("shared_output_image_path", group_output_path)
    spectralmatch_kwargs.setdefault("shared_temp_dir", os.path.join(temp_root, "spectralmatch"))
    spectralmatch_kwargs.setdefault("shared_debug_logs", args.log_to_console)
    spectralmatch_kwargs.setdefault("shared_output_dtype", args.dtype)
    _log_step_plan(
        "spectralmatch",
        inputs=child_inputs,
        outputs=[group_output_path],
        message=f"Running SpectralMatch group {output_name}",
        enabled=args.log_to_console,
    )
    spectralmatch(
        method=args.spectralmatch_method,
        log_to_console=args.log_to_console,
        **spectralmatch_kwargs,
    )
    return group_output_path


def _run_named_spectralmatch_groups(
    group_spec: Dict[str, object],
    *,
    available_paths: List[str],
    args: argparse.Namespace,
    temp_root: str,
    output_root: str,
) -> Optional[str]:
    """Run named spectralmatch groups in order and return the final output path."""
    final_output: Optional[str] = None
    for output_name, group_value in group_spec.items():
        final_output = _run_named_spectralmatch_group(
            output_name,
            group_value,
            available_paths=available_paths,
            args=args,
            temp_root=temp_root,
            output_root=output_root,
        )
    return final_output


def _run_default_spectralmatch(
    available_paths: List[str],
    *,
    args: argparse.Namespace,
    output_path: str,
    temp_root: str,
) -> str:
    """Run one default SpectralMatch over all scene outputs."""
    spectralmatch_kwargs = _build_spectralmatch_kwargs(args)
    group_output_path = str(spectralmatch_kwargs.get("shared_output_image_path") or output_path)
    spectralmatch_kwargs.setdefault("shared_input_images", available_paths)
    spectralmatch_kwargs.setdefault("shared_output_image_path", group_output_path)
    spectralmatch_kwargs.setdefault("shared_temp_dir", os.path.join(temp_root, "spectralmatch"))
    spectralmatch_kwargs.setdefault("shared_debug_logs", args.log_to_console)
    spectralmatch_kwargs.setdefault("shared_output_dtype", args.dtype)
    _log_step_plan(
        "spectralmatch",
        inputs=available_paths,
        outputs=[group_output_path],
        message="Running SpectralMatch",
        enabled=args.log_to_console,
    )
    spectralmatch(
        method=args.spectralmatch_method,
        log_to_console=args.log_to_console,
        **spectralmatch_kwargs,
    )
    return group_output_path


def _run_spectralmatch_workflow(
    scene_output_paths: List[str],
    *,
    args: argparse.Namespace,
    reference_state: SceneWorkflowState,
) -> Optional[str]:
    """Run grouped SpectralMatch.
    Args:
        scene_output_paths: Per-scene output raster paths.
        args: Parsed CLI arguments.
        reference_state: Reference scene state used for directory resolution.
    Returns:
        Final spectralmatch output path or None.
    """
    if not args.run_spectralmatch:
        return None
    available_paths = _dedupe_paths([str(path) for path in scene_output_paths if str(path)])
    if not available_paths:
        raise ValueError("No scene outputs were available for SpectralMatch.")

    with _processing_step(
        "spectralmatch", "all_scenes", available_paths, [],
        enabled=args.log_to_console, index=len(available_paths),
        total=getattr(args, "scene_total", len(available_paths)),
    ):
        group_spec = _normalize_group_by_basename_spec(args.group_by_basename)
        if group_spec is not None:
            return _run_named_spectralmatch_groups(
                group_spec,
                available_paths=available_paths,
                args=args,
                temp_root=reference_state.step_dirs["temp_root"],
                output_root=reference_state.step_dirs["output_root"],
            )

        return _run_default_spectralmatch(
            available_paths,
            args=args,
            output_path=reference_state.step_dirs["spectralmatch"],
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
    if args.run_from_existing and not args.run_from_existing_check_validity and _existing_outputs_are_reusable(
        [output_path],
        check_validity=False,
        validity_check_grid_size=args.validity_check_grid_size,
        log_to_console=args.log_to_console,
        step="seamline_metadata",
    ):
        _log_step_start("seamline_metadata", enabled=args.log_to_console)
        total = len({os.path.basename(state.current_files[0]) for state in states if state.current_files})
        _log(f"Already processed {total}/{getattr(args, 'scene_total', total)}; reusing metadata without validity checking", enabled=args.log_to_console, step="seamline_metadata")
        return output_path

    return write_seamline_metadata_gpkg(
        states,
        output_path,
        layer=args.seamline_metadata_layer,
        image_field_name=args.seamline_metadata_image_field_name,
        scene_total=getattr(args, "scene_total", len(states)),
        footprint_source=args.seamline_metadata_footprint_source,
        calculate_bounds_eight_connected=args.seamline_metadata_calculate_bounds_eight_connected,
        epsg=args.epsg,
        run_from_existing_check_validity=(
            args.run_from_existing and args.run_from_existing_check_validity
        ),
        log_to_console=args.log_to_console,
        concurrent_processing=args.concurrent_processing,
        concurrent_processing_backend=args.concurrent_processing_backend,
        dask_scheduler_file=getattr(args, "dask_scheduler_file", None),
        dask_scheduler_address=getattr(args, "dask_scheduler_address", None),
    )


def _apply_weighted_seamline_metadata_defaults(args: argparse.Namespace, seamline_metadata_output: str | None) -> None:
    """Default weighted seamline inputs from the generated WorldView metadata GPKG."""
    if not seamline_metadata_output:
        return
    if not _spectralmatch_steps_include(args, "weighted_seamline"):
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
        command_template: Command argument template to execute.
        input_image_path: Input raster path.
        output_image_path: Output raster path.
        scene_root_path: Scene root directory.
        image_basename: Scene image basename.
        log_to_console: Whether to emit console logs.
        scene_basename: Optional scene basename for log prefixes.
    Returns:
        None.
    """
    command = [
        argument.format(
            input=input_image_path,
            output=output_image_path,
            scene_root=scene_root_path,
            image_basename=image_basename,
        )
        for argument in shlex.split(command_template)
    ]
    _log(
        "Running external cloud mask command",
        enabled=log_to_console,
        step="cloud_mask",
        scene_basename=scene_basename,
    )
    subprocess.run(command, check=True)  # nosec B603


def _image_source_file_map(image: Optional[WorldViewImage], output_dir: str) -> Dict[str, str]:
    """Return source bundle files mapped to their staged output paths."""
    if image is None:
        return {}
    source_paths = _worldview_image_source_files(image)
    path_map: Dict[str, str] = {}
    for source_path in _dedupe_paths([path for path in source_paths if os.path.isfile(path)]):
        path_map[os.path.abspath(source_path)] = os.path.join(output_dir, os.path.basename(source_path))
    return path_map


def _set_image_file_paths_from_source_map(image: Optional[WorldViewImage], path_map: Mapping[str, str]) -> None:
    """Update a WorldView image to point at staged source bundle files."""
    if image is None:
        return
    for attr_name in ("tif_file", "imd_file", "til_file"):
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


@_logged_scene_step("file_source")
def _run_file_source_step(state: SceneWorkflowState, args: argparse.Namespace) -> List[str]:
    """Run the file_source staging step."""
    if not args.run_file_source or state.current_step != "file_source":
        return state.current_files
    output_dir = state.step_dirs["file_source"]
    path_map: Dict[str, str] = {}
    path_map.update(_image_source_file_map(state.scene.mul_image, output_dir))
    path_map.update(_image_source_file_map(state.scene.pan_image, output_dir))
    output_paths = list(path_map.values())
    expected_outputs = _get_expected_scene_step_outputs(state, args)

    if _prepare_step_outputs(
        output_paths,
        input_paths=list(path_map),
        args=args,
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
        _log(
            "Calculating overviews for step file_source",
            enabled=args.log_to_console,
            step="overviews",
            scene_basename=state.scene.primary_basename,
        )
        for output_path in [mul_output, pan_output]:
            if not output_path:
                continue
            calculate_raster_overviews(
                output_path, args.overview_scales,
                log_to_console=args.log_to_console,
                scene_basename=state.scene.primary_basename,
                scene_index=getattr(args, "scene_indices", {}).get(state.scene.primary_basename, 1),
                scene_total=getattr(args, "scene_total", 1),
            )
    return state.current_files


@_logged_scene_step("fetch_atmosphere")
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
    if _prepare_step_outputs(
        plan.output_paths,
        input_paths=[mul_image.tif_file],
        args=args,
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

    scene_bbox = materialize_scene_bounds(mul_image.standardized_metadata.source_metadata).bounds
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
            min_lon=scene_bbox[0],
            min_lat=scene_bbox[1],
            max_lon=scene_bbox[2],
            max_lat=scene_bbox[3],
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
            min_lon=scene_bbox[0],
            min_lat=scene_bbox[1],
            max_lon=scene_bbox[2],
            max_lat=scene_bbox[3],
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


@_logged_scene_step("atmospheric_correction")
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
    if mul_image is None or mul_image.standardized_metadata is None:
        raise ValueError("WorldView scene is missing multispectral inputs for atmospheric correction.")
    plan = plan_step_outputs(
        state.current_files,
        output_dir=state.step_dirs["atmospheric_correction"],
        suffix=args.atmospheric_correction_output_suffix,
        extension=_get_atmospheric_extension(args),
        skip_existing=False,
    )
    if _prepare_step_outputs(
        plan.output_paths,
        input_paths=plan.input_paths,
        args=args,
        step="atmospheric_correction",
        scene_basename=state.scene.primary_basename,
        remove_invalid=(
            args.atmospheric_method != "none"
            and not (args.atmospheric_method == "flaash" and args.skip_flaash)
        ),
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

    if args.atmospheric_method == "flaash" and args.skip_flaash:
        if not args.existing_flaash_input:
            raise ValueError("--existing-flaash-input is required when --skip-flaash is set.")
        state.current_files = _register_step_outputs(state, "atmospheric_correction", [args.existing_flaash_input])
        state.current_step = "atmospheric_correction"
        return state.current_files

    resolved_dem_file_path = _resolve_scene_dem_file_path(state, args)
    if resolved_dem_file_path is None:
        raise ValueError("DEM file path is required for atmospheric correction.")

    if args.atmospheric_method == "flaash":
        footprint_geometry = materialize_scene_bounds(mul_image.standardized_metadata.source_metadata)
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
            footprint_geometry=footprint_geometry,
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
        footprint_geometry = materialize_scene_bounds(mul_image.standardized_metadata.source_metadata)
        ground_elevation_m = get_image_percentile_value(
            resolved_dem_file_path,
            percentile=args.flaash_dem_ground_percentile,
            mask=footprint_geometry,
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
        _log(
            "Calculating overviews for step atmospheric_correction",
            enabled=args.log_to_console,
            step="overviews",
            scene_basename=state.scene.primary_basename,
        )
        for output_path in plan.output_paths:
            calculate_raster_overviews(
                output_path, args.overview_scales,
                log_to_console=args.log_to_console,
                scene_basename=state.scene.primary_basename,
                scene_index=getattr(args, "scene_indices", {}).get(state.scene.primary_basename, 1),
                scene_total=getattr(args, "scene_total", 1),
            )
    state.current_step = "atmospheric_correction"
    return state.current_files


@_logged_scene_step("orthorectification")
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
        if _prepare_step_outputs(
            plan.output_paths,
            input_paths=plan.input_paths + [resolved_dem_file_path],
            args=args,
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
            _log(
                "Calculating overviews for step orthorectification",
                enabled=args.log_to_console,
                step="overviews",
                scene_basename=state.scene.primary_basename,
            )
            for output_path in plan.output_paths:
                calculate_raster_overviews(
                    output_path, args.overview_scales,
                    log_to_console=args.log_to_console,
                    scene_basename=state.scene.primary_basename,
                    scene_index=getattr(args, "scene_indices", {}).get(state.scene.primary_basename, 1),
                    scene_total=getattr(args, "scene_total", 1),
                )

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
            if _prepare_step_outputs(
                pan_plan.output_paths,
                input_paths=pan_plan.input_paths + [resolved_dem_file_path],
                args=args,
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
                _log(
                    "Calculating overviews for step orthorectification_pan",
                    enabled=args.log_to_console,
                    step="overviews",
                    scene_basename=state.scene.primary_basename,
                )
                for output_path in pan_plan.output_paths:
                    calculate_raster_overviews(
                        output_path, args.overview_scales,
                        log_to_console=args.log_to_console,
                        scene_basename=state.scene.primary_basename,
                        scene_index=getattr(args, "scene_indices", {}).get(state.scene.primary_basename, 1),
                        scene_total=getattr(args, "scene_total", 1),
                    )

    state.current_step = "orthorectification"
    return state.current_files


@_logged_scene_step("pansharpen")
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
    if _prepare_step_outputs(
        plan.output_paths,
        input_paths=plan.input_paths + [state.pan_ortho_path],
        args=args,
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
        _log(
            "Calculating overviews for step pansharpen",
            enabled=args.log_to_console,
            step="overviews",
            scene_basename=state.scene.primary_basename,
        )
        for output_path in plan.output_paths:
            calculate_raster_overviews(
                output_path, args.overview_scales,
                log_to_console=args.log_to_console,
                scene_basename=state.scene.primary_basename,
                scene_index=getattr(args, "scene_indices", {}).get(state.scene.primary_basename, 1),
                scene_total=getattr(args, "scene_total", 1),
            )
    state.current_step = "pansharpen"
    return state.current_files


@_logged_scene_step("cloud_mask")
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

    if _prepare_step_outputs(
        output_plan.output_paths + mask_plan.output_paths,
        input_paths=output_plan.input_paths,
        args=args,
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
        _log(
            "Calculating overviews for step cloud_mask",
            enabled=args.log_to_console,
            step="overviews",
            scene_basename=state.scene.primary_basename,
        )
        for output_path in output_plan.output_paths:
            calculate_raster_overviews(
                output_path, args.overview_scales,
                log_to_console=args.log_to_console,
                scene_basename=state.scene.primary_basename,
                scene_index=getattr(args, "scene_indices", {}).get(state.scene.primary_basename, 1),
                scene_total=getattr(args, "scene_total", 1),
            )
    if mask_plan.output_paths:
        _register_step_outputs(state, "cloud_mask_mask", mask_plan.output_paths)
        state.cloud_mask_path = mask_plan.output_paths[0]
    state.current_step = "cloud_mask"
    return state.current_files


@_logged_scene_step("alignment")
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
    if _prepare_step_outputs(
        plan.output_paths,
        input_paths=plan.input_paths + [args.alignment_fixed_image],
        args=args,
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
                delete_temp_dir=args.delete_temp_dir,
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
        _log(
            "Calculating overviews for step alignment",
            enabled=args.log_to_console,
            step="overviews",
            scene_basename=state.scene.primary_basename,
        )
        for output_path in plan.output_paths:
            calculate_raster_overviews(
                output_path, args.overview_scales,
                log_to_console=args.log_to_console,
                scene_basename=state.scene.primary_basename,
                scene_index=getattr(args, "scene_indices", {}).get(state.scene.primary_basename, 1),
                scene_total=getattr(args, "scene_total", 1),
            )
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
    if args.run_spectralmatch:
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


def _configured_cleanup_inputs(args: argparse.Namespace) -> List[str]:
    """Collect caller-owned file inputs, including nested SpectralMatch references."""
    def configured_files(value):
        if isinstance(value, str):
            path = os.path.expanduser(value.removeprefix("file:"))
            if os.path.isfile(path):
                yield path
        elif isinstance(value, dict):
            for key, child in value.items():
                if not str(key).startswith("_"):
                    yield from configured_files(child)
        elif isinstance(value, (list, tuple)):
            for child in value:
                yield from configured_files(child)

    return _dedupe_paths([
        *configured_files(vars(args)),
        *configured_files(_parse_json_dict(args.spectralmatch_kwargs_json)),
        *configured_files(_normalize_group_by_basename_spec(args.group_by_basename)),
    ])


def _cleanup_workflow_temp_dirs(
    states: List[SceneWorkflowState], args: argparse.Namespace, saved_aggregate_outputs: List[str],
) -> None:
    """Remove workflow temp directories after success, independently of per-file cleanup."""
    if not args.delete_temp_dir:
        return
    protected = list(getattr(args, "_cleanup_source_files", []))
    protected.extend(_configured_cleanup_inputs(args))
    protected.extend(saved_aggregate_outputs)
    temp_roots = set()
    for state in states:
        temp_roots.add(os.path.realpath(state.step_dirs["temp_root"]))
        protected.extend(state.source_files)
        protected.append(state.step_dirs["output_root"])
        # Explicit non-temp save targets can live within a configured temp root.
        for step, path in state.step_dirs.items():
            if hasattr(args, f"save_{step}") and not _is_temp_save_value(getattr(args, f"save_{step}")):
                protected.append(path)
    locations = set()
    for path in protected:
        locations.add(os.path.realpath(path))
        # Preserve the source link itself as well as its resolved target.
        locations.add(os.path.join(os.path.realpath(os.path.dirname(path)), os.path.basename(path)))
    reserved_roots = {os.path.realpath(path) for path in (
        os.sep, os.path.expanduser("~"), os.getcwd(), tempfile.gettempdir(), "/tmp", "/var/tmp",
    )}
    locations.update(reserved_roots)
    for root in sorted(temp_roots):
        if not os.path.isdir(root):
            continue
        if root in reserved_roots or any(os.path.commonpath([root, path]) == root for path in locations):
            _log(
                f"Keeping temp directory {root}: it is shared or contains protected inputs or saved outputs",
                enabled=args.log_to_console, step="temp_cleanup",
            )
            continue
        _log(f"Deleting temp directory {root}", enabled=args.log_to_console, step="temp_cleanup")
        shutil.rmtree(root)


def _files_with_sidecars(paths: List[str]) -> List[str]:
    """Name owned raster sidecars without following dataset references or scanning directories."""
    files = list(paths)
    for path in paths:
        if not _is_gdal_raster_path(path):
            continue
        stem = os.path.splitext(path)[0]
        for suffix in (".aux.xml", ".ovr", ".msk", ".msk.ovr", ".hdr", ".params.txt", ".flaash_params.txt"):
            files.extend([path + suffix, path + suffix.upper()])
        for suffix in (".hdr", ".rpb", ".rpc", ".rpc.txt", ".imd", ".til", ".xml", ".att", ".eph", ".geo", ".ste"):
            files.extend([stem + suffix, stem + suffix.upper()])
    return _dedupe_paths(files)


def _cleanup_file_keys(paths: List[str]) -> tuple[set[str], set[tuple[int, int]]]:
    """Identify protected files by both resolved path and inode (including hard links)."""
    real_paths = set()
    inodes = set()
    for path in paths:
        real_paths.add(os.path.realpath(path))
        try:
            stat = os.stat(path)
            inodes.add((stat.st_dev, stat.st_ino))
        except FileNotFoundError:
            pass
    return real_paths, inodes


def _cleanup_completed_scene_temp_steps(
    state: SceneWorkflowState,
    args: argparse.Namespace,
    *,
    aggregate_outputs: List[str] | None = None,
    spectralmatch_outputs: List[str] | None = None,
) -> None:
    """Delete a scene's temp products only after all its persistent outputs are complete.

    Planned paths include leftovers from prior runs and partially populated temp
    steps. None means aggregate consumers may still need the final scene raster;
    a list is supplied only after all enabled aggregate steps finish successfully.
    """
    if not getattr(args, "delete_temp_steps_proactively", False):
        return
    expected = state.cleanup_step_outputs or _discovered_step_outputs(state, args)
    outputs = {**expected, **state.scene.step_outputs}
    # Stage-tagged source inputs are registered before file_source has copied them.
    outputs["file_source"] = expected.get("file_source", [])
    required, temporary = [], []
    for step in ["file_source", "fetch_atmosphere", *RASTER_STEP_ORDER[1:]]:
        if not getattr(args, f"run_{step}"):
            continue
        paths = _scene_step_expected_outputs(outputs, step)
        if _is_temp_save_value(getattr(args, f"save_{step}")):
            # Include planned files even if the step did not run in this invocation.
            temporary.extend(_scene_step_expected_outputs(expected, step))
            temporary.extend(paths)
        else:
            required.extend(paths)
    required.extend(aggregate_outputs or [])
    if required and not _existing_outputs_are_reusable(
        _dedupe_paths(required), check_validity=args.skip_existing_check_validity,
        validity_check_grid_size=args.validity_check_grid_size,
        log_to_console=args.log_to_console, step="temp_cleanup",
        scene_basename=state.scene.primary_basename,
    ):
        return
    if not required and not spectralmatch_outputs:
        return

    # SpectralMatch outputs are trusted only after its pipeline returns successfully.
    protected = list(state.source_files) + required + (spectralmatch_outputs or []) + _configured_cleanup_inputs(args)
    if args.dem_file_path not in (None, "", "online"):
        protected.append(state.dem_file_path or resolve_relative_to_input(
            args.dem_file_path, os.path.dirname(state.scene.mul_image.tif_file),
        ))
    final_paths = list(expected.get("final_raster", []))
    if state.current_step != "file_source":
        final_paths.extend(state.current_files)
    if aggregate_outputs is None and (args.run_seamline_metadata or args.run_spectralmatch):
        protected.extend(final_paths)
    temporary.extend(os.path.splitext(path)[0] + "_metadata.json" for path in final_paths if path in temporary)
    temporary.append(os.path.join(state.step_dirs["scene_work"], f"{state.scene.mul_image.basename}.gpkg"))
    if args.dem_file_path == "online":
        temporary.append(os.path.join(state.step_dirs["temp_root"], "dem", f"{state.scene.mul_image.basename}_dem.tif"))

    protected_paths, protected_inodes = _cleanup_file_keys(_files_with_sidecars(protected))
    shared_paths, shared_inodes = getattr(args, "_cleanup_source_keys", (set(), set()))
    # Protect companions of aliased inputs as well as the primary raster itself.
    for path in _dedupe_paths(temporary):
        if not os.path.isfile(path):
            continue
        real_path = os.path.realpath(path)
        stat = os.stat(path)
        inode = (stat.st_dev, stat.st_ino)
        if (real_path in protected_paths or real_path in shared_paths
                or inode in protected_inodes or inode in shared_inodes):
            protected_paths.update(os.path.realpath(item) for item in _files_with_sidecars([path]))
    to_delete = []
    spectralmatch_roots = [os.path.realpath(path) for path in (spectralmatch_outputs or [])]
    for path in _files_with_sidecars(temporary):
        real_path = os.path.realpath(path)
        if any(os.path.commonpath([real_path, root]) == root for root in spectralmatch_roots):
            continue
        if real_path in protected_paths or real_path in shared_paths or not os.path.isfile(path):
            continue
        stat = os.stat(path)
        inode = (stat.st_dev, stat.st_ino)
        if inode not in protected_inodes and inode not in shared_inodes:
            to_delete.append(path)
    if to_delete:
        _log(
            f"Deleting {len(to_delete)} temp files because saved scene outputs are complete",
            enabled=args.log_to_console, step="temp_cleanup",
            scene_basename=state.scene.primary_basename,
        )
        _delete_files(to_delete)


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
    if mul_image is None:
        raise ValueError("WorldView scene is missing required multispectral image for metadata reporting.")
    resolved_dem_file_path = state.dem_file_path
    if resolved_dem_file_path is None and args.dem_file_path not in (None, "", "online"):
        resolved_dem_file_path = str(args.dem_file_path)
    final_scene_path, scene_metadata_path = _final_output_paths(state, args)
    saved_output_paths = _scene_saved_output_paths(state, args)
    payload = {
        "scene": {
            "scene_id": state.scene.scene_id,
            "catalog_id": state.scene.catalog_id,
            "scene_root": state.scene.root_folder_path,
            "mul_photo_basename": mul_image.basename,
            "pan_photo_basename": pan_image.basename if pan_image else None,
            "started_utc": scene_started_utc,
            "completed_utc": datetime.utcnow().isoformat() + "Z",
        },
        "inputs": {
            "mul_imd_file": mul_image.imd_file,
            "mul_tif_file": mul_image.tif_file,
            "pan_imd_file": pan_image.imd_file if pan_image else None,
            "pan_tif_file": pan_image.tif_file if pan_image else None,
            "dem_file_path": resolved_dem_file_path,
            "dem_file_path_requested": args.dem_file_path,
        },
        "standardized_metadata": {
            "mul": mul_image.standardized_metadata.to_dict() if mul_image.standardized_metadata else None,
            "pan": pan_image.standardized_metadata.to_dict() if pan_image and pan_image.standardized_metadata else None,
        },
        "workflow": {
            "run_from_existing": args.run_from_existing,
            "run_from_existing_check_validity": args.run_from_existing_check_validity,
            "skip_existing_check_validity": args.skip_existing_check_validity,
            "delete_temp_steps_proactively": args.delete_temp_steps_proactively,
            "delete_temp_dir": args.delete_temp_dir,
            "run_file_source": args.run_file_source,
            "run_fetch_atmosphere": args.run_fetch_atmosphere,
            "run_atmospheric_correction": args.run_atmospheric_correction,
            "run_orthorectification": args.run_orthorectification,
            "run_pansharpen": args.run_pansharpen,
            "run_cloud_mask": args.run_cloud_mask,
            "run_alignment": args.run_alignment,
            "run_seamline_metadata": args.run_seamline_metadata,
            "run_spectralmatch": args.run_spectralmatch,
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
    _log_image_start(
        scene.primary_basename, state.current_files,
        _scene_skip_required_outputs(state, args) if args.log_to_console else [],
        enabled=args.log_to_console,
    )
    cloud_cover = _scene_cloud_cover_percent(scene)
    if args.max_cloud_cover_to_process is not None and cloud_cover is not None and cloud_cover > args.max_cloud_cover_to_process:
        _write_cloud_cover_skip_report(state, args, cloud_cover=cloud_cover)
        _log(
            f"Skipping because cloud cover {cloud_cover:.2f}% exceeds max {args.max_cloud_cover_to_process:.2f}%",
            enabled=args.log_to_console,
            step="workflow",
            scene_basename=scene.primary_basename,
        )
        _cleanup_completed_scene_temp_steps(state, args)
        return state
    if args.skip_existing and _scene_final_outputs_complete(state, args):
        _log(
            f"Skipped {getattr(args, 'scene_indices', {}).get(scene.primary_basename, 1)}/{getattr(args, 'scene_total', 1)}"
            " | Reason: desired scene level outputs exist | out="
            + ", ".join(os.path.basename(path) for path in _scene_skip_required_outputs(state, args)),
            enabled=args.log_to_console,
            scene_basename=scene.primary_basename,
        )
        _mark_scene_complete_from_existing_output(state, args)
        _cleanup_completed_scene_temp_steps(state, args)
        return state

    scene_started_utc = datetime.utcnow().isoformat() + "Z"
    _run_file_source_step(state, args)
    if args.run_atmospheric_correction and RASTER_STEP_ORDER.index(state.current_step) < RASTER_STEP_ORDER.index("atmospheric_correction"):
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
        _log(
            "Wrote scene outputs: " + ", ".join(saved_output_paths),
            enabled=args.log_to_console,
            step="workflow",
            scene_basename=scene.primary_basename,
        )
    _write_scene_report(state, args, scene_started_utc=scene_started_utc)
    _cleanup_completed_scene_temp_steps(state, args)
    _log_image_completed(
        scene.primary_basename,
        getattr(args, "scene_indices", {}).get(scene.primary_basename, 1),
        getattr(args, "scene_total", 1), enabled=args.log_to_console,
    )
    return state


def _discovered_step_outputs(state: SceneWorkflowState, args: argparse.Namespace) -> Dict[str, List[str]]:
    """Inspect each stage independently, including explicitly supplied stage inputs."""
    expected = _get_expected_scene_step_outputs(state, args)
    # Inputs can begin at a later stage. Do not add earlier suffixes to that input.
    current = list(state.current_files)
    start_index = RASTER_STEP_ORDER.index(state.current_step)
    for step in RASTER_STEP_ORDER:
        supplied = state.scene.step_outputs.get(step)
        if supplied and step != "file_source":
            expected[step] = list(supplied)
        if RASTER_STEP_ORDER.index(step) <= start_index or not getattr(args, f"run_{step}", False):
            continue
        if step == "cloud_mask":
            expected["cloud_mask_mask"] = plan_step_outputs(
                current, output_dir=state.step_dirs[step], suffix=args.cloud_mask_mask_suffix,
                skip_existing=False,
            ).output_paths
        current = plan_step_outputs(
            current, output_dir=state.step_dirs[step],
            suffix=getattr(args, f"{step}_output_suffix"),
            extension=_get_atmospheric_extension(args) if step == "atmospheric_correction" else None,
            skip_existing=False,
        ).output_paths
        if step == "orthorectification" and args.existing_mul_ortho_input:
            current = [args.existing_mul_ortho_input]
        if step == "atmospheric_correction" and args.skip_flaash and args.existing_flaash_input:
            current = [args.existing_flaash_input]
        expected[step] = list(current)
    if args.existing_pan_ortho_input:
        expected["orthorectification_pan"] = [args.existing_pan_ortho_input]
    expected["final_raster"] = current
    return expected


def _count_processing_steps(scenes: List[WorldViewScene], args: argparse.Namespace) -> Dict[str, Dict[str, int]]:
    """Count saved outputs and actual pending work over the complete discovered set."""
    steps = [
        "file_source", "fetch_atmosphere", "atmospheric_correction", "orthorectification",
        "pansharpen", "cloud_mask", "alignment", "seamline_metadata", "spectralmatch",
    ]
    counts = {step: {"loaded": 0, "processing": 0} for step in steps}
    quiet_args = argparse.Namespace(**vars(args))
    quiet_args.log_to_console = False
    states = [_initialize_scene_state(scene, quiet_args) for scene in scenes]
    outputs = [_discovered_step_outputs(state, quiet_args) for state in states]

    def complete(paths, step):
        return _existing_outputs_are_reusable(
            paths, check_validity=args.run_from_existing_check_validity,
            validity_check_grid_size=args.validity_check_grid_size,
            log_to_console=False, step=step,
        )

    for state, expected in zip(states, outputs):
        cloud_cover = _scene_cloud_cover_percent(state.scene)
        excluded = (args.max_cloud_cover_to_process is not None and cloud_cover is not None
                    and cloud_cover > args.max_cloud_cover_to_process)
        skipped = excluded or (args.skip_existing and _scene_final_outputs_complete(state, quiet_args))
        for step in steps[:-2]:
            if not getattr(args, f"run_{step}"):
                continue
            loaded = complete(_scene_step_expected_outputs(expected, step), step)
            counts[step]["loaded"] += int(loaded)
            if step == "file_source":
                runs = state.current_step == "file_source"
            elif step == "fetch_atmosphere":
                runs = args.run_atmospheric_correction and _step_will_run_from(
                    state.current_step, "atmospheric_correction", args)
            else:
                runs = _step_will_run_from(state.current_step, step, args)
            if step == "orthorectification" and args.existing_mul_ortho_input:
                runs = runs and args.run_pansharpen and not args.existing_pan_ortho_input
            if step == "atmospheric_correction" and args.skip_flaash:
                runs = False
            counts[step]["processing"] += int(runs and not skipped and not (args.run_from_existing and loaded))

    if not states:
        return counts
    final_paths = [expected["final_raster"][0] for expected in outputs]
    reference = states[0]
    if args.run_seamline_metadata:
        path = reference.step_dirs["seamline_metadata"]
        existing = set()
        if os.path.exists(path):
            import geopandas as gpd
            frame = gpd.read_file(path, layer=args.seamline_metadata_layer)
            if ("image_basename" in frame and args.seamline_metadata_image_field_name in frame
                    and frame.crs is not None and frame.crs.to_epsg() == args.epsg):
                existing = set(frame["image_basename"].dropna())
        loaded = sum(os.path.basename(path) in existing for path in final_paths)
        counts["seamline_metadata"]["loaded"] = loaded
        reuse_whole = args.run_from_existing and not args.run_from_existing_check_validity and os.path.exists(path)
        counts["seamline_metadata"]["processing"] = (
            0 if reuse_whole else len(scenes) - (loaded if args.run_from_existing else 0)
        )
    if args.run_spectralmatch:
        # Every enabled run reaches SpectralMatch; its pipeline owns reuse counts.
        counts["spectralmatch"]["processing"] = len(scenes)
    return counts


def _log_processing_steps(args: argparse.Namespace, counts: Dict[str, Dict[str, int]]) -> None:
    """Show configured steps and independently discovered output counts."""
    if not args.log_to_console:
        return
    steps = [
        "file_source", "fetch_atmosphere", "atmospheric_correction",
        "orthorectification", "pansharpen", "cloud_mask", "alignment",
        "seamline_metadata", "spectralmatch",
    ]
    _log("Processing steps:", enabled=True)
    for step in steps:
        enabled = getattr(args, f"run_{step}")
        message = f"{step}: {str(enabled).lower()}"
        if enabled:
            save_value = getattr(args, f"save_{step}")
            storage = "temp" if _is_temp_save_value(save_value) else "output"
            if step == "spectralmatch" and args.group_by_basename:
                groups = _normalize_group_by_basename_spec(args.group_by_basename)
                save_value = ", ".join(f"$output/{os.path.basename(name)}" for name in groups)
                storage = "output"
            message += f" | {storage}: {save_value}"
            if step == "spectralmatch":
                message += " | output reuse and validation handled by SpectralMatch"
            else:
                message += f" | loaded: {counts[step]['loaded']} | processing: {counts[step]['processing']}"
        _log(message, enabled=True)


def _run_workflow(args: argparse.Namespace) -> int:
    """Run the full WorldView preprocessing workflow.
    Args:
        args: Parsed CLI arguments.
    Returns:
        Process exit code.
    """
    _log_step_start("glob_matches", enabled=args.log_to_console, uppercase=False)
    filter_basenames = _parse_filter_basenames(args.filter_basename)
    input_files_by_stage = _collect_input_files_by_stage(args.input_file_glob)
    input_files = sorted({path for paths in input_files_by_stage.values() for path in paths})
    if not input_files:
        raise ValueError("No files matched --input-file-glob.")

    _log_step_start("load_worldview_scenes", enabled=args.log_to_console, uppercase=False)
    scenes = _load_worldview_scenes_from_stage_paths(input_files_by_stage, filter_basenames=filter_basenames)
    _log(
        f"Discovered {len(input_files)} input files across {len(scenes)} scenes",
        enabled=args.log_to_console,
        step="workflow",
    )
    if args.log_to_console:
        _log("Inspecting existing step outputs", enabled=True, step="workflow")
        _log_processing_steps(args, _count_processing_steps(scenes, args))
    processed_states = _process_scenes(scenes, args)

    saved_aggregate_outputs = []
    saved_spectralmatch_outputs = []
    seamline_metadata_output = None
    if processed_states and args.run_seamline_metadata:
        seamline_metadata_output = _run_seamline_metadata_workflow(
            processed_states,
            args=args,
            reference_state=processed_states[0],
        )
        if seamline_metadata_output:
            if not _is_temp_save_value(args.save_seamline_metadata):
                saved_aggregate_outputs.append(seamline_metadata_output)
            _log(
                f"Wrote seamline metadata {seamline_metadata_output}",
                enabled=args.log_to_console,
                step="seamline_metadata",
            )
        _apply_weighted_seamline_metadata_defaults(args, seamline_metadata_output)

    if processed_states and args.run_spectralmatch:
        _log(
            f"Preparing grouped SpectralMatch for {len([state for state in processed_states if state.current_files])} scene outputs",
            enabled=args.log_to_console,
            step="workflow",
        )
        spectralmatch_output = _run_spectralmatch_workflow(
            [state.current_files[0] for state in processed_states if state.current_files],
            args=args,
            reference_state=processed_states[0],
        )
        if spectralmatch_output:
            if (args.group_by_basename
                    or _build_spectralmatch_kwargs(args).get("shared_output_image_path")
                    or not _is_temp_save_value(args.save_spectralmatch)):
                saved_spectralmatch_outputs.append(spectralmatch_output)
            _log(
                f"Wrote SpectralMatch output {spectralmatch_output}",
                enabled=args.log_to_console,
                step="workflow",
            )

    if args.delete_temp_steps_proactively and (args.run_seamline_metadata or args.run_spectralmatch):
        for state in processed_states:
            _cleanup_completed_scene_temp_steps(
                state, args, aggregate_outputs=saved_aggregate_outputs,
                spectralmatch_outputs=saved_spectralmatch_outputs,
            )

    _cleanup_workflow_temp_dirs(processed_states, args, saved_aggregate_outputs + saved_spectralmatch_outputs)
    _log("All processing complete", enabled=args.log_to_console, step="workflow")
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
        help="Glob used to find raw source files. YAML configs should use one-key stage dictionaries.",
    )
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
    parser.add_argument("--concurrent-processing-backend", choices=["process_pool", "dask"], default="process_pool")
    parser.add_argument("--dask-scheduler-file")
    parser.add_argument("--dask-scheduler-address")
    parser.add_argument("--overview-scales", nargs="+")
    parser.add_argument("--temp-dir")
    parser.add_argument(
        "--delete-temp-dir", action=argparse.BooleanOptionalAction, default=True,
        help="Delete workflow temp directories after successful completion; independent of proactive per-file cleanup.",
    )
    parser.add_argument(
        "--delete-temp-steps-proactively", action=argparse.BooleanOptionalAction, default=False,
        help="Delete individual scene temp files once non-temp outputs are complete, including skipped scenes; leaves directories in place.",
    )
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
    parser.add_argument("--run-spectralmatch", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--save-spectralmatch", default="$temp/spectralmatch_root.tif")
    parser.add_argument(
        "--calculate-overviews-spectralmatch", action=argparse.BooleanOptionalAction, default=False,
        help="Enable overviews on the last eligible match step; conflicts with any enabled match-*-build-overviews flag.",
    )
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
    parser.add_argument("--spectralmatch-method", default="spectralmatch")
    parser.add_argument("--spectralmatch-kwargs-json")
    parser.add_argument("--group-by-basename")
    parser.add_argument("--match-steps", nargs="+")
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
        ("save_spectralmatch", "$temp/spectralmatch_root.tif"),
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
    args.concurrent_processing_backend = _resolve_concurrent_processing_backend(args.concurrent_processing_backend)
    if args.concurrent_processing_backend == "dask" and args.concurrent_processing != 1:
        parser.error("--concurrent-processing must be 1 when --concurrent-processing-backend=dask.")
    if args.concurrent_processing_backend == "dask" and bool(args.dask_scheduler_file) == bool(args.dask_scheduler_address):
        parser.error("--concurrent-processing-backend=dask requires exactly one of --dask-scheduler-file or --dask-scheduler-address.")
    if args.concurrent_processing_backend != "dask" and (args.dask_scheduler_file or args.dask_scheduler_address):
        parser.error("--dask-scheduler-file/--dask-scheduler-address require --concurrent-processing-backend=dask.")
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
    spectralmatch_kwargs_json = _parse_json_dict(args.spectralmatch_kwargs_json)
    save_spectralmatch_output_is_explicit = args.save_spectralmatch not in (None, "") and (
        "save_spectralmatch" in config_defaults
        or _explicit_cli_arg_present(raw_argv, "save_spectralmatch")
    )
    match_shared_output_is_explicit = (
        getattr(args, "match_shared_output_image_path", None) is not None
        or "shared_output_image_path" in spectralmatch_kwargs_json
    )
    if save_spectralmatch_output_is_explicit and match_shared_output_is_explicit:
        parser.error(
            "Cannot set both save_spectralmatch and "
            "match_shared_output_image_path/shared_output_image_path; use only one spectralmatch output path option."
        )
    group_by_basename_spec = _normalize_group_by_basename_spec(args.group_by_basename)
    if group_by_basename_spec is not None and save_spectralmatch_output_is_explicit:
        parser.error("Cannot set save_spectralmatch when group_by_basename is set; use group keys as output filenames.")
    if group_by_basename_spec is not None and match_shared_output_is_explicit:
        parser.error("Cannot set match_shared_output_image_path/shared_output_image_path when group_by_basename is set; use group keys as output filenames.")
    try:
        _build_spectralmatch_kwargs(args)
    except ValueError as exc:
        parser.error(str(exc))
    return _run_workflow(args)


if __name__ == "__main__":
    sys.exit(main())
