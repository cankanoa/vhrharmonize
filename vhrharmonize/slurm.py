"""Prepare file-to-file Slurm staging plans for vhrharmonize providers."""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import glob as std_glob
import hashlib
import os
import re
import shlex
import subprocess
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Tuple

import yaml

from vhrharmonize.cli import worldview
from vhrharmonize.providers.worldview import WorldViewImage, WorldViewScene, load_worldview_scenes_from_tif_files


PATH_TEMPLATE_RUN_ID = "{run_id}"
STAGING_DIR = ".vhr-slurm"


def load_yaml_file(path: str) -> Dict[str, Any]:
    """Load a YAML file as a dictionary."""
    with open(path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Expected YAML mapping in {path}")
    return data


def write_yaml_file(path: str, data: Mapping[str, Any]) -> None:
    """Write a YAML mapping to disk."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(dict(data), handle, sort_keys=False)


def _write_log_file(path: str, data: Mapping[str, Any]) -> None:
    """Write the Slurm log with raw status text at the bottom."""
    ordered = dict(data)
    raw_status_text = ordered.pop("raw_status_text", None)
    if raw_status_text is not None:
        ordered["raw_status_text"] = raw_status_text
    write_yaml_file(path, ordered)


def make_run_id(now: dt.datetime | None = None) -> str:
    """Return a compact UTC run id."""
    current = now or dt.datetime.now(dt.timezone.utc)
    return current.astimezone(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def resolve_run_template(value: str, run_id: str) -> str:
    """Resolve a Slurm path template with ``{run_id}``."""
    return value.replace(PATH_TEMPLATE_RUN_ID, run_id)


def _require_config_value(config: Mapping[str, Any], key: str) -> str:
    value = config.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"slurm config requires non-empty string: {key}")
    return value


def validate_slurm_config(config: Mapping[str, Any]) -> None:
    """Validate orchestration-level Slurm config values."""
    required_keys = (
        "provider",
        "provider_config",
        "log_file",
        "ssh_host",
        "ssh_user",
        "remote_output_dir",
        "remote_log_dir",
        "remote_temp_dir",
        "remote_reference_dir",
    )
    for key in required_keys:
        _require_config_value(config, key)
    _require_config_value(config, "slurm_start_file")
    for key in ("remote_output_dir", "remote_log_dir", "remote_temp_dir"):
        value = _require_config_value(config, key)
        if PATH_TEMPLATE_RUN_ID not in value:
            raise ValueError(f"{key} must include {{run_id}}.")
    slurm_start_file = _require_config_value(config, "slurm_start_file")
    if not os.path.isfile(slurm_start_file):
        raise ValueError(f"slurm_start_file does not exist: {slurm_start_file}")
    upload_keys = config.get("provider_upload_keys", [])
    if upload_keys is not None and not isinstance(upload_keys, list):
        raise ValueError("provider_upload_keys must be a list of provider YAML keys.")


def resolve_slurm_paths(config: Mapping[str, Any], run_id: str) -> Dict[str, str]:
    """Resolve remote output/log/temp/reference directories."""
    return {
        "remote_output_dir": resolve_run_template(_require_config_value(config, "remote_output_dir"), run_id),
        "remote_log_dir": resolve_run_template(_require_config_value(config, "remote_log_dir"), run_id),
        "remote_temp_dir": resolve_run_template(_require_config_value(config, "remote_temp_dir"), run_id),
        "remote_reference_dir": resolve_run_template(_require_config_value(config, "remote_reference_dir"), run_id),
    }


def _load_worldview_args(provider_config: str) -> argparse.Namespace:
    config_defaults = worldview._normalize_config_defaults(worldview._load_worldview_yaml_config(provider_config))
    parser = worldview._build_parser()
    parser.set_defaults(**config_defaults)
    args, _ = parser.parse_known_args(["--config-yaml", provider_config])
    return args


def _discover_worldview_input_files(args: argparse.Namespace) -> List[str]:
    return worldview._collect_input_files(args.input_file_glob)


def _iter_existing(paths: Iterable[str | None]) -> Iterable[str]:
    for path in paths:
        if path and os.path.isfile(path):
            yield os.path.abspath(path)


def _worldview_image_required_files(image: WorldViewImage | None) -> List[str]:
    if image is None:
        return []
    paths = std_glob.glob(os.path.join(os.path.dirname(image.tif_file), f"{image.basename}.*"))
    if image.shp_file:
        stem = os.path.splitext(image.shp_file)[0]
        paths.extend(std_glob.glob(f"{stem}.*"))
    return list(_iter_existing(paths))


def collect_worldview_upload_input_files(
    scenes: Iterable[WorldViewScene],
    args: argparse.Namespace,
) -> Tuple[List[str], List[str], List[str]]:
    """Collect workflow source files and remote discovery TIFF inputs."""
    paths: List[str] = []
    input_tifs: List[str] = []
    file_source_files: List[str] = []
    for scene in scenes:
        state = _make_planning_state(scene, args)
        source_step = worldview._get_scene_upload_source_step(state, args)
        if source_step == "file_source":
            scene_source_files: List[str] = []
            scene_source_files.extend(_worldview_image_required_files(scene.mul_image))
            scene_source_files.extend(_worldview_image_required_files(scene.pan_image))
            paths.extend(scene_source_files)
            file_source_files.extend(scene_source_files)
            for image in scene.iter_images():
                input_tifs.append(image.tif_file)
            continue

        step_files = list(_iter_existing(worldview._get_scene_upload_source_files(state, args)))
        paths.extend(step_files)
        input_tifs.extend(path for path in step_files if os.path.splitext(path)[1].lower() in {".tif", ".tiff"})
    return sorted(set(paths)), sorted(set(input_tifs)), sorted(set(file_source_files))


def _common_parent(paths: Iterable[str]) -> str:
    normalized = [os.path.abspath(path) for path in paths]
    if not normalized:
        return os.getcwd()
    common = os.path.commonpath(normalized)
    return common if os.path.isdir(common) else os.path.dirname(common)


def _remote_file_source_path(local_path: str, *, input_root: str, remote_output_dir: str) -> str:
    rel_path = os.path.relpath(local_path, input_root)
    return os.path.join(remote_output_dir, "file_source", rel_path)


def build_input_path_map(
    input_files: Iterable[str],
    *,
    file_source_files: Iterable[str],
    remote_output_dir: str,
) -> Dict[str, str]:
    """Map local workflow input files to remote output-tree paths."""
    local_files = sorted({os.path.abspath(path) for path in input_files})
    file_source_set = {os.path.abspath(path) for path in file_source_files}
    input_root = _common_parent(local_files)
    output_map: Dict[str, str] = {}
    for local_path in local_files:
        if local_path in file_source_set:
            output_map[local_path] = _remote_file_source_path(
                local_path,
                input_root=input_root,
                remote_output_dir=remote_output_dir,
            )
        else:
            output_map[local_path] = os.path.join(remote_output_dir, os.path.basename(local_path))
    return output_map


def _hash_path(path: str) -> str:
    return hashlib.sha1(os.path.abspath(path).encode("utf-8")).hexdigest()[:10]


def build_reference_path_map(reference_files: Iterable[str], *, remote_reference_dir: str) -> Dict[str, str]:
    """Map local reference files to remote reference file paths."""
    local_files = sorted({os.path.abspath(path) for path in reference_files})
    basenames: Dict[str, List[str]] = {}
    for path in local_files:
        basenames.setdefault(os.path.basename(path), []).append(path)

    path_map: Dict[str, str] = {}
    for local_path in local_files:
        basename = os.path.basename(local_path)
        if len(basenames[basename]) > 1:
            basename = f"{_hash_path(local_path)}_{basename}"
        path_map[local_path] = os.path.join(remote_reference_dir, basename)
    return path_map


def _is_local_file_reference(value: str) -> bool:
    if value.strip().lower() == "online":
        return False
    return os.path.isfile(value)


def _collect_values_for_key(value: Any, target_key: str) -> List[Any]:
    """Collect values for a key anywhere in a nested mapping."""
    matches: List[Any] = []
    if isinstance(value, dict):
        for key, item in value.items():
            if key == target_key:
                matches.append(item)
            matches.extend(_collect_values_for_key(item, target_key))
    elif isinstance(value, list):
        for item in value:
            matches.extend(_collect_values_for_key(item, target_key))
    return matches


def _collect_group_by_basename_file_values(value: Any) -> List[str]:
    """Collect literal file: paths from group_by_basename values."""
    refs: List[str] = []
    if isinstance(value, str):
        if value.startswith("file:") and os.path.isfile(value[len("file:"):]):
            refs.append(os.path.abspath(value[len("file:"):]))
    elif isinstance(value, list):
        for item in value:
            refs.extend(_collect_group_by_basename_file_values(item))
    elif isinstance(value, dict):
        for item in value.values():
            refs.extend(_collect_group_by_basename_file_values(item))
    return refs


def _collect_simple_path_values(value: Any) -> List[str]:
    """Collect existing file paths from a configured simple path value."""
    refs: List[str] = []
    if isinstance(value, str):
        if value.startswith("file:"):
            value = value[len("file:"):]
        if _is_local_file_reference(value):
            refs.append(os.path.abspath(value))
    elif isinstance(value, list):
        for item in value:
            refs.extend(_collect_simple_path_values(item))
    return refs


def collect_provider_reference_files(
    provider_config_data: Mapping[str, Any],
    *,
    upload_keys: Iterable[str],
    exclude_paths: Iterable[str],
) -> List[str]:
    """Collect explicit provider file paths requested by upload keys."""
    excluded = {os.path.abspath(path) for path in exclude_paths}
    refs: List[str] = []
    for key in upload_keys:
        if key == "input_file_glob":
            continue
        values = _collect_values_for_key(provider_config_data, key)
        for value in values:
            if key == "group_by_basename":
                refs.extend(_collect_group_by_basename_file_values(value))
            else:
                refs.extend(_collect_simple_path_values(value))
    return sorted({path for path in refs if path not in excluded})


def _remote_save_value(save_value: Any, *, remote_kind: str) -> Any:
    if not isinstance(save_value, str):
        return save_value
    normalized, save_kind = worldview._classify_save_target(save_value, default="$temp")
    if remote_kind == "temp":
        return normalized if save_kind in {"temp_root", "temp_child"} else f"$temp/{os.path.basename(normalized)}"
    if save_kind in {"temp_root", "temp_child"}:
        return normalized
    if save_kind in {"output_root", "output_child"}:
        return normalized
    return f"$output/{os.path.basename(normalized)}"


def _remote_save_value_for_key(save_key: str, save_value: Any, *, remote_output_dir: str) -> Any:
    if save_key in {"save_seamline_metadata", "save_radiometric_normalization"} and isinstance(save_value, str):
        normalized, save_kind = worldview._classify_save_target(save_value, default="$temp")
        if save_kind in {"temp_root", "temp_child"}:
            return normalized
        return os.path.join(remote_output_dir, os.path.basename(normalized))
    return _remote_save_value(save_value, remote_kind="output")


def _set_nested_key(config: MutableMapping[str, Any], section: str, key: str, value: Any) -> None:
    section_value = config.setdefault(section, {})
    if isinstance(section_value, dict):
        section_value[key] = value
    else:
        config[key] = value


def rewrite_worldview_config_for_remote(
    provider_config_data: Mapping[str, Any],
    *,
    input_tif_paths: Iterable[str] | None,
    path_rewrites: Mapping[str, str],
    remote_output_dir: str,
    remote_temp_dir: str,
) -> Dict[str, Any]:
    """Rewrite a WorldView provider config for remote execution."""
    rewritten = copy.deepcopy(dict(provider_config_data))
    if input_tif_paths is not None:
        _set_nested_key(rewritten, "shared", "input_file_glob", list(input_tif_paths))
    _set_nested_key(rewritten, "shared", "output_dir", remote_output_dir)
    _set_nested_key(rewritten, "shared", "temp_dir", remote_temp_dir)

    for section_name, section_value in list(rewritten.items()):
        if not isinstance(section_value, dict):
            continue
        for key, value in list(section_value.items()):
            if key.startswith("save_"):
                section_value[key] = _remote_save_value_for_key(
                    key,
                    value,
                    remote_output_dir=remote_output_dir,
                )

    return _rewrite_paths_recursive(rewritten, path_rewrites)


def _rewrite_string_path(value: str, path_rewrites: Mapping[str, str]) -> str:
    if value.startswith("file:"):
        local_path = os.path.abspath(value[len("file:"):])
        if local_path in path_rewrites:
            return f"file:{path_rewrites[local_path]}"
        return value
    local_path = os.path.abspath(value)
    return path_rewrites.get(local_path, value)


def _rewrite_paths_recursive(value: Any, path_rewrites: Mapping[str, str]) -> Any:
    if isinstance(value, str):
        return _rewrite_string_path(value, path_rewrites)
    if isinstance(value, list):
        return [_rewrite_paths_recursive(item, path_rewrites) for item in value]
    if isinstance(value, dict):
        return {key: _rewrite_paths_recursive(item, path_rewrites) for key, item in value.items()}
    return value


def _make_planning_state(scene: WorldViewScene, args: argparse.Namespace) -> worldview.SceneWorkflowState:
    mul_image = worldview._require_scene_image(scene, "mul")
    return worldview.SceneWorkflowState(
        scene=scene,
        step_dirs=worldview._resolve_scene_step_dirs(args, scene),
        scene_bbox_wgs84=(0.0, 0.0, 0.0, 0.0),
        current_files=[mul_image.tif_file],
    )


def _collect_planned_non_temp_outputs(
    scenes: Iterable[WorldViewScene],
    local_args: argparse.Namespace,
    *,
    remote_output_dir: str,
) -> Dict[str, str]:
    scene_list = list(scenes)
    output_map: Dict[str, str] = {}
    for scene in scene_list:
        local_state = _make_planning_state(scene, local_args)
        for local_path in worldview._scene_skip_required_outputs(local_state, local_args):
            output_map[os.path.abspath(local_path)] = _remote_output_file_path(local_path, remote_output_dir)

    if local_args.run_seamline_metadata and not worldview._is_temp_save_value(local_args.save_seamline_metadata):
        first_scene = next(iter(scene_list), None)
        if first_scene is not None:
            local_state = _make_planning_state(first_scene, local_args)
            local_path = local_state.step_dirs["seamline_metadata"]
            output_map[os.path.abspath(local_path)] = _remote_output_file_path(local_path, remote_output_dir)

    if local_args.run_radiometric_normalization and not getattr(local_args, "group_by_basename", None):
        first_scene = next(iter(scene_list), None)
        if first_scene is not None and not worldview._is_temp_save_value(local_args.save_radiometric_normalization):
            local_state = _make_planning_state(first_scene, local_args)
            local_path = local_state.step_dirs["radiometric_normalization"]
            output_map[os.path.abspath(local_path)] = _remote_output_file_path(local_path, remote_output_dir)

    if local_args.run_radiometric_normalization and getattr(local_args, "group_by_basename", None):
        first_scene = next(iter(scene_list), None)
        if first_scene is not None:
            local_state = _make_planning_state(first_scene, local_args)
            output_map.update(
                _collect_group_by_basename_output_downloads(
                    local_args.group_by_basename,
                    local_temp_root=local_state.step_dirs["temp_root"],
                    local_output_root=local_state.step_dirs["output_root"],
                    remote_output_dir=remote_output_dir,
                )
            )

    return dict(sorted(output_map.items()))


def _remote_output_file_path(local_path: str, remote_output_dir: str) -> str:
    """Return the remote output file path for a planned local output."""
    return os.path.join(remote_output_dir, os.path.basename(local_path))


def _collect_group_by_basename_output_downloads(
    group_by_basename: Any,
    *,
    local_temp_root: str,
    local_output_root: str,
    remote_output_dir: str,
) -> Dict[str, str]:
    output_map: Dict[str, str] = {}
    group_spec = worldview._normalize_group_by_basename_spec(group_by_basename)

    def _walk(spec: Mapping[str, Any]) -> None:
        for output_key, value in spec.items():
            local_path = worldview._resolve_radiometric_group_output_path(
                output_name=output_key,
                temp_root=local_temp_root,
                output_root=local_output_root,
            )
            output_map[
                os.path.abspath(local_path)
            ] = _remote_output_file_path(local_path, remote_output_dir)
            if isinstance(value, dict):
                _walk(value)
            elif isinstance(value, list):
                for item in value:
                    if isinstance(item, dict):
                        _walk(item)

    _walk(group_spec)
    return output_map


def _remote_provider_config_path(staged_provider_config: str, *, remote_reference_dir: str) -> str:
    return os.path.join(remote_reference_dir, os.path.basename(staged_provider_config))


def _add_reference_upload(reference_uploads: Dict[str, str], local_path: str, *, remote_reference_dir: str) -> str:
    local_abs = os.path.abspath(local_path)
    basename = os.path.basename(local_abs)
    used_remote_paths = set(reference_uploads.values())
    remote_path = os.path.join(remote_reference_dir, basename)
    if remote_path in used_remote_paths:
        remote_path = os.path.join(remote_reference_dir, f"{_hash_path(local_abs)}_{basename}")
    reference_uploads[local_abs] = remote_path
    return remote_path


def prepare_slurm_plan(config_path: str, *, run_id: str | None = None) -> Dict[str, Any]:
    """Prepare a Slurm staging log and staged provider config."""
    slurm_config = load_yaml_file(config_path)
    validate_slurm_config(slurm_config)
    provider = _require_config_value(slurm_config, "provider")
    if provider != "vhr-worldview":
        raise ValueError(f"Unsupported provider: {provider}")

    resolved_run_id = run_id or make_run_id()
    paths = resolve_slurm_paths(slurm_config, resolved_run_id)
    provider_config = _require_config_value(slurm_config, "provider_config")
    provider_config_data = load_yaml_file(provider_config)
    local_args = _load_worldview_args(provider_config)
    input_tifs = _discover_worldview_input_files(local_args)
    scenes = load_worldview_scenes_from_tif_files(input_tifs, filter_basenames=local_args.filter_basename)
    upload_keys = [str(key) for key in (slurm_config.get("provider_upload_keys") or [])]
    if "input_dir" in upload_keys and "input_file_glob" not in upload_keys:
        upload_keys.append("input_file_glob")
    input_tifs_for_remote: List[str] = []
    if "input_file_glob" in upload_keys:
        (
            upload_input_files,
            input_tifs_for_remote,
            file_source_files,
        ) = collect_worldview_upload_input_files(scenes, local_args)
        input_uploads = build_input_path_map(
            upload_input_files,
            file_source_files=file_source_files,
            remote_output_dir=paths["remote_output_dir"],
        )
    else:
        input_uploads = {}

    reference_files = collect_provider_reference_files(
        provider_config_data,
        upload_keys=upload_keys,
        exclude_paths=input_uploads.keys(),
    )
    reference_uploads = build_reference_path_map(reference_files, remote_reference_dir=paths["remote_reference_dir"])
    path_rewrites = {**input_uploads, **reference_uploads}
    remote_input_tifs = (
        [path_rewrites[os.path.abspath(path)] for path in input_tifs_for_remote if os.path.abspath(path) in path_rewrites]
        if input_uploads
        else None
    )
    staged_config = os.path.join(STAGING_DIR, resolved_run_id, os.path.basename(provider_config))
    staged_config_data = rewrite_worldview_config_for_remote(
        provider_config_data,
        input_tif_paths=remote_input_tifs,
        path_rewrites=path_rewrites,
        remote_output_dir=paths["remote_output_dir"],
        remote_temp_dir=paths["remote_temp_dir"],
    )
    write_yaml_file(staged_config, staged_config_data)
    staged_config_abs = os.path.abspath(staged_config)
    remote_provider_config = _add_reference_upload(
        reference_uploads,
        staged_config_abs,
        remote_reference_dir=paths["remote_reference_dir"],
    )

    output_downloads = _collect_planned_non_temp_outputs(
        scenes,
        local_args,
        remote_output_dir=paths["remote_output_dir"],
    )

    slurm_start_file = _require_config_value(slurm_config, "slurm_start_file")
    remote_slurm_start_file = _add_reference_upload(
        reference_uploads,
        slurm_start_file,
        remote_reference_dir=paths["remote_reference_dir"],
    )

    log_data: Dict[str, Any] = {
        "run_id": resolved_run_id,
        "provider": provider,
        "provider_config": provider_config,
        "staged_provider_config": staged_config_abs,
        "remote_provider_config": remote_provider_config,
        "ssh_host": _require_config_value(slurm_config, "ssh_host"),
        "ssh_user": _require_config_value(slurm_config, "ssh_user"),
        **paths,
        "slurm_start_file": slurm_start_file,
        "remote_slurm_start_file": remote_slurm_start_file,
        "submitted_job_id": None,
        "status": "prepared",
        "uploaded_input_paths": input_uploads,
        "uploaded_reference_paths": dict(sorted(reference_uploads.items())),
        "download_output_paths": output_downloads,
        "download_log_paths": {},
        "raw_status_text": "",
    }
    _write_log_file(_require_config_value(slurm_config, "log_file"), log_data)
    return log_data


def _ssh_target(log_data: Mapping[str, Any]) -> str:
    return f"{_require_config_value(log_data, 'ssh_user')}@{_require_config_value(log_data, 'ssh_host')}"


def _run_local_command(command: List[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=check, text=True, capture_output=True)


def _run_ssh(log_data: Mapping[str, Any], remote_command: str, *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return _run_local_command(["ssh", _ssh_target(log_data), remote_command], check=check)


def _remote_quote(path: str) -> str:
    return shlex.quote(path)


def _remote_parent(path: str) -> str:
    return os.path.dirname(path.rstrip("/")) or "."


def _remote_mtime(log_data: Mapping[str, Any], remote_path: str) -> int | None:
    command = (
        f"if [ -f {_remote_quote(remote_path)} ]; then "
        f"stat -c %Y {_remote_quote(remote_path)}; "
        "else echo MISSING; fi"
    )
    result = _run_ssh(log_data, command, check=False)
    text = (result.stdout or result.stderr).strip()
    if result.returncode != 0 or text == "MISSING":
        return None
    try:
        return int(text.splitlines()[-1])
    except ValueError:
        return None


def _should_upload(log_data: Mapping[str, Any], local_path: str, remote_path: str) -> bool:
    if not os.path.isfile(local_path):
        raise FileNotFoundError(local_path)
    remote_mtime = _remote_mtime(log_data, remote_path)
    if remote_mtime is None:
        return True
    return int(os.path.getmtime(local_path)) > remote_mtime


def _scp_upload(log_data: Mapping[str, Any], local_path: str, remote_path: str) -> None:
    _run_ssh(log_data, f"mkdir -p {_remote_quote(_remote_parent(remote_path))}")
    _run_local_command(["scp", "-p", local_path, f"{_ssh_target(log_data)}:{remote_path}"])


def _scp_download(log_data: Mapping[str, Any], remote_path: str, local_path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(local_path)) or ".", exist_ok=True)
    _run_local_command(["scp", "-p", f"{_ssh_target(log_data)}:{remote_path}", local_path])


def _iter_upload_maps(log_data: Mapping[str, Any]) -> Iterable[Tuple[str, str, str]]:
    for section in ("uploaded_input_paths", "uploaded_reference_paths"):
        mapping = log_data.get(section) or {}
        if not isinstance(mapping, dict):
            raise ValueError(f"{section} must be a mapping of local file to remote file.")
        for local_path, remote_path in mapping.items():
            yield section, str(local_path), str(remote_path)


def upload_required_files(log_data: Mapping[str, Any]) -> Dict[str, Dict[str, str]]:
    """Upload missing or stale files listed in the log."""
    results: Dict[str, Dict[str, str]] = {}
    for section, local_path, remote_path in _iter_upload_maps(log_data):
        label = f"{section}: {local_path} -> {remote_path}"
        try:
            if _should_upload(log_data, local_path, remote_path):
                print(f"uploading {label}")
                _scp_upload(log_data, local_path, remote_path)
                results[local_path] = {"remote_path": remote_path, "status": "uploaded"}
            else:
                print(f"already current {label}")
                results[local_path] = {"remote_path": remote_path, "status": "current"}
        except Exception as exc:
            print(f"upload error {label}")
            print(exc)
            results[local_path] = {"remote_path": remote_path, "status": "error", "error": str(exc)}
            raise
    return results


def _status_command(job_id: str) -> str:
    quoted_job = shlex.quote(job_id)
    return (
        f"squeue -j {quoted_job} -o '%i %T %M %D %R' || true; "
        "echo '---'; "
        f"sacct -j {quoted_job} --format=JobID,State,Elapsed,ExitCode,NodeList%30 -P || true"
    )


def fetch_status_text(log_data: Mapping[str, Any]) -> str:
    """Fetch raw Slurm status text for the submitted job."""
    job_id = str(log_data.get("submitted_job_id") or "").strip()
    if not job_id or job_id == "None":
        return "No submitted_job_id in log."
    result = _run_ssh(log_data, _status_command(job_id), check=False)
    return (result.stdout or "") + (result.stderr or "")


def _status_from_text(raw_status_text: str) -> str:
    text = raw_status_text.upper()
    if "COMPLETED" in text:
        return "completed"
    if any(value in text for value in ("FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY")):
        return "failed"
    if any(value in text for value in ("RUNNING", "PENDING", "CONFIGURING", "COMPLETING")):
        return "running"
    return "unknown"


def update_status_log(config_path: str) -> Dict[str, Any]:
    """Update job status fields in a Slurm log YAML."""
    log_data = load_yaml_file(config_path)
    raw_status_text = fetch_status_text(log_data)
    log_data["raw_status_text"] = raw_status_text
    log_data["status"] = _status_from_text(raw_status_text)
    _write_log_file(config_path, log_data)
    print(raw_status_text)
    return log_data


def start_slurm_job(config_path: str) -> Dict[str, Any]:
    """Upload staged files, submit the Slurm job, and update the log."""
    log_data = load_yaml_file(config_path)
    upload_results = upload_required_files(log_data)
    log_data["upload_results"] = upload_results

    remote_start_file = _require_config_value(log_data, "remote_slurm_start_file")
    remote_provider_config = _require_config_value(log_data, "remote_provider_config")
    remote_command = (
        f"cd {_remote_quote(_remote_parent(remote_start_file))} && "
        f"sbatch {_remote_quote(remote_start_file)} {_remote_quote(remote_provider_config)}"
    )
    result = _run_ssh(log_data, remote_command, check=True)
    start_output = (result.stdout or "") + (result.stderr or "")
    print(start_output)
    match = re.search(r"Submitted batch job\s+(\d+)", start_output)
    if not match:
        raise RuntimeError(f"Could not parse sbatch job id from output:\n{start_output}")

    log_data["submitted_job_id"] = match.group(1)
    log_data["status"] = "submitted"
    log_data["raw_start_output"] = start_output
    raw_status_text = fetch_status_text(log_data)
    log_data["raw_status_text"] = raw_status_text
    log_data["status"] = _status_from_text(raw_status_text)
    _write_log_file(config_path, log_data)
    print(raw_status_text)
    return log_data


def download_slurm_outputs(config_path: str) -> None:
    """Download files listed in download sections of the Slurm log."""
    log_data = load_yaml_file(config_path)
    for section in ("download_output_paths", "download_log_paths"):
        mapping = log_data.get(section) or {}
        if not isinstance(mapping, dict):
            print(f"{section} is not a mapping")
            continue
        for local_path, remote_path in mapping.items():
            print(f"downloading {remote_path} -> {local_path}")
            try:
                _scp_download(log_data, str(remote_path), str(local_path))
                print(f"downloaded {local_path}")
            except Exception as exc:
                print(f"download error {remote_path} -> {local_path}")
                print(exc)


def _build_cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="vhr-slurm", description="Prepare and manage vhrharmonize Slurm jobs.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare", help="Prepare local staged files and log YAML.")
    prepare_parser.add_argument("--config", required=True, help="Path to Slurm config YAML.")
    prepare_parser.add_argument("--run-id", help="Optional run id. Defaults to current UTC timestamp.")

    start_parser = subparsers.add_parser("start", help="Upload files and submit the Slurm job.")
    start_parser.add_argument("--config", required=True, help="Path to Slurm log YAML.")

    status_parser = subparsers.add_parser("status", help="Refresh Slurm job status in the log YAML.")
    status_parser.add_argument("--config", required=True, help="Path to Slurm log YAML.")

    download_parser = subparsers.add_parser("download", help="Download files listed in the log YAML.")
    download_parser.add_argument("--config", required=True, help="Path to Slurm log YAML.")
    return parser


def main(argv: List[str] | None = None) -> int:
    """CLI entrypoint for Slurm staging and job control."""
    parser = _build_cli_parser()
    args = parser.parse_args(argv)
    if args.command == "prepare":
        log_data = prepare_slurm_plan(args.config, run_id=args.run_id)
        print(f"Wrote Slurm log {load_yaml_file(args.config)['log_file']}")
        print(f"Wrote staged provider config {log_data['staged_provider_config']}")
        return 0
    if args.command == "start":
        start_slurm_job(args.config)
        return 0
    if args.command == "status":
        update_status_log(args.config)
        return 0
    if args.command == "download":
        download_slurm_outputs(args.config)
        return 0
    parser.error(f"Unsupported command: {args.command}")
    return 2


__all__ = [
    "build_input_path_map",
    "build_reference_path_map",
    "collect_provider_reference_files",
    "collect_worldview_upload_input_files",
    "download_slurm_outputs",
    "fetch_status_text",
    "load_yaml_file",
    "make_run_id",
    "prepare_slurm_plan",
    "resolve_slurm_paths",
    "rewrite_worldview_config_for_remote",
    "start_slurm_job",
    "update_status_log",
    "upload_required_files",
    "validate_slurm_config",
    "write_yaml_file",
]
