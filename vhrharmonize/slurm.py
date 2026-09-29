"""Prepare file-to-file Slurm staging plans for vhrharmonize workflows."""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile
from typing import Any, Dict, Iterable, List, Mapping, Tuple

import yaml

from vhrharmonize.workflow.config import load_config
from vhrharmonize.progress import ProgressSnapshot


PATH_TEMPLATE_RUN_ID = "{run_id}"
SLURM_PREPARE_CONFIG_KEYS = (
    "run_id",
    "workflow_config",
    "staged_workflow_file",
    "slurm_start_file",
    "staged_slurm_start_file",
    "staged_hpc_file",
    "debug_logs",
    "enable_rsync_checksum",
    "override_download_conflict",
    "ssh_host",
    "ssh_user",
    "ssh_private_key",
    "remote_output_dir",
    "remote_log_dir",
    "remote_temp_dir",
    "remote_reference_dir",
)


def _load_yaml_file(path: str) -> Dict[str, Any]:
    """Load a YAML file as a dictionary."""
    with open(path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Expected YAML mapping in {path}")
    return data


def _write_yaml_file(path: str, data: Mapping[str, Any]) -> None:
    """Write a YAML mapping to disk."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(dict(data), handle, sort_keys=False)


LOG_HEADER_BY_KEY = {
    "run_id": "# Set task ID",
    "workflow_config": "# Editable files",
    "staged_workflow_file": "# Generated files",
    "ssh_host": "# SSH login",
    "remote_output_dir": "# Remote directories.",
    "remote_workflow_config": "# Uploaded remote control files.",
    "remote_slurm_log_templates": "# Slurm log files.",
    "submitted_job_id": "# Job status.",
    "uploaded_input_paths": "# All mappings are local file: remote file.",
    "raw_slurm_log_text": "# Raw Slurm log text.",
    "raw_status_text": "# Raw scheduler status text.",
}


def _write_sectioned_yaml_file(
    path: str,
    data: Mapping[str, Any],
    *,
    header_by_key: Mapping[str, str] | None = None,
) -> None:
    """Write a YAML mapping with optional comments before selected top-level keys."""
    lines: List[str] = []
    for key, value in data.items():
        header = (header_by_key or {}).get(key)
        if header:
            if lines:
                lines.append("")
            lines.append(header)
        rendered = yaml.safe_dump({key: value}, sort_keys=False, width=4096).rstrip()
        lines.extend(rendered.splitlines())

    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines).rstrip() + "\n")


def _write_staged_hpc_file(path: str, data: Mapping[str, Any]) -> None:
    """Write the staged HPC YAML with readable section comments."""
    ordered = dict(data)
    upload_results = ordered.pop("upload_results", None)
    raw_start_output = ordered.pop("raw_start_output", None)
    raw_slurm_log_text = ordered.pop("raw_slurm_log_text", None)
    raw_status_text = ordered.pop("raw_status_text", None)
    if upload_results is not None:
        ordered["upload_results"] = upload_results
    if raw_start_output is not None:
        ordered["raw_start_output"] = raw_start_output
    if raw_slurm_log_text is not None:
        ordered["raw_slurm_log_text"] = raw_slurm_log_text
    if raw_status_text is not None:
        ordered["raw_status_text"] = raw_status_text
    _write_sectioned_yaml_file(path, ordered, header_by_key=LOG_HEADER_BY_KEY)


def _make_run_id(now: dt.datetime | None = None) -> str:
    """Return a compact UTC run id."""
    current = now or dt.datetime.now(dt.timezone.utc)
    return current.astimezone(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _resolve_run_template(value: str, run_id: str) -> str:
    """Resolve a Slurm path template with ``{run_id}``."""
    return value.replace(PATH_TEMPLATE_RUN_ID, run_id)


def _render_run_templates(text: str, variables: Mapping[str, str]) -> str:
    """Render supported template variables in staged text files."""
    rendered = text
    for key, value in variables.items():
        rendered = rendered.replace("{" + key + "}", value)
    return rendered


def _write_staged_template_file(source_path: str, staged_path: str, variables: Mapping[str, str]) -> None:
    """Copy a text template to a staged path with supported variables rendered."""
    with open(source_path, "r", encoding="utf-8") as handle:
        rendered = _render_run_templates(handle.read(), variables)
    os.makedirs(os.path.dirname(os.path.abspath(staged_path)) or ".", exist_ok=True)
    with open(staged_path, "w", encoding="utf-8") as handle:
        handle.write(rendered)


def _require_config_value(config: Mapping[str, Any], key: str) -> str:
    value = config.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"slurm config requires non-empty string: {key}")
    return value


def _parse_bool(value: Any, *, key: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes", "y", "on"}:
            return True
        if normalized in {"false", "0", "no", "n", "off"}:
            return False
    raise ValueError(f"{key} must be true or false.")


def _download_conflict_mode(value: Any = "validate") -> str:
    """Normalize the download policy, including YAML's unquoted yes/no booleans."""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, str) and value.strip().lower() in {"no", "yes", "validate"}:
        return value.strip().lower()
    raise ValueError("override_download_conflict must be no, yes, or validate.")


def _validate_slurm_config(config: Mapping[str, Any]) -> None:
    """Validate orchestration-level Slurm config values."""
    required_keys = (
        "workflow_config",
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
    slurm_start_file = _require_config_value(config, "slurm_start_file")
    if not os.path.isfile(slurm_start_file):
        raise ValueError(f"slurm_start_file does not exist: {slurm_start_file}")
    staged_workflow_file = config.get("staged_workflow_file")
    if staged_workflow_file is not None and (
        not isinstance(staged_workflow_file, str) or not staged_workflow_file.strip()
    ):
        raise ValueError("staged_workflow_file must be a non-empty string when set.")
    staged_hpc_file = config.get("staged_hpc_file")
    if staged_hpc_file is not None and (
        not isinstance(staged_hpc_file, str) or not staged_hpc_file.strip()
    ):
        raise ValueError("staged_hpc_file must be a non-empty string when set.")
    staged_slurm_start_file = config.get("staged_slurm_start_file")
    if staged_slurm_start_file is not None and (
        not isinstance(staged_slurm_start_file, str) or not staged_slurm_start_file.strip()
    ):
        raise ValueError("staged_slurm_start_file must be a non-empty string when set.")
    configured_run_id = config.get("run_id")
    if configured_run_id is not None and (
        not isinstance(configured_run_id, str) or not configured_run_id.strip()
    ):
        raise ValueError("run_id must be a non-empty string when set.")
    if "debug_logs" in config:
        _parse_bool(config["debug_logs"], key="debug_logs")
    if "enable_rsync_checksum" in config:
        _parse_bool(config["enable_rsync_checksum"], key="enable_rsync_checksum")
    _download_conflict_mode(config.get("override_download_conflict", "validate"))
    for key in ("ssh_private_key",):
        value = config.get(key)
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise ValueError(f"{key} must be a non-empty string when set.")


def _resolve_run_id(config: Mapping[str, Any]) -> str:
    """Resolve the Slurm run id from config or timestamp."""
    configured_run_id = config.get("run_id")
    if isinstance(configured_run_id, str) and configured_run_id.strip():
        return configured_run_id.strip()
    return _make_run_id()


def _resolve_local_template_path(value: str) -> str:
    return value if os.path.isabs(value) else os.path.abspath(value)


def _resolve_staged_workflow_file(
    config: Mapping[str, Any],
    *,
    config_path: str,
    workflow_config: str,
    run_id: str,
) -> str:
    """Resolve the local staged workflow YAML path."""
    del config_path
    configured_path = config.get("staged_workflow_file")
    if isinstance(configured_path, str) and configured_path.strip():
        staged_path = _resolve_run_template(configured_path.strip(), run_id)
    else:
        workflow_abs = os.path.abspath(workflow_config)
        workflow_dir = os.path.dirname(workflow_abs)
        workflow_name = os.path.basename(workflow_abs)
        stem, extension = os.path.splitext(workflow_name)
        workflow_label = stem.rsplit(".", 1)[-1]
        staged_path = os.path.join(workflow_dir, f"{run_id}.staged.{workflow_label}{extension or '.yml'}")
    if not os.path.isabs(staged_path):
        staged_path = _resolve_local_template_path(staged_path)
    if os.path.abspath(staged_path) == os.path.abspath(workflow_config):
        raise ValueError("staged_workflow_file must not overwrite workflow_config.")
    return staged_path


def _resolve_staged_hpc_file(config: Mapping[str, Any], *, config_path: str, run_id: str) -> str:
    """Resolve the local staged HPC YAML path."""
    configured_path = config.get("staged_hpc_file")
    if isinstance(configured_path, str) and configured_path.strip():
        staged_path = _resolve_run_template(configured_path.strip(), run_id)
    else:
        config_abs = os.path.abspath(config_path)
        config_dir = os.path.dirname(config_abs)
        staged_path = os.path.join(config_dir, f"{run_id}.staged.hpc.yml")
    return _resolve_local_template_path(staged_path)


def _resolve_staged_slurm_start_file(
    config: Mapping[str, Any],
    *,
    slurm_start_file: str,
    run_id: str,
) -> str:
    """Resolve the local staged sbatch path."""
    configured_path = config.get("staged_slurm_start_file")
    if isinstance(configured_path, str) and configured_path.strip():
        staged_path = _resolve_run_template(configured_path.strip(), run_id)
    else:
        start_abs = os.path.abspath(slurm_start_file)
        start_dir = os.path.dirname(start_abs)
        staged_path = os.path.join(start_dir, f"{run_id}.staged.slurm.sbatch")
    if not os.path.isabs(staged_path):
        staged_path = _resolve_local_template_path(staged_path)
    if os.path.abspath(staged_path) == os.path.abspath(slurm_start_file):
        raise ValueError("staged_slurm_start_file must not overwrite slurm_start_file.")
    return staged_path


def _resolve_slurm_paths(config: Mapping[str, Any], run_id: str) -> Dict[str, str]:
    """Resolve remote output, log, temp, and reference directories."""
    return {
        "remote_output_dir": _resolve_run_template(_require_config_value(config, "remote_output_dir"), run_id),
        "remote_log_dir": _resolve_run_template(_require_config_value(config, "remote_log_dir"), run_id),
        "remote_temp_dir": _resolve_run_template(_require_config_value(config, "remote_temp_dir"), run_id),
        "remote_reference_dir": _resolve_run_template(_require_config_value(config, "remote_reference_dir"), run_id),
    }


def _hash_path(path: str) -> str:
    return hashlib.sha1(  # nosec B324
        os.path.abspath(path).encode("utf-8"),
        usedforsecurity=False,
    ).hexdigest()[:10]






def _parse_sbatch_log_templates(sbatch_path: str) -> Dict[str, str]:
    """Return output/error log templates declared by SBATCH flags."""
    flag_to_key = {
        "-o": "output",
        "--output": "output",
        "-e": "error",
        "--error": "error",
    }
    templates: Dict[str, str] = {}
    with open(sbatch_path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped.startswith("#SBATCH"):
                continue
            try:
                args = shlex.split(stripped[len("#SBATCH"):].strip())
            except ValueError:
                continue
            index = 0
            while index < len(args):
                arg = args[index]
                key = None
                value = None
                if arg.startswith("--output="):
                    key = "output"
                    value = arg.split("=", 1)[1]
                elif arg.startswith("--error="):
                    key = "error"
                    value = arg.split("=", 1)[1]
                elif arg in flag_to_key and index + 1 < len(args):
                    key = flag_to_key[arg]
                    value = args[index + 1]
                    index += 1
                if key and value:
                    templates[key] = value
                index += 1
    return templates


def _resolve_remote_sbatch_log_templates(
    local_templates: Mapping[str, str],
    *,
    remote_slurm_start_file: str,
) -> Dict[str, str]:
    """Resolve sbatch log templates relative to the remote sbatch file directory."""
    remote_start_dir = _remote_parent(remote_slurm_start_file)
    resolved: Dict[str, str] = {}
    for key, template in local_templates.items():
        if template.startswith("/") or template.startswith("~/"):
            remote_template = template
        else:
            remote_template = os.path.normpath(os.path.join(remote_start_dir, template))
        resolved[key] = remote_template
    return resolved


def _resolve_slurm_log_paths(slurm_data: Mapping[str, Any]) -> Dict[str, str]:
    """Resolve sbatch log templates to concrete paths once a job id is known."""
    job_id = str(slurm_data.get("submitted_job_id") or "").strip()
    templates = slurm_data.get("remote_slurm_log_templates") or slurm_data.get("remote_slurm_log_paths") or {}
    if not isinstance(templates, dict):
        return {}
    resolved: Dict[str, str] = {}
    for key, value in templates.items():
        remote_path = str(value)
        if job_id and job_id != "None":
            remote_path = remote_path.replace("%j", job_id).replace("%A", job_id)
        resolved[str(key)] = remote_path
    return resolved


def _add_reference_upload(reference_uploads: Dict[str, str], local_path: str, *, remote_reference_dir: str) -> str:
    local_abs = os.path.abspath(local_path)
    basename = os.path.basename(local_abs)
    used_remote_paths = set(reference_uploads.values())
    remote_path = os.path.join(remote_reference_dir, basename)
    if remote_path in used_remote_paths:
        remote_path = os.path.join(remote_reference_dir, f"{_hash_path(local_abs)}_{basename}")
    reference_uploads[local_abs] = remote_path
    return remote_path


def _slurm_config_with_overrides(
    config_path: str,
    overrides: Mapping[str, Any] | None = None,
) -> Dict[str, Any]:
    """Load a HPC YAML config and apply explicit caller/CLI overrides."""
    config = _load_yaml_file(config_path)
    for key, value in (overrides or {}).items():
        if value is not None:
            config[key] = value
    return config


def prepare_slurm_plan(
    config: str,
    *,
    overrides: Mapping[str, Any] | None = None,
) -> Dict[str, Any]:
    """Prepare local Slurm, workflow and HPC staging files without uploading.

    Args:
        config: Local HPC YAML file.
        overrides: Optional mapping of HPC settings overriding values in the YAML.
    """
    slurm_config = _slurm_config_with_overrides(config, overrides)
    _validate_slurm_config(slurm_config)
    from vhrharmonize.workflow.staging import stage_workflow
    resolved_run_id = _resolve_run_id(slurm_config)
    paths = _resolve_slurm_paths(slurm_config, resolved_run_id)
    staged_hpc_file = _resolve_staged_hpc_file(slurm_config, config_path=config, run_id=resolved_run_id)
    workflow_config = _require_config_value(slurm_config, "workflow_config")
    workflow_config_data = load_config(workflow_config)
    staged_config_data, input_uploads, output_downloads = stage_workflow(
        workflow_config_data, config_dir=os.path.dirname(os.path.abspath(workflow_config)),
        remote_output_dir=paths["remote_output_dir"], remote_temp_dir=paths["remote_temp_dir"],
        remote_reference_dir=paths["remote_reference_dir"],
    )
    reference_uploads = {}
    staged_config = _resolve_staged_workflow_file(
        slurm_config, config_path=config, workflow_config=workflow_config, run_id=resolved_run_id,
    )
    _write_yaml_file(staged_config, staged_config_data)
    staged_config_abs = os.path.abspath(staged_config)
    remote_workflow_config = _add_reference_upload(
        reference_uploads,
        staged_config_abs,
        remote_reference_dir=paths["remote_reference_dir"],
    )

    slurm_start_file = _require_config_value(slurm_config, "slurm_start_file")
    staged_slurm_start_file = _resolve_staged_slurm_start_file(
        slurm_config,
        slurm_start_file=slurm_start_file,
        run_id=resolved_run_id,
    )
    _write_staged_template_file(
        slurm_start_file,
        staged_slurm_start_file,
        {"run_id": resolved_run_id},
    )
    remote_slurm_start_file = _add_reference_upload(
        reference_uploads,
        staged_slurm_start_file,
        remote_reference_dir=paths["remote_reference_dir"],
    )
    remote_slurm_log_templates = _resolve_remote_sbatch_log_templates(
        _parse_sbatch_log_templates(staged_slurm_start_file),
        remote_slurm_start_file=remote_slurm_start_file,
    )

    slurm_data: Dict[str, Any] = {
        "run_id": resolved_run_id,
        "workflow_config": workflow_config,
        "slurm_start_file": slurm_start_file,
        "staged_workflow_file": staged_config_abs,
        "staged_slurm_start_file": staged_slurm_start_file,
        "staged_hpc_file": staged_hpc_file,
        "debug_logs": _parse_bool(slurm_config.get("debug_logs", False), key="debug_logs"),
        "enable_rsync_checksum": _parse_bool(
            slurm_config.get("enable_rsync_checksum", False), key="enable_rsync_checksum"
        ),
        "override_download_conflict": _download_conflict_mode(
            slurm_config.get("override_download_conflict", "validate")
        ),
        "ssh_host": _require_config_value(slurm_config, "ssh_host"),
        "ssh_user": _require_config_value(slurm_config, "ssh_user"),
        **({"ssh_private_key": str(slurm_config["ssh_private_key"])} if slurm_config.get("ssh_private_key") else {}),
        **paths,
        "remote_workflow_config": remote_workflow_config,
        "remote_slurm_start_file": remote_slurm_start_file,
        "remote_slurm_log_templates": remote_slurm_log_templates,
        "remote_slurm_log_paths": {},
        "submitted_job_id": None,
        "status": "prepared",
        "uploaded_input_paths": input_uploads,
        "uploaded_reference_paths": dict(sorted(reference_uploads.items())),
        "download_output_paths": output_downloads,
        "download_log_paths": {},
        "raw_status_text": "",
    }
    _write_staged_hpc_file(staged_hpc_file, slurm_data)
    return slurm_data


def _ssh_target(slurm_data: Mapping[str, Any]) -> str:
    return f"{_require_config_value(slurm_data, 'ssh_user')}@{_require_config_value(slurm_data, 'ssh_host')}"


def _ssh_private_key(slurm_data: Mapping[str, Any]) -> str | None:
    configured = slurm_data.get("ssh_private_key")
    if not isinstance(configured, str) or not configured.strip():
        return None
    return os.path.expanduser(configured.strip())


def _ssh_option_args(slurm_data: Mapping[str, Any], *, open_master: bool = False) -> List[str]:
    del open_master
    option_args: List[str] = []
    private_key = _ssh_private_key(slurm_data)
    if private_key:
        option_args.extend(["-i", private_key])
    return option_args


def _ssh_command(slurm_data: Mapping[str, Any], *, open_master: bool = False) -> List[str]:
    return ["ssh", *_ssh_option_args(slurm_data, open_master=open_master), _ssh_target(slurm_data)]


def _ssh_close_command(slurm_data: Mapping[str, Any]) -> List[str]:
    return ["ssh", *_ssh_option_args(slurm_data), "-O", "exit", _ssh_target(slurm_data)]


def _ssh_command_string(slurm_data: Mapping[str, Any]) -> str:
    return " ".join(shlex.quote(part) for part in ["ssh", *_ssh_option_args(slurm_data)])


def _debug_enabled(slurm_data: Mapping[str, Any]) -> bool:
    return _parse_bool(slurm_data.get("debug_logs", False), key="debug_logs")


def _debug(slurm_data: Mapping[str, Any], message: str) -> None:
    if _debug_enabled(slurm_data):
        print(f"[slurm debug] {message}", flush=True)


def _run_local_command(
    command: List[str],
    *,
    check: bool = True,
    capture_output: bool = True,
    stream_output: bool = False,
) -> subprocess.CompletedProcess[str]:
    if stream_output:
        process = subprocess.Popen(  # nosec B603
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        output_parts: List[str] = []
        if process.stdout is None:
            raise RuntimeError("Subprocess stdout pipe was not created")
        while True:
            chunk = process.stdout.read(1)
            if not chunk:
                break
            print(chunk, end="", flush=True)
            output_parts.append(chunk)
        return_code = process.wait()
        output = "".join(output_parts)
        if check and return_code:
            raise subprocess.CalledProcessError(return_code, command, output=output)
        return subprocess.CompletedProcess(command, return_code, stdout=output, stderr="")
    return subprocess.run(  # nosec B603
        command,
        check=check,
        text=True,
        capture_output=capture_output,
    )


def _run_ssh(
    slurm_data: Mapping[str, Any],
    remote_command: str,
    *,
    check: bool = True,
    capture_output: bool = True,
    stream_output: bool = False,
) -> subprocess.CompletedProcess[str]:
    return _run_local_command(
        [*_ssh_command(slurm_data), remote_command],
        check=check,
        capture_output=capture_output,
        stream_output=stream_output,
    )


def _remote_quote(path: str) -> str:
    if path == "~":
        return "~"
    if path.startswith("~/"):
        return "~/" + shlex.quote(path[2:])
    return shlex.quote(path)


def _remote_parent(path: str) -> str:
    return os.path.dirname(path.rstrip("/")) or "."


def _scp_download(slurm_data: Mapping[str, Any], remote_path: str, local_path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(local_path)) or ".", exist_ok=True)
    command = ["scp", "-p"]
    option_args = _ssh_option_args(slurm_data)
    for index in range(0, len(option_args), 2):
        option, value = option_args[index:index + 2]
        command.extend(["-o", value] if option == "-o" else [option, value])
    command.extend([f"{_ssh_target(slurm_data)}:{remote_path}", local_path])
    _run_local_command(command)


def _remote_is_directory(slurm_data: Mapping[str, Any], remote_path: str) -> bool:
    result = _run_ssh(
        slurm_data, f"test -d {_remote_quote(remote_path)}", check=False, capture_output=True
    )
    if result.returncode not in (0, 1):
        raise RuntimeError(f"Unable to inspect remote output: {remote_path}")
    return result.returncode == 0


def _rsync_download_tree(
    slurm_data: Mapping[str, Any], remote_path: str, local_path: str, mode: str
) -> None:
    """Download a complete directory, preserving relative tile and VRT paths."""
    from vhrharmonize.io import validation
    os.makedirs(local_path, exist_ok=True)
    command = ["rsync", "-a", "--itemize-changes"]
    if _ssh_option_args(slurm_data):
        command.extend(["-e", _ssh_command_string(slurm_data)])
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8") as exclusions:
        if mode == "no":
            command.append("--ignore-existing")
        else:
            # Invalid files must be replaced even if size and mtime still match.
            command.append("--ignore-times")
            if mode == "validate":
                for root, _, files in os.walk(local_path):
                    for name in files:
                        path = os.path.join(root, name)
                        try:
                            valid = validation._existing_outputs_are_reusable(
                                [path], check_validity=True, validity_check_grid_size=0,
                                log_to_console=True, step="download",
                            )
                        except Exception:
                            valid = False
                        if valid:
                            relative = os.path.relpath(path, local_path)
                            # Anchor literal file names, including rsync pattern characters.
                            pattern = re.sub(r"([\\*?\[\]])", r"\\\1", relative)
                            exclusions.write("/" + pattern + "\0")
                exclusions.flush()
                command.extend(["--from0", "--exclude-from", exclusions.name])
        command.extend([
            f"{_ssh_target(slurm_data)}:{_remote_quote(remote_path.rstrip('/') + '/')}",
            os.path.abspath(local_path).rstrip('/') + '/',
        ])
        _run_local_command(command)


def _iter_upload_maps(slurm_data: Mapping[str, Any]) -> Iterable[Tuple[str, str, str]]:
    for section in ("uploaded_input_paths", "uploaded_reference_paths"):
        mapping = slurm_data.get(section) or {}
        if not isinstance(mapping, dict):
            raise ValueError(f"{section} must be a mapping of local file to remote file.")
        for local_path, remote_path in mapping.items():
            yield section, str(local_path), str(remote_path)


def _upload_root_for_section(slurm_data: Mapping[str, Any], section: str) -> str:
    if section == "uploaded_input_paths":
        return _require_config_value(slurm_data, "remote_output_dir")
    if section == "uploaded_reference_paths":
        return _require_config_value(slurm_data, "remote_reference_dir")
    raise ValueError(f"Unsupported upload section: {section}")


def _remote_relative_to_root(remote_path: str, remote_root: str) -> str | None:
    root = remote_root.rstrip("/")
    if remote_path == root:
        return os.path.basename(remote_path)
    prefix = f"{root}/"
    if remote_path.startswith(prefix):
        return remote_path[len(prefix):]
    return None


def _safe_remote_relative_path(path: str) -> str:
    normalized = os.path.normpath(path).lstrip("/")
    if normalized in {"", "."} or normalized == ".." or normalized.startswith("../") or "/../" in normalized:
        raise ValueError(f"Unsafe remote relative path: {path}")
    return normalized


def _stage_upload_tree(upload_items: Iterable[Tuple[str, str, str, str]], stage_root: str) -> None:
    seen: Dict[str, str] = {}
    def files_and_directories():
        for section, local_path, remote_path, remote_relative_path in upload_items:
            if os.path.isdir(local_path):
                for directory, _subdirs, filenames in os.walk(local_path):
                    relative = os.path.relpath(directory, local_path)
                    staged_directory = _safe_remote_relative_path(os.path.join(remote_relative_path, relative))
                    os.makedirs(os.path.join(stage_root, staged_directory), exist_ok=True)
                    for filename in filenames:
                        yield (section, os.path.join(directory, filename), os.path.join(remote_path, relative, filename),
                               os.path.join(remote_relative_path, relative, filename))
            else:
                yield section, local_path, remote_path, remote_relative_path

    for _section, local_path, remote_path, remote_relative_path in files_and_directories():
        if not os.path.isfile(local_path):
            raise FileNotFoundError(local_path)
        safe_relative_path = _safe_remote_relative_path(remote_relative_path)
        previous_source = seen.get(safe_relative_path)
        if previous_source is not None and os.path.abspath(previous_source) != os.path.abspath(local_path):
            raise ValueError(f"Multiple local files map to one remote path: {remote_path}")
        seen[safe_relative_path] = local_path
        staged_path = os.path.join(stage_root, safe_relative_path)
        os.makedirs(os.path.dirname(staged_path), exist_ok=True)
        if os.path.lexists(staged_path):
            os.unlink(staged_path)
        os.symlink(os.path.abspath(local_path), staged_path)


def _rsync_upload_tree(
    slurm_data: Mapping[str, Any],
    *,
    stage_root: str,
    remote_root: str,
) -> subprocess.CompletedProcess[str]:
    debug = _debug_enabled(slurm_data)
    command = ["rsync", "-aL", "--itemize-changes"]
    if _parse_bool(slurm_data.get("enable_rsync_checksum", False), key="enable_rsync_checksum"):
        command.append("--checksum")
    if debug:
        command.append("--info=progress2")
    command.extend(["--rsync-path", f"mkdir -p {_remote_quote(remote_root)} && rsync"])
    if _ssh_option_args(slurm_data):
        command.extend(["-e", _ssh_command_string(slurm_data)])
    command.extend([f"{stage_root.rstrip('/')}/", f"{_ssh_target(slurm_data)}:{remote_root.rstrip('/')}/"])
    _debug(slurm_data, f"starting batched rsync to {remote_root}")
    return _run_local_command(command, capture_output=not debug, stream_output=debug)


def _group_upload_items_by_remote_root(
    slurm_data: Mapping[str, Any],
    upload_items: Iterable[Tuple[str, str, str]],
) -> Dict[str, List[Tuple[str, str, str, str]]]:
    grouped: Dict[str, List[Tuple[str, str, str, str]]] = {}
    for section, local_path, remote_path in upload_items:
        remote_root = _upload_root_for_section(slurm_data, section)
        remote_relative_path = _remote_relative_to_root(remote_path, remote_root)
        if remote_relative_path is None:
            remote_root = _remote_parent(remote_path)
            remote_relative_path = os.path.basename(remote_path)
        grouped.setdefault(remote_root, []).append((section, local_path, remote_path, remote_relative_path))
    return grouped


def _upload_required_files(slurm_data: Mapping[str, Any]) -> Dict[str, Dict[str, str]]:
    """Upload missing or stale files listed in the staged HPC YAML."""
    results: Dict[str, Dict[str, str]] = {}
    upload_items = list(_iter_upload_maps(slurm_data))
    _debug(slurm_data, f"checking {len(upload_items)} upload files")
    grouped_uploads = _group_upload_items_by_remote_root(slurm_data, upload_items)
    for remote_root, grouped_items in grouped_uploads.items():
        print(f"syncing {len(grouped_items)} files -> {remote_root}", flush=True)
        try:
            with tempfile.TemporaryDirectory(prefix="vhr-hpc-upload-") as stage_root:
                _stage_upload_tree(grouped_items, stage_root)
                result = _rsync_upload_tree(slurm_data, stage_root=stage_root, remote_root=remote_root)
            rsync_output = ((result.stdout or "") + (result.stderr or "")).strip()
            status = "rsync_complete" if _debug_enabled(slurm_data) else ("synced" if rsync_output else "current")
            print(f"{status} {len(grouped_items)} files -> {remote_root}", flush=True)
            if rsync_output:
                _debug(slurm_data, f"rsync output for {remote_root}: {rsync_output}")
            for _section, local_path, remote_path, _remote_relative_path in grouped_items:
                results[local_path] = {"remote_path": remote_path, "status": status}
        except Exception as exc:
            print(f"sync error {len(grouped_items)} files -> {remote_root}", flush=True)
            print(exc)
            for _section, local_path, remote_path, _remote_relative_path in grouped_items:
                results[local_path] = {"remote_path": remote_path, "status": "error", "error": str(exc)}
            raise
    return results


def upload_slurm_files(config: str, *, overrides: Mapping[str, Any] | None = None) -> Dict[str, Any]:
    """Prepare when needed, upload mapped files, and write upload results.

    Args:
        config: Local HPC recipe or prepared HPC YAML file.
        overrides: Optional mapping of HPC settings overriding values in the YAML.
    """
    slurm_data = _load_yaml_file(config)
    if "uploaded_input_paths" not in slurm_data and "uploaded_reference_paths" not in slurm_data:
        slurm_data = prepare_slurm_plan(config, overrides=overrides)
        output_path = _require_config_value(slurm_data, "staged_hpc_file")
    else:
        if overrides:
            slurm_data.update(overrides)
        output_path = config

    _debug(slurm_data, f"loaded upload config: {config}")
    _debug(slurm_data, f"ssh target: {_ssh_target(slurm_data)}")
    _debug(slurm_data, "step: upload files")
    upload_results = _upload_required_files(slurm_data)
    slurm_data["upload_results"] = upload_results
    slurm_data["status"] = "uploaded"
    _debug(slurm_data, "step complete: upload files")
    _write_staged_hpc_file(output_path, slurm_data)
    return slurm_data


def _status_command(job_id: str) -> str:
    return f"scontrol show job {shlex.quote(job_id)} -dd"


def _fetch_status_text(slurm_data: Mapping[str, Any]) -> str:
    """Fetch raw Slurm status text for the submitted job."""
    job_id = str(slurm_data.get("submitted_job_id") or "").strip()
    if not job_id or job_id == "None":
        return "No submitted_job_id in staged HPC YAML."
    result = _run_ssh(slurm_data, _status_command(job_id), check=False)
    return (result.stdout or "") + (result.stderr or "")


def _read_remote_slurm_log(slurm_data: Mapping[str, Any], remote_path: str) -> str:
    """Read a remote Slurm log file, returning the remote error text on failure."""
    result = _run_ssh(slurm_data, f"cat {_remote_quote(remote_path)}", check=False)
    text = (result.stdout or "") + (result.stderr or "")
    if result.returncode:
        return f"Could not read {remote_path} (exit {result.returncode}).\n{text}".rstrip()
    return text


def _fetch_slurm_log_texts(slurm_data: Mapping[str, Any]) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Fetch resolved Slurm output/error logs declared by the sbatch file."""
    job_id = str(slurm_data.get("submitted_job_id") or "").strip()
    if not job_id or job_id == "None":
        return {}, {}
    log_paths = _resolve_slurm_log_paths(slurm_data)
    log_texts: Dict[str, str] = {}
    for key in ("error", "output"):
        remote_path = log_paths.get(key)
        if not remote_path:
            continue
        try:
            log_texts[key] = _read_remote_slurm_log(slurm_data, remote_path)
        except Exception as exc:
            log_texts[key] = f"Could not read {remote_path}: {exc}"
    return log_paths, log_texts


def _print_slurm_log_texts(log_paths: Mapping[str, str], log_texts: Mapping[str, str]) -> None:
    for key, label in (("output", "Slurm output log"), ("error", "Slurm error log")):
        remote_path = log_paths.get(key)
        if not remote_path:
            continue
        print(f"\n=== {label}: {remote_path} ===")
        text = log_texts.get(key, "")
        print(text if text else "(empty)")


def _print_slurm_status_text(raw_status_text: str) -> None:
    print("\n=== Slurm status ===")
    print(raw_status_text)


def _fetch_workflow_progress(slurm_data):
    from .progress import validate_progress_snapshot

    filename = slurm_data.get("remote_workflow_config")
    if not filename or not slurm_data.get("submitted_job_id"):
        return None
    result = _run_ssh(slurm_data, f"cat {_remote_quote(filename + '.progress.json')}", check=False)
    if result.returncode:
        return None
    try:
        data = validate_progress_snapshot(json.loads(result.stdout))
    except (ValueError, TypeError):
        return None
    if data.get("job_id") and str(data["job_id"]) != str(slurm_data["submitted_job_id"]):
        return None  # A resubmitted job must never show an earlier job's progress.
    return data


def get_slurm_progress(config: str | Path | Mapping[str, Any]) -> ProgressSnapshot | None:
    """Fetch a workflow progress snapshot over SSH without printing or editing files.

    Args:
        config: Staged HPC YAML filename or its already-loaded mapping.

    Returns:
        The public progress snapshot, or None when absent, malformed, unsupported
        or from another job. SSH connection errors propagate to the caller.
    """
    data = config if isinstance(config, Mapping) else _load_yaml_file(config)
    return _fetch_workflow_progress(data)


def _print_workflow_progress(data):
    if data is None:
        return
    from .workflow.progress_terminal import TerminalProgressDisplay

    print(f"\nWorkflow progress: {data.get('status', 'unknown')} · {data.get('updated_at', '')}")
    # Print every row once; a status command never starts the live terminal UI.
    display = TerminalProgressDisplay(stream=sys.stdout)
    display.update(data)
    display.print_snapshot()


def _status_from_text(raw_status_text: str) -> str:
    match = re.search(r"\bJobState=([A-Za-z_]+)", raw_status_text)
    if not match:
        return "unknown"
    state = match.group(1).upper()
    if state in {"FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY"}:
        return "failed"
    if state in {"RUNNING", "PENDING", "CONFIGURING", "COMPLETING"}:
        return "running"
    if state == "COMPLETED":
        return "completed"
    return "unknown"


def update_status_slurm_file(config: str) -> Dict[str, Any]:
    """Update job status fields in a staged HPC YAML.

    Args:
        config: Local staged HPC YAML file.
    """
    slurm_data = _load_yaml_file(config)
    raw_status_text = _fetch_status_text(slurm_data)
    log_paths, log_texts = _fetch_slurm_log_texts(slurm_data)
    slurm_data["raw_status_text"] = raw_status_text
    slurm_data["remote_slurm_log_paths"] = log_paths
    slurm_data["raw_slurm_log_text"] = log_texts
    slurm_data["status"] = _status_from_text(raw_status_text)
    progress = get_slurm_progress(slurm_data)
    slurm_data["workflow_progress"] = progress
    _write_staged_hpc_file(config, slurm_data)
    _print_slurm_log_texts(log_paths, log_texts)
    _print_slurm_status_text(raw_status_text)
    _print_workflow_progress(progress)
    return slurm_data


def start_slurm_job(config: str) -> Dict[str, Any]:
    """Submit the Slurm job and update the staged HPC YAML.

    Args:
        config: Local staged HPC YAML file.
    """
    slurm_data = _load_yaml_file(config)
    _debug(slurm_data, f"loaded start file: {config}")
    _debug(slurm_data, f"ssh target: {_ssh_target(slurm_data)}")

    remote_start_file = _require_config_value(slurm_data, "remote_slurm_start_file")
    remote_workflow_config = _require_config_value(slurm_data, "remote_workflow_config")
    _debug(slurm_data, f"remote start file: {remote_start_file}")
    _debug(slurm_data, f"remote workflow config: {remote_workflow_config}")
    remote_command = (
        f"cd {_remote_quote(_remote_parent(remote_start_file))} && "
        f"sbatch {_remote_quote(remote_start_file)} {_remote_quote(remote_workflow_config)}"
    )
    _debug(slurm_data, f"submit command: {remote_command}")
    _debug(slurm_data, "step: submit job")
    result = _run_ssh(
        slurm_data,
        remote_command,
        check=True,
        capture_output=not _debug_enabled(slurm_data),
        stream_output=_debug_enabled(slurm_data),
    )
    start_output = (result.stdout or "") + (result.stderr or "")
    if not _debug_enabled(slurm_data):
        print(start_output)
    match = re.search(r"Submitted batch job\s+(\d+)", start_output)
    if not match:
        raise RuntimeError(f"Could not parse sbatch job id from output:\n{start_output}")

    slurm_data["submitted_job_id"] = match.group(1)
    _debug(slurm_data, f"submitted job id: {slurm_data['submitted_job_id']}")
    slurm_data["status"] = "submitted"
    slurm_data["raw_start_output"] = start_output
    slurm_data["remote_slurm_log_paths"] = _resolve_slurm_log_paths(slurm_data)
    _debug(slurm_data, "step: fetch status")
    raw_status_text = _fetch_status_text(slurm_data)
    slurm_data["raw_status_text"] = raw_status_text
    slurm_data["status"] = _status_from_text(raw_status_text)
    _debug(slurm_data, f"status after submit: {slurm_data['status']}")
    _debug(slurm_data, "step: write start file")
    _write_staged_hpc_file(config, slurm_data)
    print(raw_status_text)
    return slurm_data


def _require_submitted_job_id(slurm_data: Mapping[str, Any]) -> str:
    job_id = str(slurm_data.get("submitted_job_id") or "").strip()
    if not job_id or job_id == "None":
        raise ValueError("staged HPC YAML has no submitted_job_id.")
    return job_id


def stop_slurm_job(config: str) -> subprocess.CompletedProcess[str]:
    """Cancel the submitted Slurm job listed in the staged HPC YAML.

    Args:
        config: Local staged HPC YAML file.
    """
    slurm_data = _load_yaml_file(config)
    job_id = _require_submitted_job_id(slurm_data)
    result = _run_ssh(slurm_data, f"scancel {shlex.quote(job_id)}", check=True)
    output = (result.stdout or "") + (result.stderr or "")
    if output:
        print(output)
    return result


def close_hpc_connection(config: str) -> subprocess.CompletedProcess[str]:
    """Close the SSH multiplex master for the staged HPC YAML target.

    Args:
        config: Local staged HPC YAML file.
    """
    slurm_data = _load_yaml_file(config)
    result = _run_local_command(_ssh_close_command(slurm_data), check=False)
    output = (result.stdout or "") + (result.stderr or "")
    if output:
        print(output)
    return result


def download_slurm_outputs(config: str) -> None:
    """Download files or complete directories declared in the staged HPC YAML.

    Args:
        config: Local staged HPC YAML file.
    """
    slurm_data = _load_yaml_file(config)
    mode = _download_conflict_mode(slurm_data.get("override_download_conflict", "validate"))
    from vhrharmonize.io import validation
    for section in ("download_output_paths", "download_log_paths"):
        mapping = slurm_data.get(section) or {}
        if not isinstance(mapping, dict):
            raise ValueError(f"{section} must be a mapping")
        for local_path, remote_path in mapping.items():
            try:
                if _remote_is_directory(slurm_data, str(remote_path)):
                    _rsync_download_tree(slurm_data, str(remote_path), str(local_path), mode)
                    print(f"downloaded directory {local_path}")
                    continue
            except Exception as exc:
                raise RuntimeError(f"Download failed: {remote_path} -> {local_path}") from exc
            if mode != "yes":
                try:
                    reusable = validation._existing_outputs_are_reusable(
                        [str(local_path)],
                        check_validity=mode == "validate",
                        validity_check_grid_size=0,
                        log_to_console=True,
                        step="download",
                    )
                except Exception as exc:
                    # GDAL can raise on corrupt local files before opening a dataset.
                    print(f"local validation failed {local_path}: {exc}")
                    reusable = False
                if reusable:
                    print(f"skipping existing local file {local_path}")
                    continue
            print(f"downloading {remote_path} -> {local_path}")
            try:
                _scp_download(slurm_data, str(remote_path), str(local_path))
                print(f"downloaded {local_path}")
            except Exception as exc:
                raise RuntimeError(f"Download failed: {remote_path} -> {local_path}") from exc


__all__ = [
    "close_hpc_connection",
    "download_slurm_outputs",
    "prepare_slurm_plan",
    "start_slurm_job",
    "stop_slurm_job",
    "update_status_slurm_file",
    "upload_slurm_files",
]
