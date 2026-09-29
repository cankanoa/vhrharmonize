"""Named workflow steps with explicit plugin selection and prefixed settings."""

from __future__ import annotations
from copy import deepcopy
import yaml
from .registry import plugin_names


class UniqueLoader(yaml.SafeLoader):
    pass


def _mapping(loader, node, deep=False):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in result:
            raise ValueError(f"Duplicate YAML key {key!r} on line {key_node.start_mark.line + 1}")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)

STEP_CONTROLS = {
    "run",
    "scope",
    "reuse",
    "check_validity",
    "calculate_overviews",
    "required",
    "requires",
    "processing_direction",
}
DEFAULT_SHARED = {
    "protect_source_files": True,
    "output_metadata_path": None,
    "delete_final_json_first": True,
    "run_from_existing": True,
    "check_validity": True,
    "validity_check_grid_size": 2048,
    "delete_temp_dir": False,
    "delete_temp_steps_proactively": True,
    "log_to_console": True,
    "show_progress": True,
    "report_progress": False,
    "concurrent_processing": 1,
    "concurrent_processing_backend": "process_pool",
    "processing_direction": "vertical",
}
SHARED_CONTROLS = set(DEFAULT_SHARED) | {"dask_scheduler_address", "dask_scheduler_file"}
BOOLEANS = {
    "protect_source_files",
    "delete_final_json_first",
    "run",
    "reuse",
    "required",
    "check_validity",
    "calculate_overviews",
    "run_from_existing",
    "delete_temp_dir",
    "delete_temp_steps_proactively",
    "log_to_console",
    "show_progress",
    "report_progress",
}


def _validate_settings(settings, *, shared=False, context_only=False):
    if not isinstance(settings, dict):
        raise ValueError("Plugin settings and shared must be mappings")
    for key, value in settings.items():
        if key == "plugin":
            continue
        if not isinstance(key, str) or ":" not in key:
            raise ValueError(f"Setting {key!r} needs a param:, var:, const:, or core: prefix")
        kind, name = key.split(":", 1)
        if (
            kind not in {"param", "var", "const", "core"}
            or not name
            or not all(p.isidentifier() for p in name.split("."))
            or kind not in {"var", "const"}
            and "." in name
        ):
            raise ValueError(
                f"Invalid setting {key!r}; use param:name, var:name, const:name, or core:name"
            )
        if context_only and (kind == "param" or _uses_returned(value)):
            raise ValueError("Steps without a plugin cannot use param: or returned: values")
        if kind == "core":
            if name not in (SHARED_CONTROLS | {"run"} if shared else STEP_CONTROLS):
                raise ValueError(f"Unknown {'shared ' if shared else ''}core setting: {name}")
            if name in BOOLEANS and not isinstance(value, bool):
                raise ValueError(f"{key} must be a boolean")
            if (
                name == "output_metadata_path"
                and value is not None
                and (not isinstance(value, str) or not value)
            ):
                raise ValueError(
                    "core:output_metadata_path must be a path/reference string or null"
                )
            if name == "scope" and value not in {"scene", "aggregate"}:
                raise ValueError("core:scope must be scene or aggregate")
            if name == "processing_direction" and value not in ("horizontal", "vertical"):
                raise ValueError("core:processing_direction must be horizontal or vertical")
            if name == "validity_check_grid_size" and (
                not isinstance(value, int) or isinstance(value, bool) or value < 0
            ):
                raise ValueError("core:validity_check_grid_size must be a non-negative integer")


def _uses_returned(value):
    if isinstance(value, dict):
        return any(_uses_returned(v) for v in value.values())
    if isinstance(value, (list, tuple)):
        return any(_uses_returned(v) for v in value)
    return isinstance(value, str) and value.startswith("returned:")


def validate_config(data):
    if not isinstance(data, dict) or not data:
        raise ValueError("Configuration must map unique step names to settings")
    known = {*plugin_names(), "shared"}
    for name, settings in data.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Every step must have a nonempty string name")
        if not isinstance(settings, dict):
            raise ValueError(f"Step {name!r} must be one mapping, with at most one plugin")
        plugin = settings.get("plugin")
        if "plugin" in settings and (not isinstance(plugin, str) or plugin not in known):
            raise ValueError(f"Unknown plugin {plugin!r} in step {name!r}")
        _validate_settings(settings, shared=plugin == "shared", context_only=plugin is None)
    return deepcopy(data)


def steps(config):
    """Step names identify occurrences; plugin names select implementations."""
    return [
        {
            "name": name,
            "plugin": settings.get("plugin"),
            "settings": {k: v for k, v in settings.items() if k != "plugin"},
            "run": False,
            **{key[5:]: value for key, value in settings.items() if key.startswith("core:")},
        }
        for name, settings in config.items()
        if settings.get("plugin") != "shared"
    ]


def step_settings(config, step):
    return config[step["name"]]


def shared_blocks(config):
    """Core shared plugins supply workflow-wide defaults before ordinary steps."""
    return [
        settings
        for settings in config.values()
        if settings.get("plugin") == "shared" and settings.get("core:run", False)
    ]


def shared_settings(config):
    core, params, assignments = dict(DEFAULT_SHARED), {}, []
    for settings in shared_blocks(config):
        core.update(
            {k[5:]: v for k, v in settings.items() if k.startswith("core:") and k != "core:run"}
        )
        params.update({k[6:]: v for k, v in settings.items() if k.startswith("param:")})
        assignments.append({k: v for k, v in settings.items() if k.startswith(("const:", "var:"))})
    return core, params, assignments


def load_config(filename):
    with open(filename, encoding="utf-8") as handle:
        return validate_config(yaml.load(handle, Loader=UniqueLoader))
