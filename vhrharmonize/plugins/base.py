"""Small plugin contract: receive explicit paths, parameters and shared settings."""

from __future__ import annotations
from dataclasses import asdict, is_dataclass
from importlib import import_module
import inspect
import json


INPUT_PATH_FEATURES = (
    "input_dependency_paths",
    "input_existence_check_paths",
    "input_protection_paths",
    "input_hpc_staging_paths",
)
OUTPUT_PATH_FEATURES = (
    "output_path_resolution_paths",
    "output_dependency_paths",
    "output_target_paths",
    "output_parent_creation_paths",
    "output_collision_check_paths",
    "output_reuse_paths",
    "output_validation_paths",
    "output_invalid_removal_paths",
    "output_overview_calculation_paths",
    "output_temporary_cleanup_paths",
    "output_hpc_staging_paths",
    "output_hpc_download_paths",
)


def file_parameter_names(features, direction):
    """Collect argument names for routing, without enabling any other feature."""
    return set().union(
        *(names for key, names in features.items() if key.startswith(direction + "_"))
    )


def json_value(value):
    if is_dataclass(value):
        value = asdict(value)
    elif hasattr(value, "to_dict"):
        value = value.to_dict()
    return json.loads(json.dumps(value, default=str, allow_nan=False))


class FunctionPlugin:
    scope = "var"
    var_records_return = None  # Returned field (or "$" for the whole result) containing scenes.
    var_records_mode = "replace"  # replace | merge (preserve existing fields, keyed by scene ID).
    constant_values_return = None
    temporary_directory_context_paths = ()
    output_directory_context_paths = ()
    var_id_return = None
    var_path_return = None  # Primary path within each returned scene; used by core:satisfies and HPC mappings.
    source_file_protection_paths_return = None
    directory_parameters = {}  # Declared directory context location -> function parameter accepting that directory.
    discovery_input_parameter = None  # Optional parameter name for discovery upload labels.
    # Each set contains function argument names, not filenames. Opt in separately.
    input_dependency_paths = frozenset()
    input_existence_check_paths = frozenset()
    input_protection_paths = frozenset()
    input_hpc_staging_paths = frozenset()

    output_path_resolution_paths = frozenset()
    output_dependency_paths = frozenset()
    output_target_paths = frozenset()  # Core fills this from core:require_outputs; adapters do not choose targets.
    output_parent_creation_paths = frozenset()
    output_collision_check_paths = frozenset()
    output_reuse_paths = frozenset()
    output_validation_paths = frozenset()
    output_invalid_removal_paths = frozenset()
    output_overview_calculation_paths = frozenset()
    output_temporary_cleanup_paths = frozenset()
    output_hpc_staging_paths = frozenset()
    output_hpc_download_paths = frozenset()

    def file_features(self):
        """Validate independent file-feature declarations without loading the function."""
        if self.scope not in {"var", "aggregate"}:
            raise ValueError("Plugin scope must be var or aggregate")
        for name in (
            "var_records_return",
            "var_id_return",
            "var_path_return",
            "constant_values_return",
            "source_file_protection_paths_return",
        ):
            selector = getattr(self, name)
            if selector is not None and (not isinstance(selector, str) or not selector):
                raise ValueError(f"{name} must be a returned-field selector or None")
        if self.var_records_mode not in {"replace", "merge"}:
            raise ValueError("var_records_mode must be replace or merge")
        if self.var_records_mode == "merge" and not self.var_id_return:
            raise ValueError("Merging scenes requires var_id_return for stable scene identities")
        for name in ("temporary_directory_context_paths", "output_directory_context_paths"):
            selectors = getattr(self, name)
            if not isinstance(selectors, (tuple, list, set, frozenset)) or any(
                not isinstance(selector, str)
                or selector.split(".", 1)[0] not in {"const", "var"}
                or "." not in selector
                or not all(p.isidentifier() for p in selector.split(".")[1:])
                for selector in selectors
            ):
                raise ValueError(f"{name} must contain const.field or var.field JSON locations")
        locations = {*self.temporary_directory_context_paths, *self.output_directory_context_paths}
        if not isinstance(self.directory_parameters, dict) or any(
            location not in locations or not isinstance(parameter, str) or not parameter.isidentifier()
            for location, parameter in self.directory_parameters.items()
        ):
            raise ValueError("directory_parameters must map declared directory context locations to parameter names")
        if self.discovery_input_parameter is not None and (
            not isinstance(self.discovery_input_parameter, str) or not self.discovery_input_parameter.isidentifier()
        ):
            raise ValueError("discovery_input_parameter must be a parameter name or None")
        obsolete = {
            "scene_records_return",
            "scene_records_mode",
            "scene_id_return",
            "scene_path_return",
            "path_parameters",
            "output_parameters",
            "manages_reuse",
            "manages_overviews",
            "input_path_resolution_paths",
            "scene_records_during_planning",
            "scene_base_dir_return",
            "scene_input_paths_return",
            "scene_metadata_output_return",
            "output_context_checkpoint_paths",
        }
        removed = sorted(name for name in obsolete if hasattr(self, name))
        if removed:
            raise ValueError(
                "Removed plugin declarations: "
                + ", ".join(removed)
                + ". Use the explicit file features and context directory locations."
            )
        features = {}
        for name in (*INPUT_PATH_FEATURES, *OUTPUT_PATH_FEATURES):
            value = getattr(self, name)
            if not isinstance(value, (set, frozenset, tuple, list)) or any(
                not isinstance(item, str) or not item.isidentifier() for item in value
            ):
                raise ValueError(f"{name} must select function parameter names")
            features[name] = frozenset(value)
        overlap = file_parameter_names(features, "input") & file_parameter_names(features, "output")
        if overlap:
            raise ValueError(
                f"File parameters cannot be both inputs and outputs: {sorted(overlap)}"
            )
        return features

    target = ""
    aliases = {}
    options = set()  # Additional supported kwargs for functions with **kwargs.
    passthrough = False

    def stage_settings(self, *, settings, params, returned, path_mappings, file_paths,
                       discovery_paths, config_dir):
        """Return param: overrides for discovery reruns after HPC path relocation.

        Called only for discovery invocations that actually ran during planning.
        Inputs are detached copies: rewritten YAML settings, resolved local params,
        the function result, local-to-remote path mappings, individually mapped file
        paths, all discovered source paths, and the local YAML directory. The default
        needs no overrides; context-loaded steps do not rerun discovery or this hook.
        """
        return {}

    def function(self):
        module, name = self.target.split(":")
        return getattr(import_module(module), name)

    def arguments(self, function, params, shared):
        signature = inspect.signature(function)
        from vhrharmonize.parameters import function_parameters

        accepted = set(function_parameters(function)) | self.options

        def rename(values):
            return {self.aliases.get(key, key): value for key, value in values.items()}

        # Every plugin sees all shared settings; only supported function arguments
        # are forwarded. Explicit per-step settings override shared values.
        inherited = rename(shared)
        kwargs = {key: value for key, value in inherited.items() if key in accepted}
        explicit = rename(params)
        unknown = set(explicit) - accepted
        if unknown and not self.passthrough:
            raise ValueError(
                f"Unsupported {type(self).__name__} options: {', '.join(sorted(unknown))}"
            )
        kwargs.update(explicit)
        signature.bind(**kwargs)
        return kwargs

    def run(self, *, params, shared):
        from vhrharmonize.io.progress import reports_progress

        function = self.function()
        function = reports_progress(function)
        return json_value(function(**self.arguments(function, params, shared)))
