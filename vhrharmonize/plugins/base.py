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
    "output_context_checkpoint_paths",
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
    scope = "scene"
    scene_records_return = None  # Returned field (or "$" for the whole result) containing scenes.
    scene_records_mode = "replace"  # replace | merge (preserve existing fields, keyed by scene ID).
    constant_values_return = None
    temporary_directory_context_paths = ()
    output_directory_context_paths = ()
    scene_id_return = None
    source_file_protection_paths_return = None
    # Each set contains function argument names, not filenames. Opt in separately.
    input_dependency_paths = frozenset()
    input_existence_check_paths = frozenset()
    input_protection_paths = frozenset()
    input_hpc_staging_paths = frozenset()

    output_path_resolution_paths = frozenset()
    output_dependency_paths = frozenset()
    output_target_paths = frozenset()
    output_parent_creation_paths = frozenset()
    output_collision_check_paths = frozenset()
    output_reuse_paths = frozenset()
    output_validation_paths = frozenset()
    output_invalid_removal_paths = frozenset()
    output_overview_calculation_paths = frozenset()
    output_temporary_cleanup_paths = frozenset()
    output_hpc_staging_paths = frozenset()
    output_hpc_download_paths = frozenset()
    # First supplied selected output anchors .context.json.
    output_context_checkpoint_paths = frozenset()

    def file_features(self):
        """Validate independent file-feature declarations without loading the function."""
        for name in (
            "scene_records_return",
            "scene_id_return",
            "constant_values_return",
            "source_file_protection_paths_return",
        ):
            selector = getattr(self, name)
            if selector is not None and (not isinstance(selector, str) or not selector):
                raise ValueError(f"{name} must be a returned-field selector or None")
        if self.scene_records_mode not in {"replace", "merge"}:
            raise ValueError("scene_records_mode must be replace or merge")
        if self.scene_records_mode == "merge" and not self.scene_id_return:
            raise ValueError("Merging scenes requires scene_id_return for stable scene identities")
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
        obsolete = {
            "path_parameters",
            "output_parameters",
            "manages_reuse",
            "manages_overviews",
            "input_path_resolution_paths",
            "scene_records_during_planning",
            "scene_base_dir_return",
            "scene_input_paths_return",
            "scene_metadata_output_return",
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

    def restore(self, params):
        return None

    target = ""
    aliases = {}
    options = set()  # Additional supported kwargs for functions with **kwargs.
    passthrough = False

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
