"""Plan file dependencies backward and execute ordered, prefixed plugin settings."""

from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import nullcontext
from copy import deepcopy
from dataclasses import dataclass, field
import json
import inspect
import os
from multiprocessing import get_context
from pathlib import Path
from time import monotonic, time_ns
from threading import RLock
from uuid import uuid4
import warnings

from vhrharmonize.plugins.base import (
    file_parameter_names,
    INPUT_PATH_FEATURES,
    OUTPUT_PATH_FEATURES,
)
from vhrharmonize.io.metadata import write_json
from vhrharmonize.io.logging import _log
from vhrharmonize.io.validation import _existing_output_failures
from vhrharmonize.io.workflow_utils import remove_output_files
from .timing import timed_core, timed_preflight
from .config import validate_config, steps, shared_settings
from .registry import load_plugin
from .paths import DIRECTORY_FEATURES, directory_values, directory_bindings
from .metadata import FinalMetadataWriter
from .context_io import ContextFiles, LOAD_CONTROLS, SAVE_CONTROLS
from .values import (
    Pending,
    PathResolver,
    empty_context,
    dependency_values,
    matches_reference,
    Deferred,
    available_context,
    remap_paths,
    assign,
    contains_pending,
    evaluate_settings,
    constant_settings,
    lookup,
    path,
    references,
    resolve,
    aggregate_variables,
    scene_values,
    update_var_records,
)


def _add_missing(existing, additions):
    """Build on a JSON object without overwriting existing values, including nulls/lists."""
    result = deepcopy(existing)
    for name, value in additions.items():
        if name not in result:
            result[name] = deepcopy(value)
        elif isinstance(result[name], dict) and isinstance(value, dict):
            result[name] = _add_missing(result[name], value)
    return result


def _paths(value):
    if isinstance(value, str):
        return [value]
    if isinstance(value, dict):
        return [p for v in value.values() for p in _paths(v)]
    if isinstance(value, list):
        return [p for v in value for p in _paths(v)]
    return []


def _within(filename, directory):
    if isinstance(directory, (list, tuple)):
        return any(_within(filename, root) for root in directory)
    return os.path.commonpath(
        [os.path.realpath(filename), os.path.realpath(directory)]
    ) == os.path.realpath(directory)


def _materialize(value, runtime, selector=""):
    if isinstance(value, Pending):
        try:
            return deepcopy(lookup(runtime, selector))
        except ValueError:
            return value
    if isinstance(value, dict):
        return {
            k: _materialize(v, runtime, f"{selector}.{k}" if selector else k)
            for k, v in value.items()
        }
    if isinstance(value, list):
        return [_materialize(v, runtime, f"{selector}.{i}") for i, v in enumerate(value)]
    return deepcopy(value)


def _arguments(features, params, context, settings, base_dir, *, scene_ids=None, updates=None):
    params = deepcopy(params)
    for name in file_parameter_names(features, "input") & params.keys():
        value = params[name]
        if isinstance(value, str):
            params[name] = os.path.expanduser(value)
        elif isinstance(value, list):
            params[name] = [os.path.expanduser(item) if isinstance(item, str) else item for item in value]
    resolution = features["output_path_resolution_paths"]
    for name in resolution & params.keys():
        if params[name] is None or contains_pending(params[name]):
            continue
        params[name] = path(params[name], base_dir=base_dir)
        template = settings.get("param:" + name)
        if isinstance(template, str) and template.startswith(("var:", "const:", "collect:")):
            scope, field = template.split(":", 1)
            if scope in {"var", "collect"} and scene_ids is not None:
                records = [
                    {"var": {key: values[i] for key, values in context["var"].items()}}
                    for i in range(len(scene_ids))
                ]
                update = {"var." + field: params[name]}
                update_var_records(records, update, scene_ids)
                context["var"] = aggregate_variables(records)
                if updates is not None:
                    updates.update(update)
            elif scope != "collect":
                assign(context[scope], field, params[name])
    return params


def _selected(params, names):
    return {
        name: value
        for name, value in params.items()
        if name in names and value is not None and not contains_pending(value)
    }


def _target_features(plugin, step):
    """The recipe chooses deliverables; adapters still declare file capabilities."""
    features = plugin.file_features()
    features["output_target_paths"] = frozenset()
    require_outputs = step.get("require_outputs", False)
    outputs = file_parameter_names(features, "output")
    if isinstance(require_outputs, bool):
        selected = outputs if require_outputs else set()
    else:
        selected = {value[6:] for value in ([require_outputs] if isinstance(require_outputs, str) else require_outputs)}
        if selected - outputs:
            raise ValueError(
                f"{step['name']} core:require_outputs must select declared output parameters: "
                f"{sorted(selected - outputs)}"
            )
    features["output_target_paths"] = frozenset(selected)
    return features


def _required_arguments(params, step, context, requirements, *, records=None, path_resolver=None):
    """Normalize files embedded in compound parameters using explicit requires links."""
    if not step.get("requires"):
        return params
    raw = _paths(resolve(step["requires"], context, records=records, path_resolver=path_resolver))
    return remap_paths(params, dict(zip(raw, requirements)))


@dataclass
class Node:
    index: int
    step_index: int
    step: dict
    record: int | None
    params: dict
    requirements: list
    context: dict
    pre_context: dict
    updates: dict
    dynamic_names: set
    dependencies: set[int] = field(default_factory=set)
    demanded_paths: set[str] = field(default_factory=set)
    demanded_values: set[str] = field(default_factory=set)
    value_dependencies: dict = field(default_factory=dict)
    loaded: bool = False
    needed: bool = False
    status: str = "unused"
    runtime_base: dict = field(default_factory=dict)
    parameter_overrides: dict = field(default_factory=dict)
    collection_snapshot: list = field(default_factory=list)
    collection_result: list = field(default_factory=list)
    base_dir: str = "."
    file_features: dict = field(default_factory=dict)
    directory_locations: dict = field(default_factory=dict)
    runtime_directory_context: dict | None = None
    requirements_pending: bool = False

    def file_arguments(self, *features):
        return _selected(self.params, set().union(*(self.file_features[f] for f in features)))

    def paths(self, *features):
        names = set().union(*(self.file_features[feature] for feature in features))
        return [filename for name, value in self.params.items() if name in names for filename in _paths(value)]

    @property
    def directories(self):
        context = (
            self.context
            if self.runtime_directory_context is None
            else self.runtime_directory_context
        )
        return directory_values(context, self.directory_locations, base_dir=self.base_dir)

    @property
    def directory_bindings(self):
        context = (
            self.context
            if self.runtime_directory_context is None
            else self.runtime_directory_context
        )
        return directory_bindings(context, self.directory_locations)

@dataclass
class ConstantBindings:
    settings: dict
    before: dict
    planned: dict
    resolved: dict = field(default_factory=dict)
    runtime_before: dict = field(default_factory=dict)
    records: list = field(default_factory=list)


def _execute(payload):
    plugin_name, params, shared, *reporters = payload
    if not reporters:
        return load_plugin(plugin_name).run(params=params, shared=shared)
    from vhrharmonize.io.progress import capture_messages, progress_context

    reporter = reporters[0]
    started, started_ns = monotonic(), time_ns()
    reporter.event("start", started_ns=started_ns)
    try:
        with progress_context(reporter, reporter.message), capture_messages():
            result = load_plugin(plugin_name).run(params=params, shared=shared)
    except BaseException as exc:
        reporter.event("failed", duration=monotonic() - started, started_ns=started_ns, error_type=type(exc).__name__)
        raise
    reporter.event("computed", duration=monotonic() - started, started_ns=started_ns)
    return result


class Workflow:
    @timed_core("initialization")
    def __init__(self, config, *, config_dir=".", selected_plugin=None, selected_step=None, preparing=False):
        self._timing_lock = RLock()
        self._timing_buffer = []
        self._timing_callbacks = None
        self._timing_finished = False
        self._historical_timings = {}
        self._run_id = uuid4().hex
        self._run_started_ns, self._run_started = time_ns(), monotonic()
        self.selected_plugin = selected_plugin
        self.selected_step = selected_step
        self.config = validate_config(config)
        self.steps = steps(self.config)
        self.preparing = preparing
        self.preparation_state = None
        self.preparation_index = None
        self.has_skipped_calls = any(step.get("skip_plugin_call", False) for step in self.steps)
        skipping = False
        for step in self.steps:
            if not step["run"]:
                continue
            if skipping and not step.get("skip_plugin_call", False):
                raise ValueError("Skipped plugin calls must form a suffix of the workflow")
            skipping |= step.get("skip_plugin_call", False)
        if selected_step is not None:
            matching = [i for i, step in enumerate(self.steps) if step["name"] == selected_step]
            if not matching:
                raise ValueError(f"Unknown workflow step: {selected_step}")
            self.steps = self.steps[:matching[0] + 1]
        if selected_plugin is not None:
            matching = [i for i, s in enumerate(self.steps) if s["plugin"] == selected_plugin]
            if not matching:
                raise ValueError(f"Plugin {selected_plugin!r} is not present in the configuration")
            self.steps = self.steps[: max(matching) + 1]
        core, params, assignments = shared_settings(self.config)
        self.controls = core
        self._log_core_start("workflow")
        self.shared = params
        self.config_dir = os.path.abspath(config_dir)
        self.path_resolver = PathResolver(self.config_dir)
        self.shared_context = empty_context()
        for index, block in enumerate(assignments):
            self.shared_context = self._evaluate(block, self.shared_context, path_scope=("shared", index))[2]
        if contains_pending(self.shared_context):
            raise ValueError("Shared variables cannot depend on returned values")
        self.directory_locations = {role: [] for role in DIRECTORY_FEATURES}
        self.metadata_writer = FinalMetadataWriter(delete_first=core["delete_final_json_first"])
        self._exported_nodes = set()
        self._executing = False
        self._progress = None
        self._last_progress_snapshot = None
        self.records = []
        self.initial_context = deepcopy(self.shared_context)
        self.preflight_steps = set()
        self._retained_protected_paths = set()
        self.completed_counts = {}
        self.start_index = 0
        self.context_files = ContextFiles(self)
        self.discovery_sources = set()
        self.discovery_source_steps = {}
        self.discovery_runs = {}
        self.satisfied = {}
        self._validate_satisfaction()
        self._discover_inputs()
        self.initial_records = deepcopy(self.records)
        self.nodes = []
        self._planned = False
        self._protected_keys = None
        self._build(self.start_index)

    def _evaluate(self, settings, context, **kwargs):
        return evaluate_settings(settings, context, path_resolver=self.path_resolver, **kwargs)

    def _resolve(self, value, context, **kwargs):
        return resolve(value, context, path_resolver=self.path_resolver, **kwargs)

    def _validate_satisfaction(self):
        configured = {step["name"]: step for step in steps(self.config)}
        for step in self.steps:
            if not step.get("satisfies"):
                continue
            plugin = load_plugin(step["plugin"])
            plugin.file_features()
            if not plugin.var_records_return or not plugin.var_path_return:
                raise ValueError("core:satisfies requires a plugin declaring var_records_return and var_path_return")
            for name, parameter in step["satisfies"].items():
                if name not in configured:
                    raise ValueError(f"{step['name']}: unknown satisfies step {name}")
                outputs = file_parameter_names(load_plugin(configured[name]["plugin"]).file_features(), "output")
                if parameter not in outputs:
                    raise ValueError(f"{step['name']}: {name}.{parameter} is not a declared output parameter")

    def _load_context_collection(self, step, constants, records, *, establish_scenes=False):
        if not any(step.get(control) for control in LOAD_CONTROLS):
            return constants, records
        return self.context_files.load(step, constants, records, establish_scenes=establish_scenes)

    def _save_node_context(self, node):
        if not any(node.step.get(control) for control in SAVE_CONTROLS):
            return
        if node.record is None:
            contexts = self._record_values(node, after=True)
            records = [{**record, "context": context} for record, context in zip(self.records, contexts)]
        else:
            current = _materialize(node.context, {"const": self.constant_values, "var": self.runtime_values[node.record]})
            records = [{**self.records[node.record], "context": current}]
        constants = records[0]["context"]["const"] if records else _materialize(node.context, {"const": self.constant_values})["const"]
        self.context_files.save(node.step, constants, records)

    @timed_core("discovery")
    def _discover_inputs(self):
        # Scene functions establish their records during planning whenever possible.
        self._log_core_start("discovery")
        for index, step in enumerate(self.steps):
            if not step["run"]:
                self.start_index = index + 1
                continue
            plugin = load_plugin(step["plugin"])
            if plugin.var_records_return or step["plugin"] is None and any(step.get(control) for control in LOAD_CONTROLS):
                constants, self.records = self._load_context_collection(step, self.initial_context["const"], self.records,
                                                                      establish_scenes=bool(plugin.var_records_return))
                initialized = self.records or "var" in self.initial_context or step["name"] in self.context_files.scene_loads
                self.initial_context = {"const": constants, **({"var": {}} if initialized else {})}
            if step["plugin"] is None and any(step.get(control) for control in LOAD_CONTROLS):
                # Explicit context-only loaders establish scenes before graph building.
                for record in self.records:
                    record["context"] = self._evaluate(step["settings"], record["context"], path_scope=(step["name"], record["id"]))[2]
                if self.records:
                    constants = self.records[0]["context"]["const"]
                    if any(record["context"]["const"] != constants for record in self.records):
                        raise ValueError(f"{step['name']}: context assignments have conflicting scene constants")
                    self.initial_context["const"] = deepcopy(constants)
                else:
                    self.initial_context = self._evaluate(step["settings"], self.initial_context, path_scope=(step["name"], None))[2]
                if not step.get("skip_plugin_call", False):
                    self.context_files.save(step, self.initial_context["const"], self.records)
                self.preflight_steps.add(index)
                self.start_index = index + 1
                continue
            if not plugin.var_records_return:
                break
            features = _target_features(plugin, step)
            self._register_directories(plugin)
            loaded = self.context_files.loaded.get(step["name"], {})
            declared = [key.replace(":", ".", 1) for key in step["settings"] if key.startswith("var:")]
            loaded_scenes = step["name"] in self.context_files.scene_loads and all(
                all(self._known(self.context_files.restored(step, record["id"]), name) for name in declared)
                for record in self.records if record["id"] in loaded
            )
            if loaded_scenes and file_parameter_names(features, "output"):
                params = self._evaluate(self._settings(step, plugin), self.initial_context,
                    planning=True, records=[record["context"] for record in self.records], aggregate=True,
                    scene_ids=[record["id"] for record in self.records],
                    restored=self.context_files.restored(step, None))[0]
                outputs = _paths(_selected(params, features["output_reuse_paths"]))
                outputs = path(outputs, base_dir=self.config_dir)
                loaded_scenes = (
                    bool(outputs) and not any(os.path.isdir(filename) for filename in outputs)
                    and not _existing_output_failures(outputs,
                        check_validity=step.get("check_validity", self.controls["check_validity"]),
                        validity_check_grid_size=self.controls["validity_check_grid_size"], log_to_console=False,
                        step=step["name"])
                    and step.get("reuse", self.controls["run_from_existing"])
                )
            if loaded_scenes:
                # This is explicitly requested loading, not an implicit sidecar cache.
                assignments = {key: value for key, value in step["settings"].items() if key.startswith("const:")
                               and not self._known(loaded.get(None, {}), key.replace(":", ".", 1))}
                if assignments:
                    self.initial_context = self._evaluate(assignments, self.initial_context,
                        records=[record["context"] for record in self.records], aggregate=True,
                        scene_ids=[record["id"] for record in self.records], path_scope=(step["name"], "loaded"))[2]
                    for record in self.records:
                        record["context"]["const"] = deepcopy(self.initial_context["const"])
                self._register_satisfied_outputs(step, [record for record in self.records if record["id"] in loaded], loaded=True)
                self.preflight_steps.add(index)
                self.start_index = index + 1
                continue
            if step.get("skip_plugin_call", False):
                break
            if file_parameter_names(features, "output"):
                break  # File-producing scene functions use normal planning below.
            if step.get("scope", "aggregate") != "aggregate":
                raise ValueError("Scene-setting plugins run once with aggregate scope")
            settings = self._settings(step, plugin)
            scene_records = [r["context"] for r in self.records] if "var" in self.initial_context else None
            scene_ids = [r["id"] for r in self.records]
            params, updates, current, _ = self._evaluate(
                settings, self.initial_context, path_scope=(step["name"], None), planning=True,
                aggregate=True, records=scene_records, scene_ids=scene_ids,
            )
            self._normalize_directories(current)
            frozen = {
                name: lookup(current, name)
                for name, value in updates.items()
                if not contains_pending(value)
            }
            params, _, current, _ = self._evaluate(
                settings, self.initial_context, path_scope=(step["name"], None), constants=frozen, planning=True,
                aggregate=True, records=scene_records, scene_ids=scene_ids,
            )
            # Defaults and normalized directory roots are also available to parameters.
            self._normalize_directories(current)
            accepted = (
                set(inspect.signature(plugin.function()).parameters) | plugin.options
                if plugin.target or "function" in vars(plugin)
                else set(self.shared)
            )
            shared = self._resolve(
                {
                    key: value
                    for key, value in self.shared.items()
                    if plugin.aliases.get(key, key) in accepted
                },
                available_context(current),
                aggregate=True, records=scene_records,
            )
            features = _target_features(plugin, step)
            for name in file_parameter_names(features, "input") | file_parameter_names(
                features, "output"
            ):
                if name not in params and name in shared:
                    params[name] = shared[name]
            params = _arguments(features, params, current, settings, self.config_dir)
            with timed_preflight(self, step, weight=max(1, len(self.records))):
                returned = plugin.run(params=params, shared=shared)
            self.discovery_runs[step["name"]] = {"params": {**shared, **params}, "returned": returned}
            _, _, current, _ = self._evaluate(
                settings, self.initial_context, path_scope=(step["name"], None), returned=returned, constants=frozen,
                aggregate=True, records=scene_records, scene_ids=scene_ids,
            )
            self._update_var_records(plugin, step, returned, current)
            self.discovery_sources.update(p for record in self.records for p in record["source_paths"])
            for filename in self.discovery_sources:
                self.discovery_source_steps.setdefault(filename, step["name"])
            imported_ids = {
                str(lookup(item, plugin.var_id_return)) if plugin.var_id_return else str(i)
                for i, item in enumerate(lookup(returned, plugin.var_records_return))
            }
            self.context_files.save(step, self.initial_context["const"],
                                    [record for record in self.records if record["id"] in imported_ids])
            self.preflight_steps.add(index)
            self.start_index = index + 1

    def _log_core_start(self, stage):
        """Announce core work before it starts, using the standard step log format."""
        _log("Start", enabled=self.controls["log_to_console"], step=f"core:{stage}")

    @staticmethod
    def _known(context, name):
        try:
            lookup(context, name)
            return True
        except ValueError:
            return False

    def _register_satisfied_outputs(self, step, records, *, loaded=False):
        if not step.get("satisfies"):
            return
        primary = load_plugin(step["plugin"]).var_path_return
        selectors = [primary]
        if loaded:
            selectors.extend(key[4:] for key, value in step["settings"].items()
                             if key.startswith("var:") and value == "returned:" + primary)
        for record in records:
            fields = record["context"]["var"]
            filename = next((lookup(fields, selector) for selector in selectors if self._known(fields, selector)), None)
            if not isinstance(filename, str) or not filename:
                raise ValueError(f"{step['name']}: var_path_return {primary!r} must supply one file path")
            for target, parameter in step["satisfies"].items():
                key = (target, parameter, record["id"])
                if key in self.satisfied and self.satisfied[key] != filename:
                    raise ValueError(f"Conflicting imported outputs for {target}.{parameter}, scene {record['id']}")
                self.satisfied[key] = filename

    @staticmethod
    def _discovered_constants(step):
        """Scene-setting steps can derive shared values after establishing records."""
        settings, names = {}, set()
        for key, value in step["settings"].items():
            if key.startswith("const:") and any(
                ref == "*" or ref.startswith(("var.", "collect:"))
                or any(matches_reference(ref, name) for name in names)
                for ref in references(value)
            ):
                settings[key] = value
                names.add("const." + key[6:].split(".")[0])
        return settings

    @staticmethod
    def _settings(step, plugin=None):
        plugin = plugin or load_plugin(step["plugin"])
        discovered = Workflow._discovered_constants(step) if plugin.var_records_return else {}
        return {
            key: value
            for key, value in step["settings"].items()
            if not (plugin.var_records_return and key.startswith("var:")) and key not in discovered
        }

    def _register_directories(self, plugin):
        for role, name in DIRECTORY_FEATURES.items():
            for selector in getattr(plugin, name):
                if selector not in self.directory_locations[role]:
                    self.directory_locations[role].append(selector)

    def _normalize_directories(self, context, *, required=()):
        return directory_values(
            context, self.directory_locations, base_dir=self.config_dir, required=required
        )

    def _update_var_records(self, plugin, step, returned, context):
        merge = plugin.var_records_mode == "merge"
        if self.controls["protect_source_files"]:
            self._retained_protected_paths.update(
                p for record in self.records for p in record["source_paths"]
            )
        self._register_directories(plugin)
        if plugin.constant_values_return:
            constants = lookup(returned, plugin.constant_values_return)
            if not isinstance(constants, dict):
                raise ValueError("constant_values_return must select a JSON object")
            explicit = {
                key[6:]: lookup(context["const"], key[6:])
                for key in step["settings"]
                if key.startswith("const:") and key not in self._discovered_constants(step)
            }
            if merge:
                context["const"] = _add_missing(context["const"], constants)
            else:
                context["const"].update(deepcopy(constants))
            for name, value in explicit.items():
                assign(context["const"], name, value)
        scenes = lookup(returned, plugin.var_records_return)
        if not isinstance(scenes, list) or any(not isinstance(item, dict) for item in scenes):
            raise ValueError(
                f"{step['plugin']}.{plugin.var_records_return} must return a list of plain dictionaries"
            )
        # Validate actual JSON, without silently stringifying arbitrary Python objects.
        scenes = json.loads(json.dumps(scenes, allow_nan=False))
        self._normalize_directories(context)
        records = deepcopy(self.records) if merge else []
        by_id = {record["id"]: record for record in records}
        for record in records:
            record["context"]["const"] = deepcopy(context["const"])
        seen = set()
        mappings = {key: value for key, value in step["settings"].items() if key.startswith("var:")}
        for index, item in enumerate(scenes):

            def field(selector, default):
                return lookup(item, selector) if selector else default

            scene_id = str(field(plugin.var_id_return, index))
            if scene_id in seen:
                raise ValueError("Scene identifiers must be unique")
            seen.add(scene_id)
            previous = by_id.get(scene_id)
            variables = previous["context"]["var"] if previous else {}
            current = {
                "const": deepcopy(context["const"]),
                "var": _add_missing(variables, item) if merge else deepcopy(item),
            }
            for key, template in mappings.items():
                if previous:
                    try:
                        existing = lookup(variables, key[4:])
                    except ValueError:
                        pass
                    else:
                        if not isinstance(existing, dict):
                            continue
                current = self._evaluate({key: template}, current, path_scope=(step["name"], scene_id), returned=item)[2]
                if previous:
                    # Preserve state before resolving the next dependent assignment.
                    current["var"] = _add_missing(variables, current["var"])
            self._normalize_directories(current)
            sources = field(plugin.source_file_protection_paths_return, [])
            if not isinstance(sources, list) or any(not isinstance(v, str) for v in sources):
                raise ValueError("source_file_protection_paths_return must select a list of paths")
            sources = list(
                dict.fromkeys(
                    [
                        *(previous["source_paths"] if previous else []),
                        *[path(v, base_dir=self.config_dir) for v in sources],
                    ]
                )
            )
            if merge and plugin.source_file_protection_paths_return:
                assign(current["var"], plugin.source_file_protection_paths_return, sources)
            record = {"id": scene_id, "context": current, "source_paths": sources}
            if previous is not None:
                previous.update(record)
            else:
                records.append(record)
        discovered = self._discovered_constants(step)
        if discovered:
            context = self._evaluate(
                discovered, {"const": context["const"], "var": {}}, path_scope=(step["name"], "discovered"),
                records=[r["context"] for r in records], scene_ids=[r["id"] for r in records],
                returned=returned,
            )[2]
            for record in records:
                record["context"]["const"] = deepcopy(context["const"])
                self._normalize_directories(record["context"])
        self.records = records
        if step.get("satisfies"):
            self._register_satisfied_outputs(step, [
                {"id": str(lookup(item, plugin.var_id_return)) if plugin.var_id_return else str(index),
                 "context": {"var": item}}
                for index, item in enumerate(scenes)
            ])
        self.initial_context = {"const": deepcopy(context["const"]), "var": {}}

    @timed_core("build")
    def _build(self, start_index=0):
        self._log_core_start("build")
        self.barrier_index = None
        self.contexts = [deepcopy(r["context"]) for r in self.records]
        self.runtime_values = [{} for _ in self.records]
        self.constant_values = {}
        self.constant_steps = {}
        self.var_constant_writes = {}
        self._restored_nodes = set()
        self.context = deepcopy(self.initial_context)
        producers = [{} for _ in self.records]
        constant_producers, paths, path_producers = {}, {}, {}
        self.step_names = [step["name"] for step in self.steps]
        for step_index in range(start_index, len(self.steps)):
            step = self.steps[step_index]
            if not step["run"]:
                continue
            if any(step.get(control) for control in LOAD_CONTROLS):
                constants, loaded_records = self._load_context_collection(
                    step, self.context["const"],
                    [{**record, "context": context} for record, context in zip(self.records, self.contexts)],
                )
                if [record["id"] for record in loaded_records] != [record["id"] for record in self.records]:
                    raise ValueError("Context loads that establish scenes must precede processing steps")
                self.context["const"] = constants
                self.contexts = [record["context"] for record in loaded_records]
                record_indices = {record["id"]: i for i, record in enumerate(self.records)}
                for scene_id, selected in self.context_files.loaded.get(step["name"], {}).items():
                    for name in dependency_values(selected):
                        if name.startswith("const."):
                            if self._known(self.context, name):
                                constant_producers.pop(name, None)
                        elif scene_id in record_indices:
                            i = record_indices[scene_id]
                            if self._known(self.contexts[i], name):
                                producers[i].pop(name, None)
            plugin = load_plugin(step["plugin"])
            features = _target_features(plugin, step)
            self._register_directories(plugin)
            input_names = file_parameter_names(features, "input")
            output_names = file_parameter_names(features, "output")
            source = plugin.var_records_return is not None
            if source and step.get("scope", "aggregate") != "aggregate":
                raise ValueError("Scene-setting plugins run once with aggregate scope")
            initialized = "var" in self.context
            if step.get("scope") == "var" and not initialized:
                raise ValueError(f"{step['name']}: core:scope var requires initialized var records")
            scope = "aggregate" if source or not initialized else step.get("scope", plugin.scope)
            settings = self._settings(step, plugin)
            constants = {}
            step_producers = deepcopy(constant_producers)
            if scope == "var":
                constants_settings = constant_settings(settings, owner=step["plugin"])
                before = {"const": deepcopy(self.context["const"]), "var": {}}
                _, constants, constant_context, _ = self._evaluate(
                    constants_settings, before, path_scope=(step["name"], "constants"), planning=True, records=self.contexts,
                    scene_ids=[r["id"] for r in self.records],
                )
                if constants_settings:
                    self.constant_steps[step_index] = ConstantBindings(
                        constants_settings, before, constants, records=deepcopy(self.contexts)
                    )
                # Declarative constants inherit their producers; the scene function
                # does not produce them and need not run just to define them.
                for key, template in constants_settings.items():
                    name = key.replace(":", ".", 1)
                    root = ".".join(name.split(".")[:2])
                    dependencies = {
                        i
                        for field, indices in constant_producers.items()
                        for ref in references(template)
                        if matches_reference(ref, field)
                        for i in indices
                    }
                    for ref in references(template):
                        if ref.startswith("collect:"):
                            dependencies.update(
                                index for producer in producers for field, indices in producer.items()
                                if matches_reference("var." + ref[8:], field) for index in indices
                            )
                    if name != root:
                        dependencies.update(constant_producers.get(root, set()))
                    if contains_pending(constants[name]):
                        constant_producers[root] = dependencies
                    elif name == root or not contains_pending(
                        dependency_values(constant_context)[root]
                    ):
                        constant_producers.pop(root, None)
                    step_producers.setdefault(root, set()).update(
                        constant_producers.get(root, set())
                    )
            indices = [None] if scope == "aggregate" else range(len(self.contexts))
            step_records = deepcopy(self.contexts)
            for record_index in indices:
                base = deepcopy(
                    self.context if record_index is None else self.contexts[record_index]
                )
                fields = dependency_values(base)
                available = {
                    k: v
                    for k, v in (
                        step_producers
                        if record_index is None
                        else {**step_producers, **producers[record_index]}
                    ).items()
                    if k.startswith("const.") or contains_pending(fields.get(k))
                }
                scene_id = self.records[record_index]["id"] if record_index is not None else None
                overrides = {}
                for (target, parameter, imported_id), filename in self.satisfied.items():
                    if target == step["name"] and (record_index is None or imported_id == scene_id):
                        if record_index is None:
                            overrides.setdefault(parameter, []).append(filename)
                        else:
                            overrides[parameter] = filename
                if record_index is None:
                    overrides = {name: list(dict.fromkeys(values)) for name, values in overrides.items()}
                    overrides = {name: values[0] if len(values) == 1 else values for name, values in overrides.items()}
                params, updates, current, post = self._evaluate(
                    settings,
                    base, path_scope=(step["name"], record_index),
                    records=step_records,
                    aggregate=record_index is None,
                    scene_ids=[r["id"] for r in self.records],
                    planning=True,
                    constants=constants if scope == "var" else None,
                    restored=self.context_files.restored(step, scene_id),
                    parameter_overrides=overrides,
                )
                if scope == "var":
                    updates = {k: v for k, v in updates.items() if k not in constants}
                for name in input_names | output_names:
                    if name not in params:
                        if name in self.shared:
                            try:
                                params[name] = self._resolve(
                                    self.shared[name], current, records=step_records,
                                    aggregate=record_index is None,
                                    path_key=(name,),
                                )
                            except Deferred:
                                params[name] = Pending(name)
                self._normalize_directories(current)
                base_dir = self.config_dir
                params = _arguments(
                    features, params, current, step["settings"], base_dir,
                    scene_ids=[r["id"] for r in self.records] if record_index is None else None,
                    updates=updates,
                )
                for name in features["output_target_paths"]:
                    value = params.get(name)
                    if contains_pending(value):
                        continue
                    if step.get("require_outputs") is True and value is None:
                        continue
                    paths_to_request = value if isinstance(value, list) else [value]
                    if not paths_to_request or any(
                        not isinstance(p, str) or not p.strip() for p in paths_to_request
                    ):
                        raise ValueError(
                            f"{step['name']} core:require_outputs param:{name} must resolve to "
                            "a nonempty path or flat list of paths"
                        )
                requirements_pending = False
                try:
                    requirements = (
                        _paths(
                            path(
                                self._resolve(
                                    step["requires"],
                                    current,
                                    records=self.contexts if record_index is None else None,
                                ),
                                base_dir=base_dir,
                            )
                        )
                        if step.get("requires")
                        else []
                    )
                except Deferred:
                    requirements, requirements_pending = [], True
                for name, value in _selected(params, output_names).items():
                    filenames = value if isinstance(value, list) else [value]
                    if any(not isinstance(filename, str) or not filename for filename in filenames):
                        raise ValueError(
                            "Output parameters must contain a path or flat list of paths"
                        )
                    checked = name in features["output_collision_check_paths"]
                    for filename in filenames:
                        if filename in paths and (checked or paths[filename]):
                            raise ValueError(f"Output path collision: {filename}")
                        paths[filename] = checked
                deps = set()
                value_dependencies = {}
                directory_refs = {
                    selector
                    for selectors in self.directory_locations.values()
                    for selector in selectors
                    if any(
                        matches_reference(selector, name) and contains_pending(value)
                        for name, value in fields.items()
                    )
                }
                save_refs = set()
                for control in SAVE_CONTROLS:
                    for template, selection in step.get(control, {}).items():
                        save_refs.update(references(template))
                        for name in self.context_files.fields(step, selection, current):
                            # Values assigned here are produced by this invocation;
                            # other selected values are inputs to the save operation.
                            if not any(name == written or name.startswith(written + ".")
                                       for written in updates):
                                save_refs.update(references(name.replace(".", ":", 1)))
                for ref in references([step["settings"], self.shared]) | directory_refs | save_refs:
                    if ref.startswith("collect:") or record_index is None and (ref.startswith("var.") or ref == "*"):
                        key = "var." + ref[8:] if ref.startswith("collect:") else ref
                        for producer in producers:
                            deps.update(
                                i
                                for field, indices in producer.items()
                                if matches_reference(key, field)
                                for i in indices
                            )
                    if not ref.startswith("collect:"):
                        deps.update(
                            i
                            for key, indices in available.items()
                            if matches_reference(ref, key)
                            for i in indices
                        )
                    field_ref = "var." + ref[8:] if ref.startswith("collect:") else ref
                    candidates = available if record_index is not None and not ref.startswith("collect:") else {
                        **step_producers,
                        **{field: set().union(*(producer.get(field, set()) for producer in producers)) for field in {key for producer in producers for key in producer}},
                    }
                    for field, indices in candidates.items():
                        if matches_reference(field_ref, field):
                            for parent in indices:
                                value_dependencies.setdefault(parent, set()).add(field)
                deps.update(
                    path_producers[p]
                    for p in [
                        *_paths(_selected(params, features["input_dependency_paths"])),
                        *requirements,
                    ]
                    if p in path_producers
                )
                dynamic = {name for name, value in updates.items() if contains_pending(value) or any(matches_reference(name, field) for field in post)}
                node = Node(
                    len(self.nodes),
                    step_index,
                    step,
                    record_index,
                    params,
                    requirements,
                    current,
                    base,
                    updates,
                    dynamic,
                    deps,
                    base_dir=base_dir,
                    file_features=features,
                    directory_locations=deepcopy(self.directory_locations),
                    parameter_overrides=overrides,
                    value_dependencies=value_dependencies,
                    requirements_pending=requirements_pending,
                )
                node.collection_snapshot = deepcopy(step_records)
                self.nodes.append(node)
                path_producers.update(
                    {p: node.index for p in node.paths("output_dependency_paths")}
                )
                names = {".".join(name.split(".")[:2]) for name in dynamic}
                static_roots = {
                    ".".join(name.split(".")[:2])
                    for name in updates
                    if not contains_pending(
                        dependency_values(current)[".".join(name.split(".")[:2])]
                    )
                }
                target = constant_producers if record_index is None else producers[record_index]
                for key in static_roots:
                    if record_index is None and key.startswith("var."):
                        for producer in producers:
                            producer.pop(key, None)
                        continue
                    target.pop(key, None)
                target.update({key: {node.index} for key in names if record_index is not None or key.startswith("const.")})
                if record_index is None:
                    update_var_records(
                        self.contexts, updates, [r["id"] for r in self.records]
                    )
                    for producer in producers:
                        producer.update({key: {node.index} for key in names if key.startswith("var.")})
                    self.context = current
                    for scene in self.contexts:
                        scene["const"] = deepcopy(current["const"])
                    node.collection_result = deepcopy(self.contexts)
                else:
                    self.contexts[record_index] = current
            if source:
                self.barrier_index = step_index
                break
            if scope == "var":
                step_nodes = [n for n in self.nodes if n.step_index == step_index]
                for name in {k for n in step_nodes for k in n.updates if k.startswith("const.")}:
                    contributors = [n for n in step_nodes if name in n.updates]
                    known = [n.updates[name] for n in contributors if not contains_pending(n.updates[name])]
                    if known and any(value != known[0] for value in known[1:]):
                        raise ValueError(f"{step['name']} {name} has conflicting scene values; use collect: for a shared list or var: for per-scene values")
                    dynamic = len(known) != len(contributors)
                    value = Pending(name) if dynamic else known[0]
                    assign(constant_context["const"], name[6:], value)
                    root = ".".join(name.split(".")[:2])
                    if dynamic:
                        constant_producers.setdefault(root, set()).update(n.index for n in contributors)
                        for n in contributors:
                            n.dynamic_names.add(name)
                    elif not contains_pending(dependency_values(constant_context)[root]):
                        constant_producers.pop(root, None)
                self.context["const"] = constant_context["const"]
                for scene in self.contexts:
                    scene["const"] = deepcopy(constant_context["const"])
        self.protected_paths = self._retained_protected_paths | {
            p
            for r in self.records
            for p in r["source_paths"]
            if self.controls["protect_source_files"]
        }
        produced = {p for n in self.nodes for p in n.paths(*OUTPUT_PATH_FEATURES)}
        self.protected_paths.update(
            p
            for n in self.nodes
            for p in [*n.paths("input_protection_paths"), *n.requirements]
            if p not in produced
        )

    def _runs(self, node):
        return self.selected_plugin is None or node.step["plugin"] in {None, self.selected_plugin}

    def _failures(self, node, filenames, *, check_validity=None):
        return _existing_output_failures(
            list(filenames),
            check_validity=(
                node.step.get("check_validity", self.controls["check_validity"])
                if check_validity is None
                else check_validity
            ),
            validity_check_grid_size=self.controls["validity_check_grid_size"],
            log_to_console=False,
            step=node.step["name"],
        )

    def _valid(self, node, output_paths):
        filenames = set(output_paths)
        # Reuse always needs existing products. Only selected paths get format checks.
        return (
            # Directory contents are owned by the function (for example tiled outputs).
            # An existing folder alone cannot prove that a batch is complete.
            bool(filenames)
            and not any(os.path.isdir(p) for p in filenames)
            and not self._failures(node, filenames, check_validity=False)
            and not self._failures(node, filenames & set(node.paths("output_validation_paths")))
        )

    def _is_protected(self, filename):
        # Cache source identities once; avoid stat-ing every source for every
        # output in a large batch on a network filesystem.
        def identity(pathname):
            try:
                info = os.stat(pathname)
                return info.st_dev, info.st_ino
            except FileNotFoundError:
                return None

        if self._protected_keys is None:
            paths = {os.path.realpath(p) for p in self.protected_paths}
            identities = {key for p in self.protected_paths if (key := identity(p)) is not None}
            self._protected_keys = paths, identities
        paths, identities = self._protected_keys
        return os.path.realpath(filename) in paths or identity(filename) in identities

    @timed_core("planning")
    def plan(self):
        if self._planned:
            return self
        self._log_core_start("planning")
        for node in self.nodes:
            node.loaded = self._valid(node, node.paths("output_reuse_paths"))

        def require(node, requested_paths=None, requested_values=()):
            requested = set(
                node.paths("output_target_paths") or node.paths("output_reuse_paths")
                if requested_paths is None
                else requested_paths
            )
            if node.needed and requested <= node.demanded_paths and set(requested_values) <= node.demanded_values:
                return
            node.demanded_paths.update(requested)
            node.demanded_values.update(requested_values)
            if node.status == "processing":
                return
            node.needed = True
            plugin = load_plugin(node.step["plugin"])
            reuse = node.step.get("reuse", self.controls["run_from_existing"])
            values_ready = not plugin.var_records_return and all(self._known(node.context, name) for name in node.demanded_values)
            unresolved = any(contains_pending(node.params.get(name)) for name in file_parameter_names(node.file_features, "output"))
            unresolved |= node.requirements_pending
            unresolved |= self.has_skipped_calls and not node.parameter_overrides and (contains_pending(node.params) or contains_pending(node.requirements))
            ready = not unresolved and (
                values_ready and bool(node.demanded_values) and not node.demanded_paths
                or self._valid(node, node.demanded_paths)
            )
            reusable = node.demanded_paths <= set(node.paths("output_reuse_paths"))
            if ready and values_ready and reuse and reusable:
                node.status = "loaded"
                return
            if not self._runs(node):
                if ready and values_ready:
                    node.status = "loaded"
                    return
                raise ValueError(
                    f"Unselected step {node.step['plugin']} has missing outputs or return values required by a later step"
                )
            node.status = "processing"
            for dependency in node.dependencies:
                parent = self.nodes[dependency]
                used_paths = set(parent.paths("output_dependency_paths")) & set(
                    [*node.paths("input_dependency_paths"), *node.requirements]
                )
                needed_values = node.value_dependencies.get(dependency, ())
                require(parent, used_paths if used_paths or needed_values else None, needed_values)

        for node in self.nodes:
            # Even a plugin-only run needs the scene list before its selected
            # processing step can be planned. Unselected setters may only restore.
            if not self._runs(node) and not load_plugin(node.step["plugin"]).var_records_return:
                continue
            targets = node.paths("output_target_paths")
            if (
                targets
                or node.step.get("require_outputs") is True
                or self.selected_plugin is not None and node.step["plugin"] == self.selected_plugin
                or self.selected_step == node.step["name"]
                or self.preparing and not node.step.get("skip_plugin_call", False)
                or load_plugin(node.step["plugin"]).var_records_return
                or node.step.get("skip_plugin_call", False) and (contains_pending(node.params) or node.requirements_pending)
            ):
                local_returns = node.dynamic_names if self.preparing and not node.step.get("skip_plugin_call", False) else ()
                require(node, targets or None, local_returns)
        self._planned = True
        self._discover_ready_var_records()
        return self

    def counts(self):
        self.plan()
        result = {
            name: {"loaded": 0, "processing": 0, "unused": 0}
            for index, name in enumerate(self.step_names)
            if index not in self.preflight_steps
        }
        result.update(deepcopy(self.completed_counts))
        for node in self.nodes:
            row = result[self.step_names[node.step_index]]
            row["loaded"] += int(node.loaded or node.status == "loaded")
            row["processing"] += int(node.status in {"processing", "completed"})
            row["unused"] += int(not node.needed)
        if self.barrier_index is not None:
            for index in range(self.barrier_index + 1, len(self.steps)):
                if self.steps[index]["run"]:
                    result[self.step_names[index]]["pending"] = 1
        return result

    def _runtime_context(self, node):
        bindings = self.constant_steps.get(node.step_index)
        runtime = {
            "const": bindings.runtime_before["const"] if bindings else self.constant_values,
            "var": self.runtime_values[node.record] if node.record is not None else self._aggregate_runtime_variables(),
        }
        return _materialize(node.pre_context, runtime)

    def _prepare_constants(self, step_index):
        bindings = self.constant_steps.get(step_index)
        if bindings is None:
            return
        bindings.runtime_before = _materialize(bindings.before, {"const": self.constant_values})
        # Reuse values already known during planning (including $random/$now).
        # Only references to pending aggregate results need runtime evaluation.
        frozen = {k: v for k, v in bindings.planned.items() if not contains_pending(v)}
        _, bindings.resolved, _, _ = self._evaluate(
            bindings.settings, bindings.runtime_before, path_scope=(self.steps[step_index]["name"], "constants"), constants=frozen,
            records=[_materialize(c, {"const": self.constant_values, "var": values})
                     for c, values in zip(bindings.records, self.runtime_values)],
            scene_ids=[r["id"] for r in self.records],
        )
        for name, value in bindings.resolved.items():
            if contains_pending(bindings.planned[name]):
                assign(self.constant_values, name.split(".", 1)[1], value)

    def _constants(self, node):
        bindings = self.constant_steps.get(node.step_index)
        return bindings.resolved if bindings else None

    def _record_values(self, node=None, *, after=False):
        contexts = (
            node.collection_result if after else node.collection_snapshot
        ) if node is not None else self.contexts
        return [
            _materialize(c, {"const": self.constant_values, "var": values})
            for c, values in zip(contexts, self.runtime_values)
        ]

    def _aggregate_runtime_variables(self):
        return aggregate_variables([{"var": values} for values in self.runtime_values])

    def _publish_values(self, node, values):
        for name, value in values.items():
            scope, field = name.split(".", 1)
            if scope == "const" and node.record is not None and name in node.updates:
                previous = self.var_constant_writes.setdefault(node.step_index, {})
                if name in previous and previous[name] != value:
                    raise ValueError(f"{node.step['name']} {name} has conflicting scene values; use collect: for a shared list or var: for per-scene values")
                previous[name] = deepcopy(value)
            if scope == "var" and node.record is None:
                update_var_records(
                    [{"var": item} for item in self.runtime_values],
                    {name: value}, [r["id"] for r in self.records],
                )
                continue
            target = self.constant_values if scope == "const" else self.runtime_values[node.record]
            assign(target, field, value)

    def _restore(self, node):
        """Publish only values already supplied by explicit context loads."""
        values = {name: lookup(node.context, name) for name in node.dynamic_names if self._known(node.context, name)}
        self._publish_values(node, values)
        self._restored_nodes.add(node.index)

    def _payload(self, node):
        base = self._runtime_context(node)
        node.runtime_base = deepcopy(base)
        params, _, context, _ = self._evaluate(
            self._settings(node.step),
            base, path_scope=(node.step["name"], node.record),
            records=self._record_values(node),
            aggregate=node.record is None,
            scene_ids=[r["id"] for r in self.records],
            constants=self._constants(node),
            parameter_overrides=node.parameter_overrides,
        )
        # Return assignments take effect after invocation. Keep preceding values
        # available in the context snapshot when this call will replace them.
        context = available_context(_materialize(context, base))
        directory_values(
            context,
            node.directory_locations,
            base_dir=self.config_dir,
        )
        node.runtime_directory_context = deepcopy(context)
        shared = self._resolve(
            self.shared, context, records=self._record_values(node), aggregate=node.record is None
        ) if node.step["plugin"] is not None else {}
        file_names = file_parameter_names(node.file_features, "input") | file_parameter_names(
            node.file_features, "output"
        )
        for name in file_names:
            if name not in params:
                if name in shared:
                    params[name] = shared[name]
        params = _arguments(
            node.file_features, params, context, node.step["settings"], node.base_dir,
            scene_ids=[r["id"] for r in self.records] if node.record is None else None,
        )
        output_names = file_parameter_names(node.file_features, "output")
        known_outputs = node.file_arguments(*OUTPUT_PATH_FEATURES)
        if any(params.get(name) != value for name, value in known_outputs.items()):
            raise ValueError(f"{node.step['name']} output paths changed after planning")
        for name in output_names & params.keys():
            value = params[name]
            if value is None:
                continue
            filenames = value if isinstance(value, list) else [value]
            if any(not isinstance(filename, str) or not filename for filename in filenames):
                raise ValueError(f"{node.step['name']}: output parameter {name} must resolve before execution")
            for other in self.nodes if name not in known_outputs else ():
                if other is node:
                    continue
                collisions = set(filenames) & set(other.paths(*OUTPUT_PATH_FEATURES))
                if collisions and (name in node.file_features["output_collision_check_paths"]
                                   or collisions & set(other.paths("output_collision_check_paths"))):
                    raise ValueError(f"Output path collision: {sorted(collisions)[0]}")
        if node.requirements_pending:
            node.requirements = _paths(path(self._resolve(node.step["requires"], context,
                records=self._record_values(node) if node.record is None else None), base_dir=node.base_dir))
            node.requirements_pending = False
        params = _required_arguments(
            params,
            node.step,
            context,
            node.requirements,
            records=self._record_values(node) if node.record is None else None,
            path_resolver=self.path_resolver,
        )
        node.params = params
        owned = node.paths(*OUTPUT_PATH_FEATURES)
        for filename in owned:
            if self._is_protected(filename):
                raise ValueError(
                    f"{node.step['name']} output aliases a protected input: {filename}"
                )
            for protected in set(node.paths("input_protection_paths")) | set(node.requirements):
                if os.path.realpath(filename) == os.path.realpath(protected) or (
                    os.path.exists(filename)
                    and os.path.exists(protected)
                    and os.path.samefile(filename, protected)
                ):
                    raise ValueError(
                        f"{node.step['name']} output aliases a protected input: {filename}"
                    )
        for filename in [*node.paths("input_existence_check_paths"), *node.requirements]:
            if not os.path.exists(filename):
                raise FileNotFoundError(f"{node.step['name']} input does not exist: {filename}")
        failures = self._failures(
            node, node.paths("output_invalid_removal_paths"), check_validity=True
        )
        remove_output_files(
            [p for p, reason in failures.items() if reason != "missing"],
            input_paths=list(self.protected_paths),
        )
        for filename in node.paths("output_parent_creation_paths"):
            Path(filename).parent.mkdir(parents=True, exist_ok=True)
        return node.step["plugin"], params, shared

    def _finish(self, node, returned):
        _, updates, _, _ = self._evaluate(
            self._settings(node.step),
            node.runtime_base, path_scope=(node.step["name"], node.record),
            records=self._record_values(node),
            aggregate=node.record is None,
            scene_ids=[r["id"] for r in self.records],
            returned=returned,
            constants=self._constants(node),
            parameter_overrides=node.parameter_overrides,
        )
        values = {k: updates[k] for k in node.dynamic_names}
        values = json.loads(json.dumps(values, allow_nan=False))
        self._publish_values(node, values)
        source = load_plugin(node.step["plugin"]).var_records_return
        if source:
            self._var_result = (node, returned)
        overview_paths = node.paths("output_overview_calculation_paths")
        if node.step.get("calculate_overviews", False) and overview_paths:
            from vhrharmonize.io.geospatial import calculate_raster_overviews

            scales = self._resolve(self.shared.get("window_scales"), self._runtime_context(node))
            if not scales:
                raise ValueError("core:calculate_overviews requires shared.param:window_scales")
            for filename in overview_paths:
                if Path(filename).suffix.lower() in {".tif", ".tiff"}:
                    with self._progress.overview(node, filename) if self._progress else nullcontext():
                        calculate_raster_overviews(
                            filename, scales, log_to_console=self.controls["log_to_console"]
                        )
        failures = self._failures(node, node.paths("output_validation_paths"))
        if failures:
            raise RuntimeError(
                f"{node.step['name']} did not produce valid declared outputs: {failures}"
            )
        if not load_plugin(node.step["plugin"]).var_records_return:
            self._save_node_context(node)
        node.status = "completed"
        if self._executing:
            self._save_final_metadata(node)
        if self._progress:
            self._progress.complete(node)

    @timed_core("cleanup")
    def _cleanup(self, *, final=False):
        if self.preparing:
            return  # Local intermediates may be required by the remote suffix.
        # Later scenes can reference any earlier result; retain files until their
        # readers are known. A scene reset must never invalidate its own inputs.
        if self.barrier_index is not None:
            return
        if not (
            self.controls["delete_temp_steps_proactively"]
            or final
            and self.controls["delete_temp_dir"]
        ):
            return
        if final:
            self._log_core_start("cleanup")
        for node in self.nodes:
            if node.status not in {"completed", "loaded"}:
                continue
            if node.status == "loaded" and node.index not in self._restored_nodes:
                continue
            outputs = set(node.paths(*OUTPUT_PATH_FEATURES))
            consumers = [
                n
                for n in self.nodes
                if n.needed
                and (
                    node.index in n.dependencies
                    or outputs.intersection(n.paths("input_protection_paths"))
                )
            ]
            # A temporary terminal result is retained for its caller.
            if not consumers or any(
                n.status != "completed" and (n.status != "loaded" or n.index not in self._restored_nodes)
                for n in consumers
            ):
                continue
            files = [
                p
                for p in node.paths("output_temporary_cleanup_paths")
                if _within(p, node.directories["temp_dir"])
                and p not in node.paths("output_target_paths")
            ]
            for filename in files:
                if not os.path.isfile(filename) or self._is_protected(filename):
                    continue
                for candidate in [
                    filename,
                    filename + ".aux.xml",
                    filename + ".msk",
                    filename + ".ovr",
                ]:
                    if os.path.isfile(candidate) and not self._is_protected(candidate):
                        os.unlink(candidate)
                if final and self.controls["delete_temp_dir"]:
                    parent = Path(filename).parent
                    while str(parent) not in node.directories["temp_dir"] and _within(
                        str(parent), node.directories["temp_dir"]
                    ):
                        try:
                            parent.rmdir()
                        except OSError:
                            break
                        parent = parent.parent

    def _discover_ready_var_records(self):
        if self.barrier_index is None:
            return
        node = next(n for n in self.nodes if n.step_index == self.barrier_index)
        if node.step.get("skip_plugin_call", False):
            return
        # Advancing replaces the current graph. Preserve independent unfinished work too.
        if not node.needed or any(
            previous.needed
            and previous.status == "processing"
            and previous.step["plugin"] is not None
            for previous in self.nodes
            if previous is not node
        ):
            return
        if contains_pending(node.params):
            return
        if any(
            not os.path.exists(p)
            for p in [*node.paths("input_existence_check_paths"), *node.requirements]
        ):
            return
        for step_index in range(node.step_index + 1):
            self._prepare_constants(step_index)
            for previous in self.nodes:
                if (
                    previous.step_index == step_index
                    and previous.needed
                    and previous.status == "loaded"
                ):
                    self._restore(previous)
                elif (
                    previous.step_index == step_index
                    and previous.needed
                    and previous.step["plugin"] is None
                    and previous.status == "processing"
                ):
                    self._execute_preflight(previous)
        if node.status == "processing":
            self._execute_preflight(node)
        self._advance_var_records(node.step_index)
        if not self._executing:
            self.preflight_steps.update(range(node.step_index + 1))
            self.start_index = node.step_index + 1
            self.initial_records = deepcopy(self.records)

    def _execute_preflight(self, node):
        scene = self.records[node.record]["id"] if node.record is not None else "all scenes"
        weight = 1 if node.record is not None else max(1, len(self.records))
        with timed_preflight(self, node.step, scene=scene, weight=weight):
            payload = self._payload(node)
            returned = _execute(payload)
            if load_plugin(node.step["plugin"]).var_records_return:
                self.discovery_runs[node.step["name"]] = {"params": {**payload[2], **payload[1]}, "returned": returned}
            self._finish(node, returned)

    def final_nodes(self):
        if self.barrier_index is not None:
            return []
        needed = [node for node in self.nodes if node.needed]
        if not needed:
            return []
        last_step = max(node.step_index for node in needed)
        return [node for node in needed if node.step_index == last_step]

    def _append_metadata(self, context, *, records=None):
        template = self.controls["output_metadata_path"]
        if template is None:
            return
        value = self._resolve(template, context, records=records)
        if not isinstance(value, str) or not value:
            raise ValueError("core:output_metadata_path must resolve to one nonempty path")
        filename = path(value, base_dir=self.config_dir)
        if self._is_protected(filename):
            raise ValueError(f"Final metadata destination aliases a protected input: {filename}")
        self.metadata_writer.append(filename, context)

    def _save_final_metadata(self, node):
        key = (node.step_index, node.record)
        if key in self._exported_nodes or node not in self.final_nodes():
            return
        runtime = {
            "const": self.constant_values,
            "var": self.runtime_values[node.record] if node.record is not None else self._aggregate_runtime_variables(),
        }
        context = available_context(_materialize(node.context, runtime))
        if node.record is None and self.records:
            # Aggregate outputs share the final constants, with one context per scene.
            for scene in self._record_values(node, after=True):
                scene["const"] = deepcopy(context["const"])
                self._append_metadata(available_context(scene), records=self._record_values(node, after=True))
        else:
            self._append_metadata(context)
        self._exported_nodes.add(key)

    def _advance_var_records(self, step_index):
        if step_index != self.barrier_index:
            return
        node, returned = self._var_result
        context = _materialize(node.context, {"const": self.constant_values, "var": {}})
        plugin = load_plugin(node.step["plugin"])
        if plugin.var_records_mode == "merge":
            # Carry actual completed/cached values across the planning boundary.
            for record, current in zip(self.records, self._record_values(node)):
                record["context"] = available_context(current)
        self._update_var_records(plugin, node.step, returned, context)
        self.discovery_sources.update(p for record in self.records for p in record["source_paths"])
        for filename in self.discovery_sources:
            self.discovery_source_steps.setdefault(filename, node.step["name"])
        self.context_files.save(node.step, self.initial_context["const"], self.records)
        counts = self.counts()
        self.completed_counts.update(
            {
                self.step_names[n.step_index]: counts[self.step_names[n.step_index]]
                for n in self.nodes
            }
        )
        self._retained_protected_paths.update(self.protected_paths)
        self.nodes = []
        self._planned = False
        self._protected_keys = None
        self._build(step_index + 1)
        self.plan()
        if self._progress:
            self._progress.sync()

    def _log_node_start(self, node, index, processing_total, total):
        """Report dispatch order using counts owned by the workflow scheduler."""
        if self._progress and self.controls["show_progress"]:
            return
        estimate = None
        if self._progress:
            with self._progress.state.lock:
                estimate = self._progress.state.mean_seconds(self._progress.state.rows[node.step["name"]])
            if estimate is not None:
                estimate *= self._progress.weight(node)
        _log(
            f"Start {index}/{processing_total}/{total}" + (f" | ETA ~{estimate:.0f}s" if estimate is not None else ""),
            enabled=self.controls["log_to_console"],
            step=f"core:{node.step['name']}",
            scene_basename=self.records[node.record]["id"] if node.record is not None else None,
        )

    def _run_horizontal_step(self, step_index, workers, backend):
        self._prepare_constants(step_index)
        nodes = [node for node in self.nodes if node.step_index == step_index and node.needed]
        for node in nodes:
            if node.status == "loaded":
                self._restore(node)
                self._save_final_metadata(node)
        pending = [node for node in nodes if node.status == "processing"]
        if not pending:
            self._advance_var_records(step_index)
            return
        payloads = [(node, self._payload(node)) for node in pending]
        total = sum(node.step_index == step_index for node in self.nodes)

        def dispatches(*, remote=False):
            for index, (node, payload) in enumerate(payloads, start=1):
                self._log_node_start(node, index, len(pending), total)
                yield node, self._progress.payload(node, payload, remote=remote) if self._progress else payload

        if backend == "dask" and pending[0].record is not None:
            from dask.distributed import Client, as_completed as dask_completed

            if workers != 1:
                raise ValueError("Use concurrent_processing: 1 with Dask")
            address = self.controls.get("dask_scheduler_address")
            scheduler_file = self.controls.get("dask_scheduler_file")
            if not address and not scheduler_file:
                raise ValueError("Dask requires dask_scheduler_address or dask_scheduler_file")
            with (
                Client(address) if address else Client(scheduler_file=scheduler_file)
            ) as client, self._progress.dask_client(client) if self._progress else nullcontext():
                futures = {
                    client.submit(_execute, payload, pure=False): node
                    for node, payload in dispatches(remote=True)
                }
                try:
                    for future in dask_completed(futures):
                        self._finish(futures[future], future.result())
                finally:
                    client.cancel(list(futures))
        elif workers > 1 and len(pending) > 1 and pending[0].record is not None:
            # The progress UI has active threads; fork can inherit locked output streams.
            with ProcessPoolExecutor(
                max_workers=min(workers, len(pending)), mp_context=get_context("spawn")
            ) as executor:
                futures = {
                    executor.submit(_execute, payload): node for node, payload in dispatches()
                }
                try:
                    for future in as_completed(futures):
                        self._finish(futures[future], future.result())
                finally:
                    for future in futures:
                        future.cancel()
        else:
            for node, payload in dispatches():
                self._finish(node, _execute(payload))
        self._cleanup()
        self._advance_var_records(step_index)

    def _can_run_vertical(self, step_index):
        step = self.steps[step_index]
        if step.get("skip_plugin_call", False):
            return False
        if not step["run"] or step.get("processing_direction", self.controls["processing_direction"]) != "vertical":
            return False
        nodes = [n for n in self.nodes if n.step_index == step_index and n.needed]
        if not nodes or any(n.record is None for n in nodes):
            return False
        if any(ref.startswith("collect:") for ref in references([step["settings"], self.shared])):
            return False
        bindings = self.constant_steps.get(step_index)
        if bindings is not None and contains_pending(bindings.planned):
            return False
        # All scene invocations must agree on a returned shared value before
        # another step can use it. Static constants already have planned values.
        return not any(name.startswith("const.") for n in nodes for name in n.dynamic_names)

    def _run_vertical_steps(self, start, end, workers, backend):
        nodes = [n for n in self.nodes if start <= n.step_index < end and n.needed]
        indices = {n.index for n in nodes}
        dependencies, previous = {}, {}
        for node in nodes:
            dependencies[node.index] = node.dependencies & indices
            if node.record in previous:
                dependencies[node.index].add(previous[node.record])
            previous[node.record] = node.index
        prepared, completed = set(), {n.index for n in nodes if n.status == "completed"}
        waiting = {n.index: n for n in nodes if n.index not in completed}
        running = {}
        totals = {
            index: sum(n.status == "processing" for n in nodes if n.step_index == index)
            for index in range(start, end)
        }
        started = {index: 0 for index in totals}
        all_totals = {
            index: sum(n.step_index == index for n in self.nodes)
            for index in range(start, end)
        }

        def finish(node, result):
            self._finish(node, result)
            completed.add(node.index)
            self._cleanup()

        def schedule(submit=None, next_completed=None, limit=1, remote=False):
            while waiting or running:
                ready = sorted(
                    (n for n in waiting.values() if dependencies[n.index] <= completed),
                    key=lambda n: (-n.step_index, n.record),
                )
                progressed = False
                for node in ready:
                    if submit is not None and node.status != "loaded" and len(running) >= limit:
                        continue
                    if node.step_index not in prepared:
                        self._prepare_constants(node.step_index)
                        prepared.add(node.step_index)
                    del waiting[node.index]
                    progressed = True
                    if node.status == "loaded":
                        self._restore(node)
                        self._save_final_metadata(node)
                        completed.add(node.index)
                        self._cleanup()
                    else:
                        payload = self._payload(node)
                        started[node.step_index] += 1
                        self._log_node_start(
                            node, started[node.step_index], totals[node.step_index],
                            all_totals[node.step_index],
                        )
                        if self._progress:
                            payload = self._progress.payload(node, payload, remote=remote)
                        if submit is None:
                            finish(node, _execute(payload))
                        else:
                            running[submit(node, payload)] = node
                    # Reconsider downstream nodes immediately instead of filling
                    # the queue with every scene's upstream work first.
                    break
                if progressed:
                    continue
                if running:
                    future = next_completed(running)
                    finish(running.pop(future), future.result())
                elif waiting:
                    raise RuntimeError("Vertical processing cannot satisfy the remaining scene dependencies")

        if not any(totals.values()):
            schedule()
        elif backend == "dask":
            from dask.distributed import Client, as_completed as dask_completed

            if workers != 1:
                raise ValueError("Use concurrent_processing: 1 with Dask")
            address = self.controls.get("dask_scheduler_address")
            scheduler_file = self.controls.get("dask_scheduler_file")
            if not address and not scheduler_file:
                raise ValueError("Dask requires dask_scheduler_address or dask_scheduler_file")
            with (
                Client(address) if address else Client(scheduler_file=scheduler_file)
            ) as client, self._progress.dask_client(client) if self._progress else nullcontext():
                try:
                    schedule(
                        lambda node, payload: client.submit(
                            _execute, payload, pure=False, priority=node.step_index
                        ),
                        lambda futures: next(iter(dask_completed(futures))),
                        limit=len(previous),
                        remote=True,
                    )
                finally:
                    client.cancel(list(running))
        elif workers > 1 and len(previous) > 1:
            # Start clean workers instead of inheriting the progress UI's locks.
            with ProcessPoolExecutor(
                max_workers=min(workers, len(previous)), mp_context=get_context("spawn")
            ) as executor:
                try:
                    schedule(
                        lambda node, payload: executor.submit(_execute, payload),
                        lambda futures: next(as_completed(futures)),
                        limit=workers,
                    )
                finally:
                    for future in running:
                        future.cancel()
        else:
            schedule()
        # All cached contexts in this segment have been restored before cleanup.
        self._cleanup()

    def get_progress(self):
        """Return a detached progress snapshot during/after a reported run, else None."""
        reporter = self._progress
        if reporter is not None:
            return reporter.snapshot()
        return deepcopy(self._last_progress_snapshot)

    def _emit_timing(self, kind, name, started_ns, duration, *, status="completed", attributes=None):
        values = (kind, name, started_ns, duration, status, attributes or {})
        with self._timing_lock:
            if self._timing_callbacks is None:
                if not self._timing_finished:
                    self._timing_buffer.append(values)
                return
            event = {
                "version": 1, "run_id": self._run_id, "run_started_ns": self._run_started_ns,
                "kind": kind, "name": name, "start_time_ns": started_ns,
                "duration_seconds": max(0.0, duration), "status": status,
                "attributes": {key: value for key, value in {
                    "vhr.backend": self.controls["concurrent_processing_backend"],
                    "vhr.processing_direction": self.controls["processing_direction"],
                    "vhr.concurrent_processing": self.controls["concurrent_processing"],
                    "vhr.job_id": os.environ.get("SLURM_JOB_ID"),
                    "vhr.config": getattr(self, "config_path", None), **(attributes or {}),
                }.items() if value is not None},
            }
            for callback in self._timing_callbacks[:]:
                try:
                    callback(deepcopy(event))
                except Exception as exc:
                    self._timing_callbacks.remove(callback)
                    try:
                        warnings.warn(f"Timing consumer disabled: {exc}", RuntimeWarning, stacklevel=2)
                    except RuntimeWarning:
                        pass  # A caller's warnings-as-errors filter must not abort processing.

    def run(self, *, progress_callback=None, progress_path=None, event_callback=None):
        """Execute with optional snapshot and unthrottled timing-event callbacks.

        core:save_statistics_path appends raw timings; core:load_statistics_path
        seeds runtime estimates. Both default to statistics.jsonl in config_dir.
        Callbacks stay in the parent, independently of console logging or the UI.
        """
        from contextlib import ExitStack

        for name, callback in (("progress_callback", progress_callback), ("event_callback", event_callback)):
            if callback is not None and not callable(callback):
                raise TypeError(f"{name} must be callable")
        if self._timing_finished:
            self._run_id = uuid4().hex
            self._run_started_ns, self._run_started = time_ns(), monotonic()
        self._timing_finished = False
        from vhrharmonize.statistics import load_statistics

        def statistics_path(name):
            filename = self.controls[name]
            return path(self._resolve(filename, self.initial_context), base_dir=self.config_dir) if filename is not None else None

        history_path = statistics_path("load_statistics_path")
        self._historical_timings = load_statistics(history_path) if history_path is not None else {}
        with ExitStack() as stack:
            callbacks = [event_callback] if event_callback is not None else []
            filename = statistics_path("save_statistics_path")
            if filename is not None:
                from vhrharmonize.statistics import StatisticsRecorder

                recorder = stack.enter_context(StatisticsRecorder(filename))
                callbacks.append(recorder)
            self._timing_callbacks = callbacks
            buffered, self._timing_buffer = self._timing_buffer, []
            for kind, name, start, duration, status, attributes in buffered:
                self._emit_timing(kind, name, start, duration, status=status, attributes=attributes)
            outcome, attributes = "completed", {}
            try:
                return self._run_with_progress(progress_callback=progress_callback, progress_path=progress_path)
            except BaseException as exc:
                outcome, attributes = "failed", {"error.type": type(exc).__name__}
                raise
            finally:
                snapshot = self.get_progress()
                if snapshot is not None:
                    for row in snapshot["rows"]:
                        self._emit_timing("step_summary", row["name"], time_ns(), 0,
                                          status=row["status"] if row["status"] != "running" else "incomplete",
                                          attributes={"vhr." + key: row[key] for key in
                                                      ("done", "run", "all", "unused", "reused")})
                self._emit_timing("workflow", "workflow", self._run_started_ns,
                                  monotonic() - self._run_started, status=outcome, attributes=attributes)
                self._timing_callbacks = None
                self._timing_finished = True

    def _run_with_progress(self, *, progress_callback=None, progress_path=None):
        """Execute the workflow, optionally publishing snapshots to a parent callback.

        Callbacks, explicit snapshot paths, saved/loaded timings, report_progress
        or show_progress enable reporting. Only show_progress starts the UI.
        """
        if progress_callback is not None and not callable(progress_callback):
            raise TypeError("progress_callback must be callable")
        self.plan()
        self._last_progress_snapshot = None
        if not (self.controls["show_progress"] or self.controls["report_progress"]
                or progress_callback is not None or progress_path is not None or self._timing_callbacks
                or self._historical_timings):
            return self._run()
        from contextlib import ExitStack
        from .progress import WorkflowProgress

        with ExitStack() as stack:
            callbacks = []
            if self.controls["show_progress"]:
                from .progress_terminal import TerminalProgressDisplay

                display = stack.enter_context(TerminalProgressDisplay())
                callbacks.append(display.update)
            if progress_callback is not None:
                callbacks.append(progress_callback)
            reporter = WorkflowProgress(
                self, callbacks=callbacks,
                snapshot_path=progress_path if progress_path is not None else getattr(self, "progress_path", None),
            )
            self._progress = reporter
            try:
                with reporter:
                    return self._run()
            finally:
                self._last_progress_snapshot = reporter.snapshot()
                self._progress = None

    def _prepare_remaining(self, step_index):
        """Freeze actual local results and plan the suffix without plugin calls."""
        first = next((node for node in self.nodes if node.step_index >= step_index), None)
        if first is not None:
            context = _materialize(first.pre_context, {"const": self.constant_values})
            bindings = self.constant_steps.get(first.step_index)
            if bindings is not None:
                context["const"] = _materialize(bindings.before, {"const": self.constant_values})["const"]
            records = self._record_values(first)
        else:
            context = _materialize(self.context, {"const": self.constant_values})
            records = self._record_values()
        self.preparation_state = {"const": deepcopy(context["const"])}
        if "var" in self.initial_context:
            self.preparation_state["scenes"] = {
                record["id"]: deepcopy(values["var"])
                for record, values in zip(self.records, records)
            }
        if contains_pending(self.preparation_state):
            raise ValueError("Local preparation has unresolved return values at its cutoff")
        self.preparation_index = step_index
        self.initial_context = {"const": deepcopy(context["const"])}
        if "scenes" in self.preparation_state:
            self.initial_context["var"] = {}
        for record, values in zip(self.records, records):
            record["context"] = {"const": deepcopy(context["const"]), "var": deepcopy(values["var"])}
        self.initial_records = deepcopy(self.records)
        self.nodes = []
        self._planned = False
        self._build(step_index)
        self.plan()
        if self._progress:
            self._progress.sync()

    @timed_core("execution")
    def _run(self):
        self.plan()
        self._executing = True
        self._log_core_start("execution")
        if self.controls["log_to_console"]:
            print(f"[workflow] Discovered {len(self.records)} input records")
            snapshot = self._progress.snapshot() if self._progress else None
            estimates = {row["name"]: row["eta_seconds"] for row in snapshot["rows"]} if snapshot else {}
            if snapshot and snapshot["total"]["eta_seconds"] is not None:
                print(f"[workflow] Estimated remaining runtime: ~{snapshot['total']['eta_seconds']:.0f}s")
            disabled_steps = {step["name"] for step in self.steps if not step["run"]}
            for step, counts in self.counts().items():
                if step in disabled_steps:
                    print(f"{step}: run:false")
                    continue
                if counts.get("pending"):
                    print(f"{step}: pending scene discovery")
                    continue
                print(
                    f"{step}: loaded: {counts['loaded']} | processing: {counts['processing']} | unused: {counts['unused']}"
                    + (f" | ETA ~{estimates[step]:.0f}s" if estimates.get(step) is not None else "")
                )
        workers = self.controls["concurrent_processing"]
        workers = (os.cpu_count() or 1) if workers == "num_cpu" else int(workers)
        if workers < 1:
            raise ValueError("concurrent_processing must be positive or num_cpu")
        backend = self.controls["concurrent_processing_backend"]
        if backend not in {"process_pool", "dask"}:
            raise ValueError("concurrent_processing_backend must be process_pool or dask")
        step_index = self.start_index
        while step_index < len(self.steps):
            if self.steps[step_index].get("skip_plugin_call", False):
                if any(step["run"] and not step.get("skip_plugin_call", False)
                       for step in self.steps[step_index:]):
                    raise ValueError("Skipped plugin calls must form a suffix of the workflow")
                self._prepare_remaining(step_index)
                return self.records
            if self._can_run_vertical(step_index):
                end = step_index + 1
                while end < len(self.steps):
                    if self.steps[end]["run"] and not self._can_run_vertical(end):
                        break
                    end += 1
                self._run_vertical_steps(step_index, end, workers, backend)
                step_index = end
            else:
                self._run_horizontal_step(step_index, workers, backend)
                step_index += 1
        if self.preparing:
            self._prepare_remaining(len(self.steps))
            return self.records
        # A discovery-only workflow can still export its imported scenes.
        if not self.nodes and not self._exported_nodes and self.records:
            for record in self.records:
                self._append_metadata(record["context"])
        self._cleanup(final=True)
        for record, context in zip(self.records, self._record_values()):
            record["context"] = available_context(context)
        self.context = available_context(
            _materialize(self.context, {"const": self.constant_values, "var": self._aggregate_runtime_variables()})
        )
        return self.records
