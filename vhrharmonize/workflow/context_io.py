"""Explicit, selective context snapshots. No implicit cache or plugin restoration."""

from copy import deepcopy
import json
from pathlib import Path

from filelock import FileLock

from vhrharmonize.io.metadata import write_json
from vhrharmonize.io.logging import _log
from .values import assign, lookup, references, path


LOAD_CONTROLS = ("load_context", "load_upsert_context")
SAVE_CONTROLS = ("save_context", "save_upsert_context")
CONTEXT_CONTROLS = (*LOAD_CONTROLS, *SAVE_CONTROLS)


def merge(existing, incoming):
    """Incoming leaves win; objects merge, while lists and scalars replace."""
    if not isinstance(existing, dict) or not isinstance(incoming, dict):
        return deepcopy(incoming)
    result = deepcopy(existing)
    for name, value in incoming.items():
        result[name] = merge(result.get(name), value)
    return result


def validate_operation(value, control):
    if not isinstance(value, dict):
        raise ValueError(f"core:{control} must map file paths to selectors")
    for filename, selection in value.items():
        if not isinstance(filename, str) or not filename.strip():
            raise ValueError(f"core:{control} needs nonempty file paths")
        items = [selection] if isinstance(selection, str) else selection
        if not isinstance(items, list) or not items:
            raise ValueError(f"core:{control} selectors must be a string or nonempty list")
        for item in items:
            if not isinstance(item, str) or not (
                item in {"all", "dependencies", "defined"}
                or item.startswith(("var.", "const."))
                and all(part.isidentifier() for part in item.split("."))
            ):
                raise ValueError(f"Invalid context selector {item!r}; use all, dependencies, defined, var.name or const.name")


def select_fields(selection, settings, definitions, context=None):
    """Expand shorthands using ordinary workflow references, not saved histories."""
    items = [selection] if isinstance(selection, str) else selection
    selected = set(items) - {"all", "dependencies", "defined"}
    if "all" in items:
        selected.update(definitions if context is None else (
            scope + "." + name for scope in ("const", "var") for name in context.get(scope, {})
        ))
    if "defined" in items:
        selected.update(key.replace(":", ".", 1) for key in settings if key.startswith(("var:", "const:")))
    if "dependencies" in items:
        pending = list(references({key: value for key, value in settings.items() if key.startswith("param:")}))
        seen = set()
        while pending:
            ref = pending.pop()
            ref = "var." + ref[8:] if ref.startswith("collect:") else ref
            if ref in seen:
                continue
            seen.add(ref)
            if ref == "*" or ref in {"var.*", "const.*"}:
                available = definitions if context is None else {
                    scope + "." + name for scope in ("const", "var") for name in context.get(scope, {})
                }
                pending.extend(name for name in available if ref == "*" or name.startswith(ref[:-1]))
                continue
            if ref.startswith(("var.", "const.")):
                selected.add(ref)
                for key, definition in definitions.items():
                    if key == ref or key.startswith(ref + ".") or ref.startswith(key + "."):
                        pending.extend(references(definition))
    # Selecting a parent already includes its children.
    return sorted(name for name in selected if not any(name.startswith(parent + ".") for parent in selected))


def validate_snapshot(value, filename="context"):
    if not isinstance(value, dict) or set(value) - {"const", "scenes"}:
        raise ValueError(f"{filename}: context must contain only const and scenes objects")
    if not isinstance(value.get("const", {}), dict) or not isinstance(value.get("scenes", {}), dict):
        raise ValueError(f"{filename}: const and scenes must be objects")
    if any(not isinstance(key, str) or not key or not isinstance(item, dict) for key, item in value.get("scenes", {}).items()):
        raise ValueError(f"{filename}: scenes must map nonempty scene IDs to variable objects")
    json.dumps(value, allow_nan=False)
    return value


class ContextFiles:
    """Parent-owned I/O shared by scene discovery, planning and execution."""

    def __init__(self, workflow):
        self.workflow = workflow
        self.read_paths = set()
        self.read_steps = {}
        self.missing_paths = set()
        self.write_paths = set()
        self.warned = set()
        self.written = set()
        self.loaded = {}
        self.scene_loads = set()
        self.definitions = {}
        definitions = {
            key.replace(":", ".", 1): value
            for block in workflow.config.values() if block.get("plugin") == "shared" and block.get("core:run", False)
            for key, value in block.items() if key.startswith(("var:", "const:"))
        }
        for step in workflow.steps:
            if step["run"]:
                definitions.update({key.replace(":", ".", 1): value for key, value in step["settings"].items() if key.startswith(("var:", "const:"))})
            self.definitions[step["name"]] = deepcopy(definitions)

    def fields(self, step, selection, context=None):
        # Shared parameters count only when accepted by this step's function.
        from .registry import load_plugin
        import inspect

        plugin = load_plugin(step["plugin"])
        accepted = set()
        if plugin.target or "function" in vars(plugin):
            accepted = set(inspect.signature(plugin.function()).parameters) | plugin.options
        settings = {"param:" + key: value for key, value in self.workflow.shared.items() if plugin.aliases.get(key, key) in accepted}
        settings.update(step["settings"])
        return select_fields(selection, settings, self.definitions[step["name"]], context)

    def filename(self, template, context):
        value = self.workflow._resolve(template, context)
        if not isinstance(value, str) or not value.strip():
            raise ValueError("Context file paths must resolve to one nonempty filename")
        return path(value, base_dir=self.workflow.config_dir)

    def read(self, filename, step):
        try:
            value = json.loads(Path(filename).read_text(encoding="utf-8"))
        except FileNotFoundError:
            self.missing_paths.add(filename)
            key = (step["name"], filename)
            if key not in self.warned:
                _log(f"Warning | Context file not found; continuing without loading: {filename}", enabled=self.workflow.controls["log_to_console"], step=step["name"])
                self.warned.add(key)
            return None
        except (OSError, ValueError) as exc:
            raise ValueError(f"Cannot load context {filename}: {exc}") from exc
        self.read_paths.add(filename)
        self.read_steps.setdefault(filename, step["name"])
        return validate_snapshot(value, filename)

    def load(self, step, constants, records, *, establish_scenes=False):
        """Load at a step boundary; return the updated collection and selected values."""
        constants, records = deepcopy(constants), deepcopy(records)
        selected = {}
        for control in LOAD_CONTROLS:
            for template, selection in step.get(control, {}).items():
                sources = records or [{"id": None, "context": {"const": constants}}]
                filenames = list(dict.fromkeys(self.filename(template, {**record["context"], "const": constants}) for record in sources))
                for filename in filenames:
                    snapshot = self.read(filename, step)
                    if snapshot is None:
                        continue
                    available = {"const": {**constants, **snapshot.get("const", {})}, "var": {
                        name: None for variables in [*[record["context"].get("var", {}) for record in records],
                                                     *snapshot.get("scenes", {}).values()] for name in variables
                    }}
                    items = [selection] if isinstance(selection, str) else selection
                    load_all = "all" in items
                    fields = self.fields(step, [item for item in items if item != "all"], available)
                    loads_scenes = any(name.startswith("var.") for name in fields) or load_all and "scenes" in snapshot
                    if loads_scenes:
                        self.scene_loads.add(step["name"])
                    if (not records or establish_scenes) and loads_scenes:
                        existing_ids = {record["id"] for record in records}
                        records.extend({"id": scene_id, "context": {"const": deepcopy(constants), "var": {}}, "source_paths": []}
                                       for scene_id in snapshot.get("scenes", {}) if scene_id not in existing_ids)
                    pairs = [(None, "const", constants, snapshot.get("const", {}))]
                    pairs.extend((record["id"], "var", record["context"]["var"], snapshot["scenes"][record["id"]])
                                 for record in records if record["id"] in snapshot.get("scenes", {}))
                    for scene_id, scope, destination, source in pairs:
                        names = [name for name in fields if name.startswith(scope + ".")]
                        if load_all:
                            # Expand against each saved object, not the destination or
                            # other scenes: their fields need not have the same shape.
                            names = select_fields(["all", *names], {}, {}, {scope: source})
                            selected.setdefault(scene_id, {}).setdefault(scope, {})
                        for name in names:
                            field = name.split(".", 1)[1]
                            try:
                                incoming = lookup(source, field)
                            except ValueError as exc:
                                raise ValueError(f"{step['name']}: {filename} is missing {name}" + (f" for scene {scene_id}" if scene_id is not None else "")) from exc
                            if control == "load_upsert_context":
                                try:
                                    incoming = merge(lookup(destination, field), incoming)
                                except ValueError:
                                    pass
                            assign(destination, field, incoming)
                            assign(selected.setdefault(scene_id, {}), name, incoming)
        for record in records:
            record["context"]["const"] = deepcopy(constants)
        self.loaded[step["name"]] = selected
        return constants, records

    def restored(self, step, scene_id):
        selected = self.loaded.get(step["name"], {})
        if scene_id is None:
            return {**deepcopy(selected.get(None, {})), "scenes": selected}
        return merge(selected.get(None, {}), selected.get(scene_id, {}))

    def save(self, step, constants, records):
        for control in SAVE_CONTROLS:
            grouped = {}
            for template, selection in step.get(control, {}).items():
                for record in records or [{"id": None, "context": {"const": constants}}]:
                    context = {**record["context"], "const": constants}
                    fields = self.fields(step, selection, context)
                    filename = self.filename(template, context)
                    snapshot = grouped.setdefault(filename, {"const": {}, "scenes": {}})
                    if "all" in ([selection] if isinstance(selection, str) else selection) and record["id"] is not None:
                        snapshot["scenes"].setdefault(record["id"], {})
                    for name in fields:
                        if name.startswith("var.") and record["id"] is None:
                            continue  # An explicitly empty scene collection stays empty.
                        try:
                            value = lookup(context, name)
                        except ValueError as exc:
                            raise ValueError(f"{step['name']}: cannot save unresolved/missing {name} to {filename}") from exc
                        scope, field = name.split(".", 1)
                        destination = snapshot["const"] if scope == "const" else snapshot["scenes"].setdefault(record["id"], {})
                        assign(destination, field, value)
            for filename, snapshot in grouped.items():
                validate_snapshot(snapshot, filename)
                Path(filename).parent.mkdir(parents=True, exist_ok=True)
                with FileLock(filename + ".lock"):
                    key = (step["name"], control, filename)
                    # Parent execution can finish scenes individually. A save replaces
                    # the old run once, then retains scenes already written this run.
                    previous = {}
                    if (control == "save_upsert_context" or key in self.written) and Path(filename).exists():
                        previous = validate_snapshot(json.loads(Path(filename).read_text()), filename)
                    if control == "save_upsert_context":
                        value = merge(previous, snapshot)
                    else:
                        value = deepcopy(previous)
                        value.setdefault("const", {}).update(snapshot["const"])
                        value.setdefault("scenes", {}).update(snapshot["scenes"])
                    write_json(filename, value)
                    self.written.add(key)
                    self.write_paths.add(filename)
