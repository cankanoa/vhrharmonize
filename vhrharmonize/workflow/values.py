"""Typed YAML values evaluated against const/var context with JSONata."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import json
import os
import tempfile
from functools import lru_cache

from jsonata import Jsonata
from jsonata.parser import Parser
from jsonata.utils import Utils


UNSET = object()
VAR_UNAVAILABLE = "var is unavailable until a plugin initializes scenes with var_records_return"
AGGREGATE_VAR = "Aggregate parameters have no single scene: use collect: to gather var values"


class _NoScene(dict):
    """JSON-compatible marker for an unavailable aggregate scene scope."""


_NO_SCENE = _NoScene()


class PathResolver:
    """Resolve explicit path values and keep system directories stable per binding."""

    def __init__(self, base_dir="."):
        self.base_dir = os.path.abspath(base_dir)
        self.system_directories = {}

    def __call__(self, value, key=()):
        if isinstance(value, list) and any(not isinstance(v, str) or not v.strip() for v in value):
            raise ValueError("path: must resolve to a path or flat list of paths")
        if value == "sys":
            if key not in self.system_directories:
                self.system_directories[key] = tempfile.mkdtemp(prefix="vhr-")
            return self.system_directories[key]
        return path(value, base_dir=self.base_dir)


def require_var_records(context):
    if "var" not in context:
        raise ValueError(VAR_UNAVAILABLE)


class _BeforeScenes(dict):
    """Reject JSONata lookups of var, including dynamically computed keys."""

    def get(self, key, default=None):
        if key == "var":
            raise ValueError(VAR_UNAVAILABLE)
        return super().get(key, default)


class _AggregateContext(dict):
    def get(self, key, default=None):
        if key == "var":
            raise ValueError(AGGREGATE_VAR)
        return super().get(key, default)


def _contains_no_scene(value):
    if isinstance(value, _NoScene):
        return True
    if isinstance(value, dict):
        return any(_contains_no_scene(v) for v in value.values())
    if isinstance(value, list):
        return any(_contains_no_scene(v) for v in value)
    return False


@dataclass(frozen=True)
class Pending:
    """A variable value that a preceding function has not returned yet."""

    name: str


class Deferred(ValueError):
    pass


def lookup(data, selector):
    if selector in {"", "$"}:
        value = data
    else:
        value = data
        for part in selector.split("."):
            if isinstance(value, Pending):
                raise Deferred(f"Variable is not available until execution: {value.name}")
            try:
                if isinstance(value, list):
                    value = value[int(part)] if part.isdigit() else [lookup(v, part) for v in value]
                else:
                    value = value[part]
            except (KeyError, IndexError, TypeError, ValueError) as exc:
                raise ValueError(f"Undefined variable field: {selector}") from exc
    if contains_pending(value):
        raise Deferred(f"Variable is not available until execution: {selector}")
    return value


def assign(data, name, value):
    """Assign a resolved JSON value, including an update to an existing dotted key."""
    target = data
    parts = name.split(".")
    if not all(part.isidentifier() for part in parts):
        raise ValueError(f"Invalid variable name: {name!r}")
    for part in parts[:-1]:
        target = target.setdefault(part, {})
        if not isinstance(target, dict):
            raise ValueError(f"Cannot assign variable field {name}: parent is not an object")
    target[parts[-1]] = deepcopy(value)


def contains_pending(value):
    if isinstance(value, Pending):
        return True
    if isinstance(value, dict):
        return any(contains_pending(v) for v in value.values())
    if isinstance(value, (list, tuple)):
        return any(contains_pending(v) for v in value)
    return False


def available_context(value):
    """Omit unavailable fields while retaining available siblings in JSON objects."""
    if isinstance(value, dict):
        return {
            key: available_context(item)
            for key, item in value.items()
            if isinstance(item, dict) or not contains_pending(item)
        }
    return deepcopy(value)


def remap_paths(value, mapping):
    """Rebase exact file paths and descendants of declared directories in JSON."""
    if isinstance(value, str):
        if value in mapping:
            return mapping[value]
        for old in sorted(mapping, key=len, reverse=True):
            if value.startswith(old.rstrip("/") + "/"):
                return mapping[old].rstrip("/") + value[len(old.rstrip("/")) :]
        return value
    if isinstance(value, list):
        return [remap_paths(v, mapping) for v in value]
    if isinstance(value, dict):
        return {k: remap_paths(v, mapping) for k, v in value.items()}
    return value


def _compile(source):
    try:
        compiled = Jsonata(source)
        compiled.set_output_convert_nulls(False)
        return compiled
    except Exception as exc:
        raise ValueError(f"Invalid JSONata expression {source!r}: {exc}") from exc


def empty_context():
    return {"const": {}}


def dependency_values(context):
    return {
        f"{scope}.{key}": value
        for scope in ("const", "var")
        for key, value in context.get(scope, {}).items()
    }


def matches_reference(reference, name):
    return (
        reference == "*"
        or reference == name
        or reference.endswith(".*")
        and name.startswith(reference[:-1])
    )


@lru_cache(maxsize=512)
def expression_names(source):
    """Extract namespaced dependencies using JSONata's own parsed expression."""
    names = set()

    def visit(node):
        if isinstance(node, Parser.Symbol):
            if (
                node.type == "function" and node.procedure.type == "variable"
                and node.procedure.value == "lookup" and len(node.arguments) == 2
                and node.arguments[0].type == "variable" and node.arguments[0].value in {"", "$"}
                and node.arguments[1].type == "string" and node.arguments[1].value in {"const", "var"}
            ):
                names.add(node.arguments[1].value + ".*")
                return
            if node.type == "path":
                steps = node.steps
                offset = int(
                    bool(steps) and steps[0].type == "variable" and steps[0].value in {"", "$"}
                )
                first = steps[offset] if len(steps) > offset else None
                if first is not None and first.type == "name" and first.value in {"const", "var"}:
                    second = steps[offset + 1] if len(steps) > offset + 1 else None
                    names.add(
                        first.value
                        + "."
                        + (second.value if second is not None and second.type == "name" else "*")
                    )
                elif offset:
                    names.add("*")
                # Look inside filters, function arguments and nested expressions;
                # don't interpret path segments as standalone root references.
                for step in steps:
                    if step.type == "function":
                        visit(step)
                        continue
                    for key, value in vars(step).items():
                        if key not in {"_outer_instance", "environment"}:
                            visit(value)
                for key, value in vars(node).items():
                    if key not in {"steps", "_outer_instance", "environment"}:
                        visit(value)
                return
            if (
                node.type in {"wildcard", "descendant"}
                or node.type == "variable"
                and node.value in {"", "$"}
            ):
                names.add("*")
            for key, value in vars(node).items():
                if key not in {"_outer_instance", "environment"}:
                    visit(value)
        elif isinstance(node, (list, tuple)):
            for value in node:
                visit(value)
        elif isinstance(node, dict):
            for value in node.values():
                visit(value)

    visit(_compile(source).ast)
    return frozenset(names)


def expression(source, context, *, aggregate=False):
    names = expression_names(source)
    if any(name.startswith("var.") for name in names):
        require_var_records(context)
        if aggregate:
            raise ValueError(AGGREGATE_VAR)
    if any(
        contains_pending(value) and any(matches_reference(ref, key) for ref in names)
        for key, value in dependency_values(context).items()
        if not (aggregate and key.startswith("var."))
    ):
        raise Deferred(f"JSONata expression needs variables from a preceding function: {source}")
    try:
        data = available_context(context)
        if aggregate and "var" in context:
            data = _AggregateContext({"const": data["const"], "var": _NO_SCENE})
        value = _compile(source).evaluate(data if "var" in context else _BeforeScenes(data))
        if _contains_no_scene(value):
            raise ValueError(AGGREGATE_VAR)
    except Exception as exc:
        raise ValueError(f"Cannot evaluate JSONata expression {source!r}: {exc}") from exc
    if value is None:
        raise ValueError(f"JSONata expression returned undefined: {source!r}")
    value = Utils.convert_nulls(value)
    # Whole-context expressions must return plain JSON, not the lookup guard.
    return value if "var" in context else json.loads(json.dumps(value, allow_nan=False))


def resolve(value, context, *, returned=UNSET, records=None, aggregate=False, returned_fields=None,
            path_resolver=None, path_key=()):
    if isinstance(value, list):
        return [resolve(v, context, returned=returned, records=records,
                        aggregate=aggregate, returned_fields=returned_fields,
                        path_resolver=path_resolver, path_key=(*path_key, i)) for i, v in enumerate(value)]
    if isinstance(value, dict):
        return {
            k: resolve(v, context, returned=returned, records=records,
                       aggregate=aggregate, returned_fields=returned_fields,
                       path_resolver=path_resolver, path_key=(*path_key, k)) for k, v in value.items()
        }
    if not isinstance(value, str):
        return value
    kind, colon, text = value.partition(":")
    if not colon:
        return value
    if kind == "literal":
        return text
    if kind == "path":
        resolved = resolve(text, context, returned=returned, records=records,
                           aggregate=aggregate, returned_fields=returned_fields,
                           path_resolver=path_resolver, path_key=path_key)
        return (path_resolver or PathResolver())(resolved, path_key)
    if kind in {"var", "const"}:
        if kind == "var":
            require_var_records(context)
            if aggregate:
                raise ValueError(AGGREGATE_VAR)
        return deepcopy(lookup(context[kind], text))
    if kind == "returned":
        if returned_fields is not None:
            return deepcopy(returned_fields[text])
        if returned is UNSET:
            raise Deferred("The current function has not returned yet")
        return deepcopy(lookup(returned, text))
    if kind == "collect":
        require_var_records(context)
        if records is None:
            raise ValueError("collect: values require initialized scene records")
        return [deepcopy(lookup(record["var"], text)) for record in records]
    if kind == "expr":
        return expression(text, context, aggregate=aggregate)
    return value


def references(value):
    """Variable dependencies, with collect: references kept separate."""
    if isinstance(value, (list, tuple)):
        return set().union(*(references(v) for v in value))
    if isinstance(value, dict):
        return set().union(*(references(v) for v in value.values()))
    if not isinstance(value, str):
        return set()
    kind, _, text = value.partition(":")
    if kind == "path":
        return references(text)
    if kind in {"var", "const"}:
        return {kind + "." + (text.split(".")[0] if text not in {"", "$"} else "*")}
    if kind == "collect":
        return {"collect:" + (text.split(".")[0] if text not in {"", "$"} else "*")}
    if kind == "expr":
        return set(expression_names(text))
    return set()


def uses_returned(value):
    if isinstance(value, (list, tuple)):
        return any(uses_returned(v) for v in value)
    if isinstance(value, dict):
        return any(uses_returned(v) for v in value.values())
    if isinstance(value, str) and value.startswith("path:"):
        return uses_returned(value[5:])
    return isinstance(value, str) and value.startswith("returned:")


def constant_settings(settings, *, owner):
    """Select assignments that can be evaluated once rather than per scene."""
    constants, scene_names = {}, set()
    for key, value in settings.items():
        if not key.startswith("const:"):
            continue
        refs = references(value)
        if uses_returned(value) or any(
            ref.startswith("var.")
            or ref == "*" and refs == {"*"}
            or any(matches_reference(ref, name) for name in scene_names)
            for ref in refs
        ):
            scene_names.add("const." + key[6:].split(".")[0])
        else:
            constants[key] = value
    return constants


def scene_values(value, scene_ids, *, name):
    """Map an aggregate assignment to exactly one value per scene."""
    if isinstance(value, Pending):
        return [deepcopy(value) for _ in scene_ids]
    if isinstance(value, list) and len(value) == len(scene_ids):
        return deepcopy(value)
    if isinstance(value, dict) and set(value) == set(scene_ids):
        return [deepcopy(value[scene_id]) for scene_id in scene_ids]
    raise ValueError(
        f"{name} must map exactly {len(scene_ids)} scenes: supply a list in scene order "
        "or a dictionary with exactly the scene IDs as keys"
    )


def _returned_scene_fields(template, returned, scene_ids, name):
    """Distribute each selected batch return once, including nested assignments."""
    if returned is UNSET:
        return {}
    if isinstance(template, (list, dict)):
        fields = {}
        for item in template.values() if isinstance(template, dict) else template:
            fields.update(_returned_scene_fields(item, returned, scene_ids, name))
        return fields
    if isinstance(template, str) and template.startswith("path:"):
        return _returned_scene_fields(template[5:], returned, scene_ids, name)
    if isinstance(template, str) and template.startswith("returned:"):
        selector = template[9:]
        return {selector: scene_values(lookup(returned, selector), scene_ids, name=name)}
    return {}


def aggregate_variables(records):
    """Internal columns for planning and directory bookkeeping."""
    if not records:
        return {}
    fields = set.intersection(*(set(record["var"]) for record in records))
    return {name: [deepcopy(record["var"][name]) for record in records] for name in fields}


def update_var_records(records, updates, scene_ids):
    """Apply aggregate var assignments without touching the shared constants."""
    for name, value in updates.items():
        if name.startswith("var."):
            for record, item in zip(records, scene_values(value, scene_ids, name=name)):
                assign(record["var"], name[4:], item)


def path(value, *, base_dir):
    """Normalize a resolved path or list of paths; file roles come from plugins."""
    if isinstance(value, list):
        return [path(v, base_dir=base_dir) for v in value]
    if not isinstance(value, str) or not value:
        raise ValueError(f"Expected a non-empty path, got {value!r}")
    return os.path.abspath(os.path.join(base_dir, os.path.expanduser(value)))


def evaluate_settings(
    settings, context, *, records=None, planning=False, returned=UNSET, constants=None,
    scene_ids=None, aggregate=None, path_resolver=None, path_scope=(), restored=None,
    parameter_overrides=None,
):
    """Resolve ordered assignments and arguments without mutating the input context.

    Returned assignments and their dependents are deferred until after invocation.
    No function argument may depend on that same invocation's return value.
    """
    current = deepcopy(context)
    aggregate = records is not None if aggregate is None else aggregate
    if records is not None and "var" not in current:
        records = None
    if records is not None:
        records = deepcopy(records)
        scene_ids = list(map(str, range(len(records)))) if scene_ids is None else scene_ids
        require_var_records(current)
        if aggregate:
            current["var"] = aggregate_variables(records)
    params, updates, post = {}, {}, set()
    parameter_overrides = parameter_overrides or {}
    overridden_variables = {
        settings["param:" + name]: value for name, value in parameter_overrides.items()
        if isinstance(settings.get("param:" + name), str)
        and settings["param:" + name].startswith(("var:", "const:"))
    }
    for key, template in settings.items():
        kind, name = key.split(":", 1)
        if kind == "core":
            continue
        if kind == "var":
            require_var_records(current)
        qualified = kind + "." + name
        if key in overridden_variables:
            assign(current[kind], name, overridden_variables[key])
            updates[qualified] = deepcopy(overridden_variables[key])
            continue
        if kind == "param" and name in parameter_overrides:
            params[name] = deepcopy(parameter_overrides[name])
            if isinstance(template, str) and template.startswith(("var:", "const:")):
                scope, field = template.split(":", 1)
                assign(current[scope], field, params[name])
            continue
        if kind == "const" and constants is not None and qualified in constants:
            # The workflow resolves scene-independent constants once per step.
            # Insert them at their original position to preserve argument ordering.
            value = deepcopy(constants[qualified])
            assign(current[kind], name, value)
            updates[qualified] = value
            continue
        template_refs = references(template)
        refs = {"var." + ref[8:] if ref.startswith("collect:") else ref for ref in template_refs}
        after = uses_returned(template) or any(
            matches_reference(ref, name) for ref in refs for name in post
        )
        if kind == "param" and after:
            raise ValueError(f"param:{name} cannot depend on this function's returned values")

        def evaluate(scope, *, returned_fields=None, batch=False, scene=None):
            try:
                return resolve(template, scope, records=records, returned=returned,
                               aggregate=batch, returned_fields=returned_fields,
                               path_resolver=path_resolver, path_key=(*path_scope, key, scene))
            except ValueError as exc:
                deferred = isinstance(exc, Deferred) or planning and (
                    str(exc).startswith("Undefined variable field:")
                    or str(exc).startswith("JSONata expression returned undefined:")
                )
                if not deferred:
                    raise
                if restored is not None:
                    try:
                        saved = restored.get("scenes", {}).get(scene, restored) if scene is not None else restored
                        return deepcopy(lookup(saved, qualified))
                    except ValueError:
                        pass
                if kind in {"var", "const"} or planning:
                    return Pending(qualified)
                raise

        if aggregate and records is not None and kind == "var":
            fields = _returned_scene_fields(template, returned, scene_ids, key)
            value = [
                evaluate({"const": current["const"], "var": record["var"]},
                         returned_fields={field: items[i] for field, items in fields.items()} if fields else None,
                         scene=scene_ids[i])
                for i, record in enumerate(records)
            ]
        elif aggregate and kind == "const" and (
            any(ref.startswith("var.") for ref in template_refs)
            or "var" in current and template_refs == {"*"}
        ):
            require_var_records(current)
            if not records:
                raise ValueError(f"{key} needs a scene; use collect: for an empty collection")
            values = [evaluate({"const": current["const"], "var": r["var"]}) for r in records]
            if contains_pending(values):
                value = Pending(qualified)
            elif any(item != values[0] for item in values[1:]):
                raise ValueError(f"{key} has conflicting scene values; use collect: for a shared list")
            else:
                value = values[0]
        else:
            value = evaluate(current, batch=bool(aggregate))
        if kind in {"var", "const"}:
            if after:
                post.add(kind + "." + name.split(".")[0])
            if kind == "var" and records is not None and aggregate:
                update_var_records(records, {qualified: value}, scene_ids)
                current["var"] = aggregate_variables(records)
            else:
                assign(current[kind], name, value)
            updates[kind + "." + name] = value
        else:
            params[name] = value
    return params, updates, current, post
