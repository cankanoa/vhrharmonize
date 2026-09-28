"""Typed YAML values evaluated against const/var context with JSONata."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import json
import os
from functools import lru_cache

from jsonata import Jsonata
from jsonata.parser import Parser
from jsonata.utils import Utils


UNSET = object()
VAR_UNAVAILABLE = "var is unavailable until a plugin initializes scenes with scene_records_return"


def require_scene_variables(context):
    if "var" not in context:
        raise ValueError(VAR_UNAVAILABLE)


class _BeforeScenes(dict):
    """Reject JSONata lookups of var, including dynamically computed keys."""

    def get(self, key, default=None):
        if key == "var":
            raise ValueError(VAR_UNAVAILABLE)
        return super().get(key, default)


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


def expression(source, context):
    names = expression_names(source)
    if any(name.startswith("var.") for name in names):
        require_scene_variables(context)
    if any(
        contains_pending(value) and any(matches_reference(ref, key) for ref in names)
        for key, value in dependency_values(context).items()
    ):
        raise Deferred(f"JSONata expression needs variables from a preceding function: {source}")
    try:
        data = available_context(context)
        value = _compile(source).evaluate(data if "var" in context else _BeforeScenes(data))
    except Exception as exc:
        raise ValueError(f"Cannot evaluate JSONata expression {source!r}: {exc}") from exc
    if value is None:
        raise ValueError(f"JSONata expression returned undefined: {source!r}")
    value = Utils.convert_nulls(value)
    # Whole-context expressions must return plain JSON, not the lookup guard.
    return value if "var" in context else json.loads(json.dumps(value, allow_nan=False))


def resolve(value, context, *, returned=UNSET, records=None):
    if isinstance(value, list):
        return [resolve(v, context, returned=returned, records=records) for v in value]
    if isinstance(value, dict):
        return {
            k: resolve(v, context, returned=returned, records=records) for k, v in value.items()
        }
    if not isinstance(value, str):
        return value
    kind, colon, text = value.partition(":")
    if not colon:
        return value
    if kind == "literal":
        return text
    if kind in {"var", "const"}:
        if kind == "var":
            require_scene_variables(context)
        return deepcopy(lookup(context[kind], text))
    if kind == "returned":
        if returned is UNSET:
            raise Deferred("The current function has not returned yet")
        return deepcopy(lookup(returned, text))
    if kind == "collect":
        require_scene_variables(context)
        if records is None:
            raise ValueError("collect: values require an aggregate plugin")
        return [deepcopy(lookup(record["var"], text)) for record in records]
    if kind == "expr":
        return expression(text, context)
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
    return isinstance(value, str) and value.startswith("returned:")


def constant_settings(settings, *, owner):
    """Select scene-independent assignments, including references nested in JSON."""
    constants = {key: value for key, value in settings.items() if key.startswith("const:")}
    for key, value in constants.items():
        if uses_returned(value) or any(
            ref.startswith(("var.", "collect:")) for ref in references(value)
        ):
            raise ValueError(
                f"{owner} {key} cannot depend on a scene: use literals, const: references, "
                "or JSONata expressions without var, collect:, or returned: values"
            )
    return constants


def path(value, *, base_dir):
    """Normalize a resolved path or list of paths; file roles come from plugins."""
    if isinstance(value, list):
        return [path(v, base_dir=base_dir) for v in value]
    if not isinstance(value, str) or not value:
        raise ValueError(f"Expected a non-empty path, got {value!r}")
    return os.path.abspath(os.path.join(base_dir, os.path.expanduser(value)))


def evaluate_settings(
    settings, context, *, records=None, planning=False, returned=UNSET, constants=None
):
    """Resolve ordered assignments and arguments without mutating the input context.

    Returned assignments and their dependents are deferred until after invocation.
    No function argument may depend on that same invocation's return value.
    """
    current = deepcopy(context)
    params, updates, post = {}, {}, set()
    for key, template in settings.items():
        kind, name = key.split(":", 1)
        if kind == "core":
            continue
        if kind == "var":
            require_scene_variables(current)
        qualified = kind + "." + name
        if kind == "const" and constants is not None and qualified in constants:
            # The workflow resolves scene-independent constants once per step.
            # Insert them at their original position to preserve argument ordering.
            value = deepcopy(constants[qualified])
            assign(current[kind], name, value)
            updates[qualified] = value
            continue
        refs = references(template)
        after = uses_returned(template) or any(
            matches_reference(ref, name) for ref in refs for name in post
        )
        if kind == "param" and after:
            raise ValueError(f"param:{name} cannot depend on this function's returned values")
        try:
            value = resolve(template, current, records=records, returned=returned)
        except Deferred:
            if kind in {"var", "const"}:
                value = Pending(kind + "." + name)
            elif planning:
                value = Pending(kind + "." + name)
            else:
                raise
        if kind in {"var", "const"}:
            if after:
                post.add(kind + "." + name.split(".")[0])
            assign(current[kind], name, value)
            updates[kind + "." + name] = value
        else:
            params[name] = value
    return params, updates, current, post
