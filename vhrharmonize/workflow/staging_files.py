"""Selective path assignments using declared discovery identifiers."""

import json
import os
import posixpath

from .values import lookup, path, remap_paths, resolve
from .registry import load_plugin


def remote_directory(template, context, selector):
    """Evaluate a destination in its original context, leaving ~ for remote SSH."""
    try:
        value = resolve(template.removeprefix("path:"), context)
    except ValueError as exc:
        raise ValueError(f"HPC path_mappings {selector}: cannot resolve destination: {exc}") from exc
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"HPC path_mappings {selector}: destination must resolve to one nonempty directory string")
    return posixpath.normpath(value)


def selected_values(selector, groups):
    """Read the first available assignment across the complete scene collection."""
    for group in groups:
        found = []
        for context in group:
            try:
                value = lookup(context, selector.replace(":", ".", 1))
            except ValueError:
                continue
            values = value if isinstance(value, list) else [value]
            if any(not isinstance(item, str) or not item.strip() for item in values):
                raise ValueError(f"HPC path_mappings {selector} must contain a path or flat list of paths")
            found.append((context, value))
        if found:
            return found  # An explicitly empty companion list is also resolved.
    raise ValueError(f"Undefined HPC path_mappings reference (or unresolved value): {selector}")


def per_scene_value(pairs, *, selector=None, path_key=False, literal=False):
    """Keep common values small; otherwise select by a scene ID or path suffix.

    Suffix keys also work when the remote home behind ~/ differs from the local
    home. Only path strings/lists are embedded, never scene metadata.
    """
    values = [value for _, value in pairs]
    if all(value == values[0] for value in values):
        value = values[0]
        return value if literal else (["path:" + item for item in value] if isinstance(value, list) else "path:" + value)
    if not selector:
        raise ValueError("Per-scene HPC mappings require a declared scene identifier")
    if not path_key:
        table = {}
        for key, value in pairs:
            if not isinstance(key, (str, int)) or isinstance(key, bool):
                raise ValueError("Per-scene HPC identifiers must be strings or integers")
            key = str(key)
            if key in table and table[key] != value:
                raise ValueError("Conflicting HPC values for the same scene identifier")
            table[key] = value
        expression = f"$lookup({json.dumps(table, ensure_ascii=False)}, $string({selector}))"
        return ("literal:expr:" if literal else "path:expr:") + expression
    if any(not isinstance(filename, str) or not filename for filename, _ in pairs):
        raise ValueError("Per-scene HPC path identifiers must be nonempty paths")
    for depth in range(1, max(len(filename.split("/")) for filename, _ in pairs) + 1):
        table = {}
        for filename, value in pairs:
            key = "/".join(filename.split("/")[-depth:])
            if key in table and table[key] != value:
                break
            table[key] = value
        else:
            break
    else:
        raise ValueError("Conflicting HPC file values for the same imported file")
    key = f"$split({selector}, '/')[-1]" if depth == 1 else (
        f"$join($filter($split({selector}, '/'), function($v, $i, $a) {{"
        f"$i >= $count($a) - {depth}" + "}), '/')")
    return ("literal:expr:" if literal else "path:expr:") + f"$lookup({json.dumps(table, ensure_ascii=False)}, {key})"


def rewrite_file_assignments(staged, workflow, bindings, mapping):
    """Replace only a mapped value's first explicit assignment, when necessary."""
    for selector, values in bindings.items():
        name = next((name for name, block in staged.items() if selector in block and block.get("core:run", False)), None)
        if name is None:
            continue  # Imported fields or explicitly loaded JSON are remapped at their source.
        template = staged[name][selector]
        pairs, unchanged = [], True
        for context, value in values:
            expected = remap_paths(path(value, base_dir=workflow.config_dir), mapping)
            remote = remap_paths(context, mapping)
            pairs.append((context, expected))
            try:
                actual = workflow._resolve(template, remote, returned=remote.get("var", {}))
                unchanged &= path(actual, base_dir=workflow.config_dir) == path(expected, base_dir=workflow.config_dir)
            except (ValueError, TypeError):
                unchanged = False
        if not unchanged:
            staged[name][selector] = context_path_value(pairs, workflow, mapping)


def context_path_value(pairs, workflow, mapping, *, literal=False, plugin=None):
    """Build a path assignment using a discovery plugin's declared scene fields."""
    values = [value for _, value in pairs]
    if all(value == values[0] for value in values):
        return per_scene_value([(None, value) for value in values], literal=literal)
    plugins = [plugin] if plugin is not None else [
        load_plugin(step["plugin"]) for step in workflow.steps if step["run"]
    ]
    for candidate in plugins:
        if not candidate.scene_records_return:
            continue
        for field, path_key in ((candidate.scene_path_return, True), (candidate.scene_id_return, False)):
            if not field:
                continue
            selector = "var." + field
            try:
                keyed = [(remap_paths(lookup(context, selector), mapping), value) for context, value in pairs]
                return per_scene_value(keyed, selector=selector, path_key=path_key, literal=literal)
            except ValueError:
                continue
    raise ValueError("Per-scene HPC mappings need a declared scene_path_return or scene_id_return with unambiguous values")
