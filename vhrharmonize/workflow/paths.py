"""Resolve directory roots selected by core or plugins without inventing variables."""

from .values import assign, Deferred, lookup, path

DIRECTORY_FEATURES = {
    "temp_dir": "temporary_directory_context_paths",
    "output_dir": "output_directory_context_paths",
}


def validate_path_mappings(mappings):
    """Map workflow references to directory strings/expressions; sharing is allowed."""
    if not isinstance(mappings, dict):
        raise ValueError("path_mappings must map var:name or const:name to remote directories")
    for selector, remote in mappings.items():
        if not isinstance(selector, str):
            raise ValueError("path_mappings keys must be var:name or const:name")
        scope, colon, name = selector.partition(":")
        if scope not in {"var", "const"} or not colon or not all(p.isidentifier() for p in name.split(".")):
            raise ValueError(f"Invalid path_mappings reference: {selector!r}")
        if not isinstance(remote, str) or not remote.strip():
            raise ValueError(f"path_mappings {selector} needs a remote directory")
    return mappings


def directory_values(context, locations, *, base_dir, required=()):
    """Return normalized roots by role and update only existing declared fields."""
    roots, pending = {role: [] for role in DIRECTORY_FEATURES}, set()
    for role, selectors in locations.items():
        for selector in selectors:
            scope, field = selector.split(".", 1)
            if scope not in context:
                continue
            try:
                value = lookup(context, selector)
            except Deferred:
                pending.add(role)
                continue
            except ValueError:
                continue
            if value is None:
                continue
            resolved = path(value, base_dir=base_dir)
            assign(context[scope], field, resolved)
            for root in resolved if isinstance(resolved, list) else [resolved]:
                if root not in roots[role]:
                    roots[role].append(root)
    for role in required:
        if not roots[role] and role not in pending:
            raise ValueError(
                f"{role} is required by the plugin's file features; an earlier plugin must "
                f"declare {DIRECTORY_FEATURES[role]} and populate that const/var location"
            )
    return roots


def directory_bindings(context, locations):
    """Capture scope-qualified roots for directory bookkeeping and HPC transfers."""
    result = {}
    for selectors in locations.values():
        for selector in selectors:
            try:
                value = lookup(context, selector)
            except ValueError:
                continue
            if isinstance(value, str):
                result[selector] = value
            elif isinstance(value, list):
                result.update({f"{selector}.{index}": root for index, root in enumerate(value)})
    return result
