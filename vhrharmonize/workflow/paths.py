"""Resolve directories declared by plugins without inventing context variables."""

from .values import assign, Deferred, lookup, path

DIRECTORY_FEATURES = {
    "temp_dir": "temporary_directory_context_paths",
    "output_dir": "output_directory_context_paths",
}


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
    """Capture scope-qualified roots for checkpoint rebasing and HPC transfers."""
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
