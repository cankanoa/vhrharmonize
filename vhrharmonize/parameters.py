"""Read public function signatures, including explicitly delegated parameters."""

from __future__ import annotations

import inspect
import re
from typing import get_type_hints


def parameter_sources(function):
    """Yield a function and the Python APIs to which its **kwargs are delegated."""
    yield function
    for factory in getattr(function, "__parameter_sources__", ()):
        yield from parameter_sources(factory())


def function_parameters(function):
    """Return named parameters and annotations, with the wrapper taking precedence."""
    parameters = {}
    for source in parameter_sources(function):
        try:
            hints = get_type_hints(source)
        except (NameError, TypeError):
            hints = {}
        for name, parameter in inspect.signature(source).parameters.items():
            if parameter.kind in (parameter.VAR_POSITIONAL, parameter.VAR_KEYWORD):
                continue
            parameters.setdefault(
                name, parameter.replace(annotation=hints.get(name, parameter.annotation))
            )
    return parameters


def parameter_descriptions(function):
    """Read concise argument descriptions from Google-style Args docstrings."""
    descriptions = {}
    for source in parameter_sources(function):
        active = False
        current = None
        for line in (inspect.getdoc(source) or "").splitlines():
            if line.strip() in {"Args:", "Arguments:", "Parameters:"}:
                active = True
                continue
            if not active:
                continue
            if line and not line.startswith(" "):
                break
            match = re.match(r"\s{4}(\w+)(?:\s*\([^)]*\))?:\s*(.*)", line)
            if match:
                current = match.group(1)
                if current in descriptions:
                    current = None
                else:
                    descriptions[current] = match.group(2)
            elif current and line.strip():
                descriptions[current] += " " + line.strip()
    return descriptions
