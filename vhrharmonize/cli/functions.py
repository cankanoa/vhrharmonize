"""Generate thin command-line wrappers from Python signatures and docstrings.

Only command syntax and value decoding belong here. Defaults, processing,
validation, and exceptions belong to the invoked Python function.
"""
from __future__ import annotations

import argparse
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, is_dataclass
from datetime import date, datetime
import inspect
import json
from pathlib import Path
from types import UnionType
from typing import Union, get_args, get_origin

from vhrharmonize.parameters import function_parameters, parameter_descriptions


def _decode(value, parameter):
    annotation = parameter.annotation
    variants = get_args(annotation) if get_origin(annotation) in (Union, UnionType) else (annotation,)
    origins = {get_origin(t) or t for t in variants}
    containers = {dict, list, tuple, Iterable, Mapping, Sequence}
    if value.startswith("@") and (str not in variants or origins & containers):
        return json.loads(Path(value[1:]).read_text())
    if annotation is str:
        return value
    if value == "null" and (type(None) in variants or parameter.default is None):
        return None
    if datetime in variants:
        return datetime.fromisoformat(value.replace("Z", "+00:00"))
    if date in variants:
        return date.fromisoformat(value)
    if str in variants and all(t in (str, Path, type(None)) for t in variants):
        return value
    try:
        decoded = json.loads(value)
    except json.JSONDecodeError:
        decoded = value
    # A float-or-string argument may intentionally contain a JSON expression
    # string. Decode containers only when the header accepts that container.
    if str in variants and isinstance(decoded, (dict, list)):
        accepted = {dict, Mapping} if isinstance(decoded, dict) else {list, tuple, Iterable, Sequence}
        if not origins & accepted:
            return value
    if annotation is tuple or getattr(annotation, "__origin__", None) is tuple:
        return tuple(decoded) if isinstance(decoded, list) else decoded
    return decoded


def build_parser(function, *, parser=None, prog=None):
    """Create options directly from a function's names, types, defaults and docs."""
    description = (inspect.getdoc(function) or function.__name__).splitlines()[0]
    parser = parser or argparse.ArgumentParser(prog=prog, description=description, allow_abbrev=False)
    descriptions = parameter_descriptions(function)
    for name, parameter in function_parameters(function).items():
        if name in {"progress_callback", "event_callback"}:
            continue  # Python callable; the workflow supplies this at runtime.
        option = "--" + name.replace("_", "-")
        default = "required" if parameter.default is parameter.empty else f"default: {parameter.default!r}"
        help_text = f"{descriptions.get(name, name.replace('_', ' '))} ({default})"
        if parameter.annotation is bool or isinstance(parameter.default, bool):
            parser.add_argument(option, action=argparse.BooleanOptionalAction, default=argparse.SUPPRESS, help=help_text)
        else:
            parser.add_argument(option, default=argparse.SUPPRESS, help=help_text,
                                type=lambda value, p=parameter: _decode(value, p))
    parser.set_defaults(_function=function)
    return parser


def _print_result(result):
    if result is None:
        return
    if isinstance(result, str):
        print(result)
    else:
        print(json.dumps(result, indent=2, default=lambda value: asdict(value) if is_dataclass(value) else str(value)))


def invoke(parser, argv=None):
    """Decode CLI arguments and call the Python function without intercepting errors."""
    arguments = vars(parser.parse_args(argv))
    function = arguments.pop("_function")
    arguments.pop("_command", None)
    result = function(**arguments)
    _print_result(result)
    return 0


def function_cli(function, argv=None, *, prog=None):
    """Run one Python function as a generated CLI."""
    return invoke(build_parser(function, prog=prog), argv)


def commands_cli(commands, argv=None, *, prog=None, description=None):
    """Generate subcommands for a mapping of names to public Python functions."""
    parser = argparse.ArgumentParser(prog=prog, description=description, allow_abbrev=False)
    subparsers = parser.add_subparsers(dest="_command", required=True)
    for name, function in commands.items():
        description = inspect.getdoc(function) or name
        child = subparsers.add_parser(name, help=description.splitlines()[0], description=description.splitlines()[0], allow_abbrev=False)
        build_parser(function, parser=child)
    return invoke(parser, argv)
