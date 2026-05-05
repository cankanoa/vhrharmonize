"""Shared workflow path and skip helpers."""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass
from typing import Iterable, List, Optional


@dataclass(frozen=True)
class StepOutputPlan:
    """Resolved output paths and pending work for a workflow step."""

    input_paths: List[str]
    output_paths: List[str]
    pending_input_paths: List[str]
    pending_output_paths: List[str]


def resolve_relative_to_input(path: str, input_folder: str) -> str:
    """Resolve a path against an input folder.
    Args:
        path: Candidate relative or absolute path.
        input_folder: Input folder used for relative resolution.
    Returns:
        Resolved absolute or normalized path.
    """
    if os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(input_folder, path))


def resolve_temp_dir(configured_dir: Optional[str], *, input_folder: str) -> str:
    """Resolve a workflow temp directory.
    Args:
        configured_dir: Configured temp directory or None.
        input_folder: Input folder used for relative resolution.
    Returns:
        Resolved temp directory path.
    """
    if configured_dir in (None, ""):
        return tempfile.mkdtemp(prefix="vhrharmonize_")
    resolved_dir = resolve_relative_to_input(str(configured_dir), input_folder)
    os.makedirs(resolved_dir, exist_ok=True)
    return resolved_dir


def resolve_output_dir(
    configured_dir: Optional[str],
    *,
    temp_dir: str,
    step_name: str,
) -> str:
    """Resolve a workflow output directory.
    Args:
        configured_dir: Configured output directory or None.
        temp_dir: Temp directory used when no explicit output directory is provided.
        step_name: Workflow step name.
    Returns:
        Resolved output directory path.
    """
    output_dir = (
        os.path.join(temp_dir, step_name)
        if configured_dir in (None, "")
        else configured_dir
    )
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def build_output_path_from_input(
    input_path: str,
    output_dir: str,
    *,
    suffix: str = "",
    extension: Optional[str] = None,
) -> str:
    """Build an output path from an input basename.
    Args:
        input_path: Source file path.
        output_dir: Output directory.
        suffix: Suffix to append before the extension.
        extension: Optional extension override.
    Returns:
        Resolved output file path.
    """
    basename = os.path.splitext(os.path.basename(input_path))[0]
    ext = extension if extension is not None else os.path.splitext(input_path)[1]
    return os.path.join(output_dir, f"{basename}{suffix}{ext}")


def plan_step_outputs(
    input_paths: Iterable[str],
    *,
    output_dir: str,
    suffix: str = "",
    extension: Optional[str] = None,
    skip_existing: bool = False,
) -> StepOutputPlan:
    """Plan step outputs and pending work.
    Args:
        input_paths: Input file paths for the step.
        output_dir: Output directory for generated files.
        suffix: Suffix to append before the extension.
        extension: Optional extension override.
        skip_existing: Whether to exclude already existing outputs from pending work.
    Returns:
        Planned output paths and pending work details.
    """
    normalized_input_paths = [str(path) for path in input_paths]
    output_paths = [
        build_output_path_from_input(path, output_dir, suffix=suffix, extension=extension)
        for path in normalized_input_paths
    ]
    if not skip_existing:
        return StepOutputPlan(
            input_paths=normalized_input_paths,
            output_paths=output_paths,
            pending_input_paths=list(normalized_input_paths),
            pending_output_paths=list(output_paths),
        )

    pending_input_paths: List[str] = []
    pending_output_paths: List[str] = []
    for input_path, output_path in zip(normalized_input_paths, output_paths):
        if os.path.exists(output_path):
            continue
        pending_input_paths.append(input_path)
        pending_output_paths.append(output_path)

    return StepOutputPlan(
        input_paths=normalized_input_paths,
        output_paths=output_paths,
        pending_input_paths=pending_input_paths,
        pending_output_paths=pending_output_paths,
    )


__all__ = [
    "StepOutputPlan",
    "build_output_path_from_input",
    "plan_step_outputs",
    "resolve_output_dir",
    "resolve_temp_dir",
    "resolve_relative_to_input",
]
