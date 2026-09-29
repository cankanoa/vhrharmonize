"""Public Python entry points for loading, planning and running workflows."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from .config import load_config
from vhrharmonize.progress import ProgressCallback, TimingCallback


def load_workflow(
    config: str | Path | Mapping, *, config_dir: str | None = None, plugin: str | None = None
):
    """Build a workflow from a YAML filename or an in-memory mapping.

    Args:
        config: YAML filename or workflow mapping.
        config_dir: Base for core-managed relative paths; defaults to the YAML directory or cwd.
        plugin: Select steps by plugin name; other processing steps expose existing files only.
    """
    from .engine import Workflow

    if isinstance(config, (str, Path)):
        filename = Path(config).expanduser().resolve()
        data = load_config(filename)
        config_dir = config_dir if config_dir is not None else str(filename.parent)
    elif isinstance(config, Mapping):
        data = dict(config)
    else:
        raise TypeError("config must be a YAML filename or a mapping")
    workflow = Workflow(data, config_dir=config_dir or ".", selected_plugin=plugin)
    if isinstance(config, (str, Path)):
        workflow.config_path = str(filename)
        workflow.progress_path = str(filename) + ".progress.json"
    return workflow


def run_workflow(
    config: str | Path | Mapping, *, config_dir: str | None = None, dry_run: bool = False,
    progress_callback: ProgressCallback | None = None, progress_path: str | Path | None = None,
    event_callback: TimingCallback | None = None,
) -> dict[str, dict[str, int]]:
    """Plan a workflow and execute its required plugins in YAML order.

    Args:
        config: YAML filename or workflow mapping.
        config_dir: Base for core-managed relative paths; defaults to the YAML directory or cwd.
        dry_run: Discover scenes and return counts without running ordinary processing steps.
        progress_callback: Receive detached progress snapshots in the parent process; Python only.
        progress_path: Optional snapshot JSON destination, also enables reporting without Rich.
        event_callback: Receive every core timing event in the parent process; Python only.

    Returns:
        Per-step loaded, processing and unused counts from the execution plan.
        Use load_workflow() for access to records, const/var context and individual nodes.
    """
    if not isinstance(dry_run, bool):
        raise TypeError("dry_run must be a boolean")
    workflow = load_workflow(config, config_dir=config_dir)
    counts = workflow.counts()
    if not dry_run:
        workflow.run(progress_callback=progress_callback, progress_path=progress_path, event_callback=event_callback)
        counts = workflow.counts()
    return counts


def run_plugin(
    plugin: str,
    config: str | Path | Mapping,
    *,
    config_dir: str | None = None,
    dry_run: bool = False,
    progress_callback: ProgressCallback | None = None,
    progress_path: str | Path | None = None,
    event_callback: TimingCallback | None = None,
) -> dict[str, dict[str, int]]:
    """Run enabled named steps selecting one plugin in a workflow recipe.

    Args:
        plugin: Registered plugin name. Selected steps still require core:run: true.
        config: YAML filename or workflow mapping.
        config_dir: Base for core-managed relative paths; defaults to the YAML directory or cwd.
        dry_run: Discover scenes and return counts without running ordinary processing steps.
        progress_callback: Receive detached progress snapshots in the parent process; Python only.
        progress_path: Optional snapshot JSON destination, also enables reporting without Rich.
        event_callback: Receive every core timing event in the parent process; Python only.

    Other enabled processing steps supply existing outputs only and are never computed.
    Missing required upstream outputs raise ValueError. Disabled steps are ignored entirely.
    """
    if not isinstance(dry_run, bool):
        raise TypeError("dry_run must be a boolean")
    workflow = load_workflow(config, config_dir=config_dir, plugin=plugin)
    counts = workflow.counts()
    if not dry_run:
        workflow.run(progress_callback=progress_callback, progress_path=progress_path, event_callback=event_callback)
        counts = workflow.counts()
    return counts
