"""Public, JSON-compatible workflow progress contract, independent of any UI."""

from __future__ import annotations

from collections.abc import Callable, Mapping
import json
import math
from pathlib import Path
from typing import TypedDict


PROGRESS_VERSION = 2


class ProgressRow(TypedDict):
    name: str
    unused: int
    done: int
    run: int
    all: int
    reused: int
    percentages: dict[str, float]
    fraction_done: float | None
    active: int
    pending: bool
    worker_progress: bool
    status: str
    eta_seconds: float | None


class ActiveOperation(TypedDict):
    task_id: str
    step: str
    scene: str
    stats: dict
    eta_seconds: float | None


class ProgressSnapshot(TypedDict):
    version: int
    run_id: str
    updated_at: str
    job_id: str | None
    status: str
    total: ProgressRow
    rows: list[ProgressRow]
    active: list[ActiveOperation]
    messages: list[str]


ProgressCallback = Callable[[ProgressSnapshot], None]


class TimingEvent(TypedDict):
    """One completed core measurement, delivered independently of UI refreshes."""

    version: int
    run_id: str
    run_started_ns: int
    kind: str
    name: str
    start_time_ns: int
    duration_seconds: float
    status: str
    attributes: dict


TimingCallback = Callable[[TimingEvent], None]


def validate_progress_snapshot(data: Mapping) -> ProgressSnapshot:
    """Validate and detach a versioned snapshot received from another process.

    Unknown versions and malformed data raise ValueError. Extra fields are
    preserved so compatible additions do not require consumer changes.
    """
    if not isinstance(data, Mapping) or data.get("version") != PROGRESS_VERSION:
        raise ValueError("Unsupported workflow progress snapshot version")
    try:
        result = json.loads(json.dumps(dict(data), allow_nan=False))

        def number(value, *, nullable=False):
            if nullable and value is None:
                return
            if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value < 0:
                raise ValueError("Progress numbers must be finite and non-negative")

        def string(value):
            if not isinstance(value, str):
                raise ValueError("Expected progress text")

        for key in ("run_id", "updated_at", "status"):
            string(result[key])
        if result["job_id"] is not None:
            string(result["job_id"])
        if result["status"] not in {"running", "completed", "failed"}:
            raise ValueError("Invalid workflow progress status")
        for key in ("rows", "active", "messages"):
            if not isinstance(result[key], list):
                raise ValueError(f"Expected progress {key} list")
        for row in [result["total"], *result["rows"]]:
            string(row["name"])
            string(row["status"])
            for key in ("unused", "done", "run", "all", "reused", "active"):
                number(row[key])
                if not isinstance(row[key], int):
                    raise ValueError("Progress counts must be integers")
            for key in ("pending", "worker_progress"):
                if not isinstance(row[key], bool):
                    raise ValueError("Expected progress boolean")
            for key in ("unused", "done", "run", "all", "reused"):
                number(row["percentages"][key])
            number(row["fraction_done"], nullable=True)
            number(row["eta_seconds"], nullable=True)
        for task in result["active"]:
            for key in ("task_id", "step", "scene"):
                string(task[key])
            number(task["eta_seconds"], nullable=True)
            stats = task["stats"]
            for key in ("n", "elapsed"):
                number(stats[key])
            for key in ("total", "rate"):
                number(stats[key], nullable=True)
            for key in ("prefix", "unit"):
                string(stats[key])
        for message in result["messages"]:
            string(message)
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"Malformed workflow progress snapshot: {exc}") from exc
    return result


def read_progress_snapshot(path: str | Path) -> ProgressSnapshot:
    """Read one atomic snapshot file without importing Rich or starting a UI.

    FileNotFoundError means no snapshot has been published at this path yet.
    Invalid JSON, malformed snapshots and unsupported versions raise ValueError.
    """
    with Path(path).open(encoding="utf-8") as stream:
        return validate_progress_snapshot(json.load(stream))


def render_progress(snapshot: ProgressSnapshot, *, console=None):
    """Build a static Rich renderable from a snapshot; never starts a live UI.

    Rich is imported only when this optional frontend is requested.
    """
    from .workflow.progress_rich import RichProgressDisplay

    display = RichProgressDisplay(console=console)
    display.update(snapshot)
    return display.render(live=False)


__all__ = ["PROGRESS_VERSION", "ProgressRow", "ActiveOperation", "ProgressSnapshot",
           "ProgressCallback", "TimingEvent", "TimingCallback", "read_progress_snapshot",
           "validate_progress_snapshot", "render_progress"]
