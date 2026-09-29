"""Append OpenTelemetry timing spans and summarize them with pandas."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import platform
from tempfile import NamedTemporaryFile
from threading import RLock
from time import time_ns

from .progress import TimingEvent


class _AppendStream:
    """Serialize complete JSON lines across threads and cooperating processes."""

    def __init__(self, path):
        from filelock import FileLock

        self.path = Path(path).expanduser().resolve()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.lock = FileLock(str(self.path) + ".lock", timeout=30)
        self.error = None
        with self.lock, self.path.open("a", encoding="utf-8"):
            pass  # Fail before processing if the configured destination is unwritable.

    def write(self, text):
        try:
            with self.lock, self.path.open("a+", encoding="utf-8") as stream:
                # Preserve a truncated record from an interrupted writer on its own line.
                # The analysis API will identify it by line number, never silently discard it.
                stream.seek(0, os.SEEK_END)
                if stream.tell():
                    with self.path.open("rb") as previous:
                        previous.seek(-1, os.SEEK_END)
                        if previous.read(1) != b"\n":
                            stream.write("\n")
                stream.write(text)
                stream.flush()
        except Exception as exc:
            self.error = exc
            raise

    def flush(self):
        pass  # Every write opens, flushes and closes its file while holding the lock.


class StatisticsRecorder:
    """Consume TimingEvent dictionaries and append the SDK's span JSON as JSONL.

    This uses a private OpenTelemetry provider and never configures an application's
    global provider or contacts a telemetry service. Use as a context manager.
    """

    def __init__(self, path: str | Path):
        from opentelemetry.sdk.resources import Resource
        from opentelemetry.sdk.trace import TracerProvider
        from opentelemetry.sdk.trace.export import ConsoleSpanExporter, SimpleSpanProcessor
        from opentelemetry.sdk.trace.sampling import ALWAYS_ON

        self.stream = _AppendStream(path)
        self.provider = TracerProvider(
            resource=Resource({"service.name": "vhrharmonize", "host.name": platform.node(),
                               "process.pid": os.getpid()}),
            sampler=ALWAYS_ON, shutdown_on_exit=False,
        )
        self.provider.add_span_processor(SimpleSpanProcessor(ConsoleSpanExporter(
            out=self.stream, formatter=lambda span: span.to_json(indent=None) + "\n",
        )))
        self.tracer = self.provider.get_tracer("vhrharmonize.statistics", "1")
        self.roots = {}
        self.lock = RLock()
        self.closed = False

    def __call__(self, event: TimingEvent):
        from opentelemetry.context import Context
        from opentelemetry.trace import Status, StatusCode, set_span_in_context

        if event["version"] != 1:
            raise ValueError("Unsupported timing event version")
        duration = event["duration_seconds"]
        if not math.isfinite(duration) or duration < 0:
            raise ValueError("Timing durations must be finite and non-negative")
        with self.lock:
            if self.closed:
                raise RuntimeError("Statistics recorder is closed")
            run_id = event["run_id"]
            root = self.roots.get(run_id)
            if root is None:
                root = self.tracer.start_span("workflow", context=Context(),
                                              start_time=event["run_started_ns"])
                root.set_attributes({"vhr.schema_version": 1, "vhr.run_id": run_id, "vhr.kind": "workflow"})
                self.roots[run_id] = root
            attributes = {**event["attributes"], "vhr.schema_version": 1, "vhr.run_id": run_id,
                          "vhr.kind": event["kind"], "vhr.name": event["name"],
                          "vhr.status": event["status"], "vhr.duration_seconds": duration}
            span = root if event["kind"] == "workflow" else self.tracer.start_span(
                event["name"], context=set_span_in_context(root), start_time=event["start_time_ns"],
            )
            span.set_attributes(attributes)
            span.set_status(Status(StatusCode.OK if event["status"] in {"completed", "reused", "skipped"}
                                   else StatusCode.ERROR))
            span.end(end_time=event["start_time_ns"] + round(duration * 1_000_000_000))
            if event["kind"] == "workflow":
                del self.roots[run_id]
            if self.stream.error is not None:
                raise OSError(f"Cannot append statistics to {self.stream.path}") from self.stream.error

    def __enter__(self):
        return self

    def __exit__(self, *args):
        from opentelemetry.trace import Status, StatusCode

        with self.lock:
            for root in self.roots.values():
                root.set_attributes({"vhr.status": "incomplete", "vhr.name": "workflow"})
                root.set_status(Status(StatusCode.ERROR, "No workflow completion event received"))
                root.end(end_time=time_ns())
            self.roots.clear()
            self.provider.shutdown()
            self.closed = True


def _raw_lines(source):
    from filelock import FileLock

    with source.open("rb") as stream:
        # Freeze the append boundary, then release the writer lock before parsing
        # a potentially large history. Later runs can keep recording meanwhile.
        with FileLock(str(source) + ".lock", timeout=30):
            remaining = os.fstat(stream.fileno()).st_size
        while remaining:
            line = stream.readline(remaining)
            if not line:
                raise ValueError(f"Statistics input was truncated during analysis: {source}")
            remaining -= len(line)
            yield line


def summarize_statistics(
    input_path: str | Path,
    output_path: str | Path,
    *,
    run_id: str | None = None,
    per_run: bool = False,
) -> dict:
    """Convert raw OpenTelemetry span JSONL into a pandas Table Schema JSON file.

    Args:
        input_path: Append-only raw statistics file created by core or StatisticsRecorder.
        output_path: Separate summary JSON file; replaced atomically on each analysis.
        run_id: Optional workflow run ID to select.
        per_run: Keep runs separate instead of pooling measurements across runs.

    Returns:
        pandas table-oriented JSON data, including its schema. Durations are seconds;
        timings group by kind, step, plugin, backend and outcome. Cached work appears
        in step_summary counts and never contributes zero-time processing samples.
    """
    import pandas as pd

    source = Path(input_path).expanduser().resolve()
    destination = Path(output_path).expanduser().resolve()
    if source == destination or destination.exists() and os.path.samefile(source, destination):
        raise ValueError("Statistics output_path must differ from the raw input_path")
    if not isinstance(per_run, bool):
        raise TypeError("per_run must be a boolean")
    records, seen, finished_runs = [], set(), set()
    for number, line in enumerate(_raw_lines(source), 1):
        if not line.strip():
            continue
        try:
            span = json.loads(line)
            attributes = span["attributes"]
            if attributes["vhr.schema_version"] != 1:
                raise ValueError("unsupported statistics schema")
            identity = (span["context"]["trace_id"], span["context"]["span_id"])
            if identity in seen:
                continue
            seen.add(identity)
            identity_run = attributes["vhr.run_id"]
            if run_id is not None and identity_run != run_id:
                continue
            kind, status = attributes["vhr.kind"], attributes["vhr.status"]
            duration = attributes.get("vhr.duration_seconds")
            if duration is not None and (isinstance(duration, bool) or not isinstance(duration, (int, float))
                                         or not math.isfinite(duration) or duration < 0):
                raise ValueError("invalid timing duration")
            if kind == "workflow" and status in {"completed", "failed"}:
                finished_runs.add(identity_run)
            record = {
                "run_id": identity_run, "kind": kind, "step": attributes["vhr.name"],
                "plugin": attributes.get("vhr.plugin", ""),
                "backend": attributes.get("vhr.backend", ""), "status": status,
                "seconds": duration if kind != "step_summary" and status != "incomplete" else None,
                "scene_units": attributes.get("vhr.scene_units", 0),
                **{key: attributes.get("vhr." + key, 0) for key in ("done", "run", "all", "reused", "unused")},
            }
            records.append(record)
        except (ValueError, TypeError, KeyError, AttributeError) as exc:
            raise ValueError(f"Invalid statistics record at {source}:{number}: {exc}") from exc
    if run_id is not None and not records:
        raise ValueError(f"No statistics for run_id {run_id!r}")
    columns = (["run_id"] if per_run else []) + ["kind", "step", "plugin", "backend", "status"]
    summaries = []
    if records:
        frame = pd.DataFrame.from_records(records)
        for key, group in frame.groupby(columns, sort=True, dropna=False):
            values = group["seconds"].dropna().astype(float)
            quantiles = values.quantile([0.5, 0.9, 0.95])
            summaries.append({
                **dict(zip(columns, key)), "runs": int(group["run_id"].nunique()),
                "unfinished_runs": len(set(group["run_id"]) - finished_runs),
                "records": len(group), "samples": len(values),
                "total_seconds": values.sum() if len(values) else None,
                "mean_seconds": values.mean(), "std_seconds": values.std(),
                "min_seconds": values.min(), "p50_seconds": quantiles.loc[0.5],
                "p90_seconds": quantiles.loc[0.9], "p95_seconds": quantiles.loc[0.95],
                "max_seconds": values.max(),
                **{key: int(group[key].sum()) for key in ("scene_units", "done", "run", "all", "reused", "unused")},
            })
    result_columns = columns + ["runs", "unfinished_runs", "records", "samples", "total_seconds",
                                "mean_seconds", "std_seconds", "min_seconds", "p50_seconds",
                                "p90_seconds", "p95_seconds", "max_seconds", "scene_units",
                                "done", "run", "all", "reused", "unused"]
    result = json.loads(pd.DataFrame(summaries, columns=result_columns).to_json(orient="table", index=False))
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with NamedTemporaryFile(mode="w", encoding="utf-8", dir=destination.parent,
                                prefix=destination.name + ".", suffix=".tmp", delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(result, stream, allow_nan=False, indent=2)
            stream.write("\n")
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return result


__all__ = ["StatisticsRecorder", "summarize_statistics"]
