"""Workflow progress accounting, worker event transport and snapshot publication."""

from __future__ import annotations

from collections import deque
from contextlib import contextmanager
from dataclasses import dataclass, field
from copy import deepcopy
from datetime import datetime, timezone
import json
import os
import re
from multiprocessing import Manager
from pathlib import Path
from queue import Empty, Queue
from tempfile import NamedTemporaryFile
from threading import RLock, Thread
from time import monotonic
from uuid import uuid4

from vhrharmonize.progress import PROGRESS_VERSION, ProgressCallback, ProgressSnapshot

from vhrharmonize.io.progress import capture_messages, progress_context


@dataclass
class EventReporter:
    queue: object
    task: str

    def event(self, kind, **values):
        self.queue.put({"task": self.task, "kind": kind, **values})

    def __call__(self, **stats):
        fields = ("n", "total", "prefix", "unit", "elapsed", "rate", "operation", "status", "scene")
        self.event("progress", stats={key: stats[key] for key in fields if key in stats})

    def message(self, text):
        self.event("message", text=text)


@dataclass
class _DaskEvents:
    topic: str

    def put(self, event):
        from distributed import get_worker

        get_worker().log_event(self.topic, event)


@dataclass
class StepProgress:
    name: str
    run: int = 0
    all: int = 0
    reused: int = 0
    done: int = 0
    pending: bool = False
    peak_active: int = 0
    unused: int = 0
    worker_progress: bool = True


@dataclass
class TaskProgress:
    step: str
    scene: str
    weight: int = 1
    dispatched: float = field(default_factory=monotonic)
    started: float | None = None
    duration: float | None = None
    done: bool = False
    failed: bool = False
    computed: bool = False
    stats: dict = field(default_factory=dict)
    stats_updated: float | None = None

    @property
    def active(self):
        return self.started is not None and not (self.done or self.failed or self.computed)


class ProgressState:
    """Completion accounting is independent of worker event delivery order."""

    def __init__(self, workers=1):
        self.rows = {}
        self.tasks = {}
        self.messages = deque(maxlen=100)
        self.workers = workers
        self.lock = RLock()
        self.failed = False

    def handle(self, event):
        with self.lock:
            kind = event["kind"]
            if kind == "message":
                self.messages.append(_plain_text(event["text"]))
                return
            task = self.tasks.get(event["task"])
            if task is None:
                return
            if kind == "computed":
                task.duration = max(0.0, event["duration"])
                task.computed = True
            elif kind == "failed":
                task.failed = True
                self.failed = True
            elif not (task.done or task.failed or task.computed):
                if kind == "start":
                    task.started = task.started or monotonic()
                    row = self.rows[task.step]
                    row.peak_active = max(row.peak_active, self.active(task.step))
                elif kind == "progress":
                    task.stats = event["stats"]
                    task.stats_updated = monotonic()

    def active(self, step=None):
        return sum(t.active for t in self.tasks.values() if step is None or t.step == step)

    def complete(self, identity):
        with self.lock:
            task = self.tasks[identity]
            if not task.done:
                task.done = True
                if task.duration is None:
                    task.duration = monotonic() - (task.started or task.dispatched)
                self.rows[task.step].done += task.weight

    def remaining_work(self, row):
        if row.pending:
            return None
        if row.done >= row.run:
            return 0.0
        samples = [t for t in self.tasks.values() if t.step == row.name and t.done]
        if not samples:
            return None
        mean = sum(t.duration for t in samples) / sum(t.weight for t in samples)
        elapsed = sum(
            min(monotonic() - t.started, mean * t.weight)
            for t in self.tasks.values() if t.step == row.name and t.active
        )
        return max(0, mean * (row.run - row.done) - elapsed)

    def total(self):
        rows = list(self.rows.values())
        return StepProgress(
            "total", run=sum(r.run for r in rows), all=sum(r.all for r in rows),
            reused=sum(r.reused for r in rows), done=sum(r.done for r in rows),
            unused=sum(r.unused for r in rows),
            pending=any(r.pending for r in rows),
        )


def _plain_text(value):
    # Keep the transport free of terminal styles, including CSI and OSC escapes.
    return re.sub(r"\x1b(?:\[[0-?]*[ -/]*[@-~]|\][^\x07]*(?:\x07|\x1b\\))", "", str(value))


class WorkflowProgress:
    """Collect worker events and publish detached snapshots to parent consumers."""

    def __init__(self, workflow, *, callbacks=(), snapshot_path=None):
        self.workflow = workflow
        workers = workflow.controls["concurrent_processing"]
        if workers == "num_cpu":
            workers = os.cpu_count() or 1
        self.state = ProgressState(max(1, int(workers)))
        self.queue = Queue()
        self.manager = None
        self.dask_topic = None
        self.snapshot_path = snapshot_path
        self.snapshot_status = "running"
        self.run_id = uuid4().hex
        self.callbacks: list[ProgressCallback] = list(callbacks)
        self._last_snapshot = 0
        self._last_publish = 0
        self.sync()

    def sync(self):
        """Add newly discovered steps without resetting finished step counters."""
        workflow = self.workflow
        with self.state.lock:
            for index, step in enumerate(workflow.steps):
                if not step["run"] or index in workflow.preflight_steps:
                    continue
                name = step["name"]
                existing = self.state.rows.get(name)
                if existing is not None and not existing.pending:
                    continue
                nodes = [n for n in workflow.nodes if n.step_index == index]
                pending = workflow.barrier_index is not None and index > workflow.barrier_index
                weights = [self.weight(node) for node in nodes]
                self.state.rows[name] = StepProgress(
                    name, all=sum(weights),
                    run=sum(w for n, w in zip(nodes, weights) if n.status == "processing"),
                    reused=sum(w for n, w in zip(nodes, weights) if n.loaded and n.status != "processing"),
                    unused=sum(w for n, w in zip(nodes, weights) if not n.needed),
                    worker_progress=self.supports_worker_progress(step),
                    pending=pending,
                )
                if step.get("calculate_overviews") and nodes:
                    counts = [len(self.overview_paths(n)) for n in nodes]
                    if any(counts):
                        label = self.overview_name(step)
                        self.state.rows[label] = StepProgress(
                            label, all=sum(counts),
                            run=sum(c for n, c in zip(nodes, counts) if n.status == "processing"),
                            reused=sum(c for n, c in zip(nodes, counts) if n.loaded and n.status != "processing"),
                            unused=sum(c for n, c in zip(nodes, counts) if not n.needed),
                            worker_progress=False,
                        )
            # Runtime scene discovery can add overview rows after later steps.
            order = [name for step in workflow.steps
                     for name in (step["name"], self.overview_name(step))]
            self.state.rows = {name: self.state.rows[name] for name in order if name in self.state.rows}

    @staticmethod
    def supports_worker_progress(step):
        from inspect import signature
        from vhrharmonize.plugins.base import FunctionPlugin
        from .registry import load_plugin

        if step["plugin"] is None:
            return False  # Context-only steps bind variables without calling a function.
        plugin = load_plugin(step["plugin"])
        if (getattr(plugin.function, "__func__", None) is FunctionPlugin.function
                and not plugin.target):
            return False  # Plugins may implement run() directly without a function target.
        function = plugin.function()
        sources = getattr(function, "__parameter_sources__", ())
        if sources:
            return all(getattr(source(), "__worker_progress__", False) for source in sources)
        declared = getattr(function, "__worker_progress__", None)
        if declared is not None:
            return declared
        return "progress_callback" in signature(function).parameters

    def weight(self, node):
        return max(1, len(self.workflow.records)) if node.record is None else 1

    @staticmethod
    def overview_name(step):
        return f"{step['name']} / overviews"

    @staticmethod
    def overview_paths(node):
        return [p for p in node.paths("output_overview_calculation_paths")
                if Path(p).suffix.lower() in {".tif", ".tiff"}]

    @staticmethod
    def identity(node):
        return f"{node.step_index}:{node.record}"

    def reporter(self, node, *, remote=False, overview=None):
        identity = self.identity(node)
        name = node.step["name"]
        scene = (self.workflow.records[node.record]["id"] if node.record is not None else "all scenes")
        weight = self.weight(node)
        if overview is not None:
            identity += f":overview:{overview}"
            name = self.overview_name(node.step)
            scene, weight = Path(overview).name, 1
        with self.state.lock:
            self.state.tasks[identity] = TaskProgress(name, scene, weight)
        return EventReporter(_DaskEvents(self.dask_topic) if remote else self.queue, identity)

    def payload(self, node, payload, *, remote=False):
        return (*payload, self.reporter(node, remote=remote))

    def complete(self, node):
        identity = self.identity(node)
        if identity in self.state.tasks:
            self.state.complete(identity)

    @contextmanager
    def overview(self, node, filename):
        reporter = self.reporter(node, overview=filename)
        started = monotonic()
        reporter.event("start")
        try:
            with progress_context(reporter, reporter.message):
                yield
        except BaseException:
            reporter.event("failed")
            raise
        else:
            reporter.event("computed", duration=monotonic() - started)
            self.state.complete(reporter.task)

    @contextmanager
    def dask_client(self, client):
        topic = f"vhr-progress-{uuid4().hex}"
        client.subscribe_topic(topic, lambda event: self.queue.put(event[1]))
        self.dask_topic = topic
        self.state.workers = max(1, sum(client.nthreads().values()))
        try:
            yield
        finally:
            self.dask_topic = None
            client.unsubscribe_topic(topic)

    def _consume(self):
        while True:
            try:
                event = self.queue.get(timeout=0.25)
            except Empty:
                self.publish()
                continue
            if event is None:
                return
            self.state.handle(event)
            self.publish()

    def _row_snapshot(self, row, *, total=False):
        failed = self.state.failed if total else any(t.failed for t in self.state.tasks.values() if t.step == row.name)
        eta = None
        if failed:
            status = "failed"
        elif row.pending:
            status = "waiting"
        elif row.done >= row.run:
            status = "completed" if row.run else "reused" if row.reused else "skipped"
            eta = 0.0
        else:
            status = "running"
            if total:
                estimates = [self.state.remaining_work(r) for r in self.state.rows.values()]
                eta = None if any(e is None for e in estimates) else sum(estimates) / self.state.workers
            else:
                work = self.state.remaining_work(row)
                eta = None if work is None else work / max(1, row.peak_active)
        counts = {key: getattr(row, key) for key in ("unused", "done", "run", "all", "reused")}
        return {
            "name": row.name, **counts,
            "percentages": {key: 100 * value / row.all if row.all else 0.0 for key, value in counts.items()},
            "fraction_done": row.done / row.run if row.run else None,
            "active": self.state.active(None if total else row.name),
            "pending": row.pending, "worker_progress": row.worker_progress,
            "status": status, "eta_seconds": eta,
        }

    def snapshot(self) -> ProgressSnapshot:
        """Return detached JSON data; never callbacks, Rich objects or monotonic clocks."""
        with self.state.lock:
            now = monotonic()
            active = []
            for identity, task in self.state.tasks.items():
                if not task.active:
                    continue
                stats = {"n": 0, "total": None, "prefix": "starting", "unit": "it",
                         "elapsed": 0, "rate": None, **deepcopy(task.stats)}
                stats["elapsed"] = ((stats.get("elapsed") or 0) + now - task.stats_updated
                                    if task.stats_updated is not None else now - task.started)
                n, total, rate = stats["n"], stats["total"], stats["rate"]
                eta = max(0.0, (total - n) / rate) if total is not None and rate and rate > 0 else None
                active.append({"task_id": identity, "step": task.step, "scene": str(task.scene),
                               "stats": stats, "eta_seconds": eta})
            return {
                "version": PROGRESS_VERSION, "run_id": self.run_id,
                "updated_at": datetime.now(timezone.utc).isoformat(),
                "job_id": os.environ.get("SLURM_JOB_ID"), "status": self.snapshot_status,
                "total": self._row_snapshot(self.state.total(), total=True),
                "rows": [self._row_snapshot(row) for row in self.state.rows.values()],
                "active": active, "messages": list(self.state.messages)[-5:],
            }

    def publish(self, *, force=False):
        """Publish at most four updates/second, with immediate initial/final updates."""
        if not force and monotonic() - self._last_publish < 0.25:
            return
        self._last_publish = monotonic()
        data = self.snapshot()
        self.save_snapshot(data, force=force)
        for callback in self.callbacks[:]:
            try:
                callback(deepcopy(data))
            except Exception as exc:
                # UI failures must not terminate processing or silence other consumers.
                self.callbacks.remove(callback)
                with self.state.lock:
                    self.state.messages.append(_plain_text(f"Progress consumer disabled: {type(exc).__name__}: {exc}"))

    def save_snapshot(self, data, *, force=False):
        if self.snapshot_path is None or (not force and monotonic() - self._last_snapshot < 1):
            return
        self._last_snapshot = monotonic()
        filename = Path(self.snapshot_path)
        temporary = None
        try:
            with NamedTemporaryFile(mode="w", encoding="utf-8", dir=filename.parent,
                                    prefix=filename.name + ".", suffix=".tmp", delete=False) as stream:
                temporary = Path(stream.name)
                json.dump(data, stream, allow_nan=False)
                stream.write("\n")
            os.replace(temporary, filename)
        except (OSError, ValueError, TypeError) as exc:
            with self.state.lock:
                message = _plain_text(f"Cannot save progress snapshot: {exc}")
                if message not in self.state.messages:
                    self.state.messages.append(message)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    def __enter__(self):
        if self.state.workers > 1 and self.workflow.controls["concurrent_processing_backend"] == "process_pool":
            self.manager = Manager()
            self.queue = self.manager.Queue()
        self.context = progress_context(messages=lambda text: self.queue.put({"kind": "message", "text": text}))
        self.context.__enter__()
        self.capture = capture_messages()
        self.capture.__enter__()
        self.publish(force=True)
        self.reader = Thread(target=self._consume, name="vhr-progress", daemon=True)
        self.reader.start()
        return self

    def __exit__(self, exc_type, exc, traceback):
        try:
            self.capture.__exit__(exc_type, exc, traceback)
            self.context.__exit__(exc_type, exc, traceback)
            self.queue.put(None)
            self.reader.join()
            with self.state.lock:
                if exc is not None:
                    self.state.failed = True
                    self.state.messages.append(_plain_text(f"{type(exc).__name__}: {exc}"))
                    for task in self.state.tasks.values():
                        if not task.done:
                            task.failed = True
                self.snapshot_status = "failed" if self.state.failed else "completed"
            self.publish(force=True)
        finally:
            if self.manager is not None:
                self.manager.shutdown()
