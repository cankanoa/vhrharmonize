"""Shared concurrency settings for scene and footprint processing."""

from __future__ import annotations

import os
import re
from concurrent.futures import Future, ProcessPoolExecutor
from contextlib import contextmanager
from multiprocessing import get_context


class DaskExecutor:
    """Borrow a workflow client and expose concurrent.futures-compatible jobs."""

    def __init__(self, client):
        self._client = client
        self._jobs = []

    def submit(self, function, /, *args, **kwargs):
        result = Future()
        # A running remote job must finish before its progress transport closes.
        result.set_running_or_notify_cancel()
        job = self._client.submit(function, *args, pure=False, **kwargs)
        self._jobs.append(job)
        def complete(job):
            try:
                value = job.result()
            except BaseException as exc:
                result.set_exception(exc)
            else:
                result.set_result(value)
        job.add_done_callback(complete)
        return result


@contextmanager
def image_executor(workers, backend, scheduler=None):
    """Own an image pool while the aggregate coordinator runs in the parent."""
    if backend == "dask":
        if workers != 1:
            raise ValueError("Use concurrent_processing: 1 with Dask; the cluster controls worker capacity")
        with _make_dask_client(scheduler) as client:
            executor = DaskExecutor(client)
            try:
                yield executor
            finally:
                client.cancel(executor._jobs)
    else:
        with ProcessPoolExecutor(max_workers=workers, mp_context=get_context("spawn")) as executor:
            yield executor


def _resolve_concurrent_processing(value: object) -> int:
    """Resolve the shared concurrency setting.
    Args:
        value: Raw concurrency config value.
    Returns:
        Resolved worker count.
    """
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized == "num_cpu":
            return max(1, os.cpu_count() or 1)
        if not re.fullmatch(r"[-+]?\d+", normalized):
            raise ValueError("concurrent_processing must be an integer or 'num_cpu'.")
        resolved = int(normalized)
    else:
        if not isinstance(value, int):
            raise ValueError("concurrent_processing must be an integer or 'num_cpu'.")
        resolved = int(value)
    if resolved < 1:
        raise ValueError("concurrent_processing must be >= 1.")
    return resolved


def _resolve_concurrent_processing_backend(value: object) -> str:
    """Resolve the shared concurrency backend setting."""
    if value is None:
        return "process_pool"
    if not isinstance(value, str):
        raise ValueError("concurrent_processing_backend must be 'process_pool' or 'dask'.")
    normalized = value.strip().lower().replace("-", "_")
    if normalized not in {"process_pool", "dask"}:
        raise ValueError("concurrent_processing_backend must be 'process_pool' or 'dask'.")
    return normalized


def _validate_dask_scheduler(value):
    """Accept the same file/address connection format as processing plugins."""
    if value is not None and (
        not isinstance(value, (list, tuple))
        or len(value) != 2
        or value[0] not in ("file", "address")
        or not isinstance(value[1], str)
        or not value[1].strip()
    ):
        raise ValueError("dask_scheduler must be ['file', path], ['address', URL], or None.")
    return value


def _make_dask_client(dask_scheduler):
    """Connect to an existing scheduler using the common connection format."""
    _validate_dask_scheduler(dask_scheduler)
    if dask_scheduler is None:
        raise ValueError("Dask requires dask_scheduler: ['file', path] or ['address', URL].")
    try:
        from dask.distributed import Client
    except ImportError as exc:
        raise ImportError(
            "Dask concurrency requires dask.distributed. Install dask[distributed] in the runtime environment."
        ) from exc
    kind, target = dask_scheduler
    if kind == "file":
        return Client(scheduler_file=os.path.abspath(os.path.expanduser(target)))
    return Client(target.strip())
