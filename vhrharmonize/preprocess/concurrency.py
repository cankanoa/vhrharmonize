"""Shared concurrency settings for scene and footprint processing."""

from __future__ import annotations

import argparse
import os
import re



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


def _make_dask_client(args: argparse.Namespace):
    """Create a Dask client from generic scheduler connection settings."""
    scheduler_file = getattr(args, "dask_scheduler_file", None)
    scheduler_address = getattr(args, "dask_scheduler_address", None)
    if bool(scheduler_file) == bool(scheduler_address):
        raise ValueError("Dask concurrency requires exactly one of dask_scheduler_file or dask_scheduler_address.")
    try:
        from dask.distributed import Client
    except ImportError as exc:
        raise ImportError(
            "Dask concurrency requires dask.distributed. Install dask[distributed] in the runtime environment."
        ) from exc
    if scheduler_file:
        return Client(scheduler_file=scheduler_file)
    return Client(scheduler_address)
