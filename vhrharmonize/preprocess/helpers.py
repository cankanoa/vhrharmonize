"""Shared lightweight logging helpers."""

from __future__ import annotations


def log(message: str, *, enabled: bool = False, step: str | None = None) -> None:
    """Emit a concise console log message when logging is enabled."""
    if not enabled:
        return
    prefix = f"[{step}] " if step else ""
    print(f"{prefix}{message}")


__all__ = ["log"]
