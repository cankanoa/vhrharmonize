"""Shared lightweight logging helpers."""

from __future__ import annotations


def log(
    message: str,
    *,
    enabled: bool = False,
    step: str | None = None,
    scene_basename: str | None = None,
) -> None:
    """Emit a concise console log message when logging is enabled."""
    if not enabled:
        return
    scene_basename = (scene_basename or "").strip()
    if step and scene_basename:
        prefix = f"[{step}_{scene_basename}] "
    elif step:
        prefix = f"[{step}] "
    else:
        prefix = ""
    print(f"{prefix}{message}", flush=True)


__all__ = ["log"]
