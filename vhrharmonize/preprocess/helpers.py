"""Shared lightweight logging helpers."""

from __future__ import annotations

import os
from contextlib import contextmanager
from contextvars import ContextVar
from functools import wraps
from inspect import signature


_active_scene = ContextVar("processing_scene", default=None)
_active_step = ContextVar("processing_step", default=None)


def _log_step_start(step: str, *, enabled: bool = False, uppercase: bool = True) -> None:
    """Retain discovery announcements; scene steps announce their own lifecycle."""
    if not uppercase:
        _log("Start", enabled=enabled, step=step)


def _log_image_start(scene_basename: str, inputs, outputs, *, enabled: bool = False, step: str | None = None) -> None:
    parts = ["Start"]
    for label, paths in (("in", inputs), ("out", outputs)):
        if paths:
            parts.append(f"{label}=" + ", ".join(os.path.basename(str(path)) for path in paths))
    _log(" | ".join(parts), enabled=enabled, step=step, scene_basename=scene_basename)


def _log_image_completed(scene_basename: str, index: int, total: int, *, enabled: bool = False, step: str | None = None) -> None:
    _log(f"Completed {index}/{total}", enabled=enabled, step=step, scene_basename=scene_basename)


@contextmanager
def _processing_step(step, scene_basename, inputs, outputs, *, enabled=False, index=1, total=1, announce_step=True, allow_nested=False):
    """Log successful completion only; nested operations share the scene context."""
    if _active_scene.get() is not None and not allow_nested:
        yield
        return
    token = _active_scene.set(scene_basename)
    step_token = _active_step.set(step)
    try:
        _log_image_start(scene_basename, inputs, outputs, enabled=enabled)
        yield
    except BaseException:
        _log("Failed", enabled=enabled, scene_basename=scene_basename)
        raise
    else:
        _log_image_completed(scene_basename, index, total, enabled=enabled)
    finally:
        _active_step.reset(step_token)
        _active_scene.reset(token)


def _logged_operation(step: str, *, inputs=(), outputs=(), allow_nested=False):
    """Give standalone preprocessing calls the same lifecycle as workflow steps."""
    def _decorate(function):
        parameters = signature(function)

        @wraps(function)
        def _wrapped(*args, **kwargs):
            bound = parameters.bind(*args, **kwargs)
            bound.apply_defaults()
            values = bound.arguments
            input_paths = [values[name] for name in inputs if values.get(name) is not None]
            output_paths = [values[name] for name in outputs if values.get(name) is not None]
            scene_basename = values.get("scene_basename") or os.path.splitext(os.path.basename(str((input_paths or output_paths or [step])[0])))[0]
            with _processing_step(
                step, scene_basename, input_paths, output_paths,
                enabled=values.get("log_to_console", False),
                index=values.get("scene_index", 1), total=values.get("scene_total", 1),
                allow_nested=allow_nested,
            ):
                return function(*args, **kwargs)
        return _wrapped
    return _decorate


def _log(
    message: str,
    *,
    enabled: bool = False,
    step: str | None = None,
    scene_basename: str | None = None,
) -> None:
    """Emit a concise console log message.
    Args:
        message: Log message text.
        enabled: Whether logging is enabled.
        step: Optional step name for the log prefix.
        scene_basename: Optional scene basename for the log prefix.
    Returns:
        None.
    """
    if not enabled:
        return
    scene_basename = (scene_basename or _active_scene.get() or "").strip()
    step = step or _active_step.get() or ("workflow" if scene_basename else None)
    if scene_basename:
        prefix = f"[{scene_basename} {step}] "
    elif step:
        prefix = f"[{step}] "
    else:
        prefix = ""
    print(f"{prefix}{message}", flush=True)


__all__ = []
