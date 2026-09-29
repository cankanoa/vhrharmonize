"""Core-only timing hooks; no processing plugin instrumentation is required."""

from functools import wraps
from contextlib import contextmanager
from time import monotonic, time_ns


@contextmanager
def timed_preflight(workflow, step, *, scene="all scenes", weight=1):
    """Record calls needed to discover scenes before normal worker reporting starts."""
    started_ns, started = time_ns(), monotonic()
    status = "completed"
    attributes = {"vhr.plugin": step["plugin"] or "", "vhr.scene": str(scene),
                  "vhr.scene_units": weight, "vhr.phase": "preflight", "vhr.measured": True}
    try:
        yield
    except BaseException as exc:
        status = "failed"
        attributes["error.type"] = type(exc).__name__
        raise
    finally:
        workflow._emit_timing("task", step["name"], started_ns, monotonic() - started,
                              status=status, attributes=attributes)


def timed_core(stage):
    def decorate(function):
        @wraps(function)
        def measured(self, *args, **kwargs):
            if (stage == "planning" and self._planned
                    or stage == "cleanup" and not kwargs.get("final", False)):
                return function(self, *args, **kwargs)
            started_ns, started = time_ns(), monotonic()
            status, attributes = "completed", {}
            try:
                return function(self, *args, **kwargs)
            except BaseException as exc:
                status = "failed"
                attributes = {"error.type": type(exc).__name__}
                raise
            finally:
                if hasattr(self, "_timing_lock"):
                    self._emit_timing("core", stage, started_ns, monotonic() - started,
                                      status=status, attributes=attributes)
        return measured
    return decorate
