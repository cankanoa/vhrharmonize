"""Publish completed plugin outputs without exposing in-progress final filenames."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from contextlib import contextmanager, nullcontext
from time import monotonic, time_ns
import os

from vhrharmonize.io.validation import _existing_output_failures
from .values import remap_paths


def map_paths(value, mapping):
    """Remap worker argument containers, including positional tuples and Path objects."""
    if isinstance(value, tuple):
        return tuple(map_paths(v, mapping) for v in value)
    if isinstance(value, list):
        return [map_paths(v, mapping) for v in value]
    if isinstance(value, dict):
        return {k: map_paths(v, mapping) for k, v in value.items()}
    if isinstance(value, Path):
        return Path(remap_paths(str(value), mapping))
    return remap_paths(value, mapping)


def paths(value):
    if isinstance(value, (str, Path)):
        return [str(value)]
    if isinstance(value, (list, tuple)):
        return [p for item in value for p in paths(item)]
    return []


def writing_path(filename, postfix):
    path = Path(filename)
    return str(path.with_name(path.stem + postfix + path.suffix))


def companions(filename):
    """GDAL sidecars whose names follow the raster destination."""
    if Path(filename).suffix.lower() not in {".tif", ".tiff", ".dat", ".img", ".vrt"}:
        return []
    result = [filename + suffix for suffix in (".aux.xml", ".msk", ".ovr", ".msk.ovr")]
    if Path(filename).suffix.lower() in {".dat", ".img"}:
        result.append(str(Path(filename).with_suffix(".hdr")))
    return result


@dataclass(frozen=True)
class OutputPolicy:
    postfix: str | None = "_writing"
    check_validity: str | bool = "before_reuse"
    tile_threads: int | str | None = 1
    raster_grid: int | bool = 255
    check_stop: str = "readable"
    validation_paths: tuple[str, ...] = ()
    overview_paths: tuple[str, ...] = ()
    overview_scales: object = None
    log_to_console: bool = False
    step: str = "output"
    overview_reporters: dict = field(default_factory=dict)

    @contextmanager
    def overview(self, filename):
        from vhrharmonize.io.progress import progress_context
        reporter = self.overview_reporters.get(filename)
        start, start_ns = monotonic(), time_ns()
        if reporter:
            reporter.event("start", started_ns=start_ns)
        try:
            with progress_context(reporter, reporter.message) if reporter else nullcontext():
                yield
        except BaseException as exc:
            if reporter:
                reporter.event("failed", duration=monotonic() - start, started_ns=start_ns,
                               error_type=type(exc).__name__)
            raise
        else:
            if reporter:
                reporter.event("computed", duration=monotonic() - start, started_ns=start_ns)
                reporter.event("accepted")

    def finish(self, outputs, mapping):
        """Build overviews, validate, then publish each closed output at its final name."""
        from vhrharmonize.io.geospatial import calculate_raster_overviews

        for final in outputs:
            work = mapping.get(final, final)
            if final in self.overview_paths and Path(final).suffix.lower() in {".tif", ".tiff"}:
                with self.overview(final):
                    calculate_raster_overviews(work, self.overview_scales, log_to_console=self.log_to_console)
        checked = [mapping.get(p, p) for p in outputs if p in self.validation_paths]
        failures = _existing_output_failures(
            checked, check_validity=self.check_validity == "around",
            tile_threads=self.tile_threads,
            validity_raster_grid=self.raster_grid, validity_raster_check_stop=self.check_stop,
            log_to_console=False, step=self.step,
        )
        # Even unvalidated declared outputs must exist before publishing.
        failures.update({mapping.get(p, p): "missing" for p in outputs
                         if not os.path.exists(mapping.get(p, p))})
        if failures:
            raise RuntimeError(f"{self.step} did not produce valid declared outputs: {failures}")
        for final, work in mapping.items():
            if os.path.exists(work + "-wal") and os.path.getsize(work + "-wal"):
                raise RuntimeError(f"Output database must be closed/checkpointed before publishing: {work}")
            # Publish companions first, the main file last. A missing final file
            # remains unambiguously incomplete if a move fails partway through.
            for src, dst in zip(companions(work), companions(final)):
                if os.path.isfile(src):
                    os.replace(src, dst)
                elif os.path.isfile(dst):
                    os.unlink(dst)  # Do not attach stale masks/overviews to the new image.
            os.replace(work, final)

    def run(self, function, args=(), kwargs=None, *, outputs=()):
        outputs = tuple(dict.fromkeys(outputs))
        mapping = {p: writing_path(p, self.postfix) for p in outputs if not os.path.isdir(p)} if self.postfix else {}
        if set(mapping.values()) & set(outputs):
            raise ValueError("Writing filenames collide with declared final outputs")
        for work in mapping.values():
            # An interrupted attempt is never reused under the writing name.
            for candidate in (work, *companions(work), work + "-wal", work + "-shm", work + "-journal"):
                if os.path.isfile(candidate):
                    os.unlink(candidate)
        result = function(*map_paths(args, mapping), **map_paths(kwargs or {}, mapping))
        self.finish(outputs, mapping)
        return map_paths(result, {work: final for final, work in mapping.items()})


def _run_output_job(function, args, kwargs, policy, outputs):
    return policy.run(function, args, kwargs, outputs=outputs)


class OutputExecutor:
    """Finalize each image in its worker before the future reports success."""
    def __init__(self, executor, plugin, policy, outputs):
        self.executor, self.plugin, self.policy = executor, plugin, policy
        self.outputs = set(outputs)

    def __getattr__(self, name):
        return getattr(self.executor, name)

    def submit(self, function, /, *args, **kwargs):
        outputs = self.outputs.intersection(self.plugin.worker_output_paths(function, args, kwargs))
        if not outputs:
            return self.executor.submit(function, *args, **kwargs)
        # Keep the original function visible to failure attribution in the parent.
        if hasattr(self.executor, "submit_output"):
            return self.executor.submit_output(function, args, kwargs, self.policy, tuple(sorted(outputs)))
        return self.executor.submit(_run_output_job, function, args, kwargs, self.policy, tuple(sorted(outputs)))
