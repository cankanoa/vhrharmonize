"""Renderer-independent progress callbacks using tqdm's ``format_dict`` fields."""

from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from functools import wraps
from inspect import Parameter, signature
from threading import RLock, local
from time import monotonic
from uuid import uuid4
import sys

from tqdm import tqdm


_callback = ContextVar("progress_callback", default=None)
_messages = ContextVar("progress_messages", default=None)
_stream_lock = RLock()
_stream_users = 0


def current_callback():
    return _callback.get()


def call_with_progress(function, *args, **kwargs):
    """Forward supported callbacks and silence unmanaged backend worker bars."""
    callback = current_callback()
    if callback is not None:
        parameters = signature(function).parameters
        if "progress_callback" in parameters:
            kwargs.setdefault("progress_callback", callback)
        if "log_to_console" in parameters:
            kwargs["log_to_console"] = False
    return function(*args, **kwargs)


class CallbackTqdm(tqdm):
    """Count and estimate rates with tqdm, forwarding snapshots without printing."""

    def __init__(self, *args, callback, **kwargs):
        self.callback = callback
        self.operation = uuid4().hex
        kwargs.update(gui=True, disable=False)
        super().__init__(*args, **kwargs)
        self.display()

    def display(self, *args, **kwargs):
        self.callback(**self.format_dict, operation=self.operation)

    def close(self):
        if not self.disable:
            self.display()
        super().close()


def progress(iterable=None, *, callback=None, **kwargs):
    """Use an inherited callback, or a regular (silent by default) tqdm bar.

    Callbacks receive keyword fields from ``tqdm.format_dict``: in particular
    ``n``, ``total``, ``prefix``, ``unit``, ``elapsed`` and ``rate``. ``operation``
    identifies a particular bar, allowing its count to restart for a new phase.
    """
    callback = callback if callback is not None else current_callback()
    if callback is not None:
        return CallbackTqdm(iterable, callback=callback, **kwargs)
    kwargs.setdefault("disable", True)
    return tqdm(iterable, **kwargs)


@contextmanager
def progress_context(callback=None, messages=None):
    token = _callback.set(callback)
    message_token = _messages.set(messages)
    try:
        yield
    finally:
        _messages.reset(message_token)
        _callback.reset(token)


@contextmanager
def operation(description, callback=None):
    """Report lifecycle for opaque operations without inventing a percentage."""
    callback = callback if callback is not None else current_callback()
    if callback is None:
        yield
        return
    identity, started = uuid4().hex, monotonic()
    fields = dict(prefix=description, unit="call", operation=identity, rate=None)
    callback(n=0, total=None, elapsed=0, **fields)
    try:
        yield
    except BaseException:
        callback(n=0, total=None, elapsed=monotonic() - started, status="failed", **fields)
        raise
    else:
        callback(n=1, total=1, elapsed=monotonic() - started, **fields)


def reports_progress(function=None, *, worker_progress=False):
    """Add an optional Python callback, preserving the function's public API."""
    if function is None:
        return lambda target: reports_progress(target, worker_progress=worker_progress)
    if getattr(function, "__reports_progress__", False):
        return function
    parameters = signature(function)
    accepts_callback = "progress_callback" in parameters.parameters

    @wraps(function)
    def wrapped(*args, progress_callback=None, **kwargs):
        callback = progress_callback if progress_callback is not None else current_callback()
        if callback is not None and not callable(callback):
            raise TypeError("progress_callback must be callable or None")
        token = _callback.set(callback)
        try:
            with operation(function.__name__, callback):
                if accepts_callback:
                    kwargs["progress_callback"] = callback
                return function(*args, **kwargs)
        finally:
            _callback.reset(token)

    if not accepts_callback:
        items = list(parameters.parameters.values())
        position = next((i for i, p in enumerate(items) if p.kind == p.VAR_KEYWORD), len(items))
        items.insert(position, Parameter("progress_callback", Parameter.KEYWORD_ONLY, default=None))
        wrapped.__signature__ = parameters.replace(parameters=items)
    wrapped.__reports_progress__ = True
    wrapped.__worker_progress__ = worker_progress or accepts_callback
    return wrapped


def gdal_progress(description):
    """Adapt GDAL's fractional callback to the same tqdm field contract."""
    callback = current_callback()
    if callback is None:
        return None
    started, identity = monotonic(), uuid4().hex
    last = [-1.0]

    def report(fraction, message, data):
        elapsed = monotonic() - started
        if fraction in (0, 1) or elapsed - last[0] >= 0.1:
            last[0] = elapsed
            callback(n=fraction * 100, total=100, prefix=description, unit="%",
                     elapsed=elapsed, rate=fraction * 100 / elapsed if elapsed else None,
                     operation=identity)
        return 1

    return report


def raster_windows(dataset, band=1, *, desc="Processing tiles"):
    """Track raster blocks without materializing a potentially large window list."""
    height, width = dataset.block_shapes[band - 1]
    total = ((dataset.height + height - 1) // height) * ((dataset.width + width - 1) // width)
    return progress(dataset.block_windows(band), total=total, desc=desc, unit="tiles")


class _MessageStream:
    """Route Python output by execution context, including concurrent threads."""

    def __init__(self, original):
        self.original = original
        self.buffers = local()

    def write(self, text):
        sink = _messages.get()
        if sink is None:
            return self.original.write(text)
        value = getattr(self.buffers, "text", "") + text
        lines = value.replace("\r", "\n").split("\n")
        self.buffers.text = lines.pop()[-4096:]
        for line in lines:
            if line.strip():
                sink(line[-4096:])
        return len(text)

    def flush(self):
        sink = _messages.get()
        value = getattr(self.buffers, "text", "")
        if sink is not None and value.strip():
            self.buffers.text = ""
            sink(value)
        self.original.flush()

    def isatty(self):
        return False if _messages.get() is not None else self.original.isatty()

    def __getattr__(self, name):
        return getattr(self.original, name)


@contextmanager
def capture_messages():
    """Install one context-aware stream router per process, restoring on exit."""
    global _stream_users
    with _stream_lock:
        if _stream_users == 0:
            sys.stdout = _MessageStream(sys.stdout)
            sys.stderr = _MessageStream(sys.stderr)
        _stream_users += 1
    try:
        yield
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        with _stream_lock:
            _stream_users -= 1
            if _stream_users == 0:
                sys.stdout = sys.stdout.original
                sys.stderr = sys.stderr.original
