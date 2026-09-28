"""create_paths: thin adapter for the installed SpectralMatch function."""

from .base import FunctionPlugin


def _upstream():
    from spectralmatch import create_paths

    return create_paths


def create_paths(*args, **kwargs):
    """Call SpectralMatch's create_paths with its native parameters."""

    return _upstream()(*args, **kwargs)


create_paths.__parameter_sources__ = (_upstream,)


class CreatePaths(FunctionPlugin):
    target = "vhrharmonize.plugins.create_paths:create_paths"
    scope = "aggregate"
