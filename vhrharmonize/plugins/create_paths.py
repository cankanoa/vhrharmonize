"""create_paths: thin adapter for the installed SpectralMatch function."""

from .base import FunctionPlugin
from vhrharmonize.io.progress import call_with_progress, reports_progress


def _upstream():
    from spectralmatch import create_paths

    return create_paths


@reports_progress
def create_paths(*args, **kwargs):
    """Call SpectralMatch's create_paths with its native parameters."""

    return call_with_progress(_upstream(), *args, **kwargs)


create_paths.__parameter_sources__ = (_upstream,)


class CreatePaths(FunctionPlugin):
    target = "vhrharmonize.plugins.create_paths:create_paths"
    scope = "aggregate"
