"""match_paths: thin adapter for the installed SpectralMatch function."""

from .base import FunctionPlugin
from vhrharmonize.io.progress import call_with_progress, reports_progress


def _upstream():
    from spectralmatch import match_paths

    return match_paths


@reports_progress
def match_paths(*args, **kwargs):
    """Call SpectralMatch's match_paths with its native parameters."""

    return call_with_progress(_upstream(), *args, **kwargs)


match_paths.__parameter_sources__ = (_upstream,)


class MatchPaths(FunctionPlugin):
    target = "vhrharmonize.plugins.match_paths:match_paths"
    scope = "aggregate"
