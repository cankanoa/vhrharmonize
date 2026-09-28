"""search_paths: thin adapter for the installed SpectralMatch function."""

from .base import FunctionPlugin


def _upstream():
    from spectralmatch import search_paths

    return search_paths


def search_paths(*args, **kwargs):
    """Call SpectralMatch's search_paths with its native parameters."""

    for name in ("match_to_paths",):
        if isinstance(kwargs.get(name), list):
            kwargs[name] = tuple(kwargs[name])
    return _upstream()(*args, **kwargs)


search_paths.__parameter_sources__ = (_upstream,)


class SearchPaths(FunctionPlugin):
    target = "vhrharmonize.plugins.search_paths:search_paths"
    scope = "aggregate"
