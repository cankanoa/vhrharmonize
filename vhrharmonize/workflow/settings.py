"""Core implementation of steps that only assign context values."""

from vhrharmonize.plugins.base import FunctionPlugin


class ContextStep(FunctionPlugin):
    def run(self, *, params, shared):
        return None
