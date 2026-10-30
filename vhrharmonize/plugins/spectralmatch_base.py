"""SpectralMatch-specific configuration for workflow-owned concurrency."""

from .base import FunctionPlugin


class SpectralMatchPlugin(FunctionPlugin):
    def worker_image_paths(self, function, args, kwargs):
        if (getattr(function, "__module__", None) == "spectralmatch.utils_logging"
                and getattr(function, "__name__", None) == "_run_logged_image" and len(args) >= 2):
            return args[1]
        return ()

    def runtime_parameters(self, *, tile_threads):
        values = super().runtime_parameters(tile_threads=tile_threads)
        if self.executor_parameter is not None:
            # Scene calls run serially inside the workflow worker; aggregate
            # calls borrow its executor. Neither creates a nested image pool.
            values.update(image_threads=None, dask_scheduler=None,
                          concurrent_processing_backend="process_pool")
        return values
