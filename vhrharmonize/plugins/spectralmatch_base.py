"""SpectralMatch-specific configuration for workflow-owned concurrency."""

from .base import FunctionPlugin


class SpectralMatchPlugin(FunctionPlugin):
    def runtime_parameters(self, *, tile_threads):
        values = super().runtime_parameters(tile_threads=tile_threads)
        if self.executor_parameter is not None:
            # Scene calls run serially inside the workflow worker; aggregate
            # calls borrow its executor. Neither creates a nested image pool.
            values.update(image_threads=None, dask_scheduler=None,
                          concurrent_processing_backend="process_pool")
        return values
