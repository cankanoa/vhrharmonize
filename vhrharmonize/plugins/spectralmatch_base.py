"""SpectralMatch-specific configuration for workflow-owned concurrency."""

from .base import FunctionPlugin


class SpectralMatchPlugin(FunctionPlugin):
    options = {"skip_single_image"}
    writing_postfix_parameter = "validity_writing_postfix"

    def output_finalizer(self, callback):
        """Keep core validation/overviews at the native writer's publication point."""
        from spectralmatch.utils_output import output_finalizer
        return output_finalizer(callback)

    def run(self, *, params, shared):
        """Optionally pass a one-image group through a matching-only step."""
        from pathlib import Path
        import shutil

        params = dict(params)
        skip = params.pop("skip_single_image", False)
        if not isinstance(skip, bool):
            raise ValueError("skip_single_image must be a boolean")
        images = params.get("input_images")
        if skip and isinstance(images, list) and len(images) == 1 and not params.get("fixed_images"):
            self.arguments(self.function(), params, shared)
            outputs = params.get("output_images")
            if not isinstance(outputs, list) or len(outputs) != 1:
                raise ValueError("skip_single_image requires one explicit output_images path")
            source, output = Path(images[0]).resolve(), Path(outputs[0]).resolve()
            if source == output:
                raise ValueError("Single-image output must differ from its input")
            output.parent.mkdir(parents=True, exist_ok=True)
            from spectralmatch.utils_output import output_context, current_output_context, writing_outputs
            _, finalizer = current_output_context()
            with output_context(params.get("validity_writing_postfix"), finalizer):
                with writing_outputs(str(output)) as (destination,):
                    shutil.copy2(source, destination)
            return [str(output)]
        return super().run(params=params, shared=shared)

    def worker_output_parameters(self, params):
        if params.get("skip_single_image") and len(params.get("input_images", [])) == 1 and not params.get("fixed_images"):
            return frozenset()
        return frozenset({"output_images", "output_vectors", "output_image_paths"})

    def worker_output_paths(self, function, args, kwargs):
        if (getattr(function, "__module__", None) == "spectralmatch.utils_logging"
                and getattr(function, "__name__", None) == "_run_logged_image" and len(args) >= 3):
            from vhrharmonize.workflow.output_files import paths
            return paths(args[2])
        return ()

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
