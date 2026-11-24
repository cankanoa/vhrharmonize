"""mask_rasters: thin adapter for the installed SpectralMatch function."""

from functools import lru_cache
import os

from .spectralmatch_base import SpectralMatchPlugin
from vhrharmonize.io.progress import call_with_progress, reports_progress


def _upstream():
    from spectralmatch import mask_rasters

    return mask_rasters


@reports_progress
def mask_rasters(*args, **kwargs):
    """Call SpectralMatch's mask_rasters with its native parameters."""

    for name in ("dask_scheduler", "vector_mask"):
        if isinstance(kwargs.get(name), list):
            kwargs[name] = tuple(kwargs[name])
    return call_with_progress(_upstream(), *args, **kwargs)


mask_rasters.__parameter_sources__ = (_upstream,)


@lru_cache(maxsize=32)
def _numbered_outputs(inputs, outputs, vector_mask, fingerprints):
    from spectralmatch.utils import _mask_raster_jobs

    return tuple(job[1] for job in _mask_raster_jobs(
        inputs, outputs, vector_mask, seperate_multipolygons=True))


class MaskRasters(SpectralMatchPlugin):
    def output_files(self, *, params):
        """Track numbered rasters instead of their unnumbered naming bases."""
        from vhrharmonize.workflow.values import contains_pending

        if params.get("seperate_multipolygons") is not True or contains_pending(params):
            return {}
        inputs, outputs = params.get("input_images"), params.get("output_images")
        if not isinstance(inputs, list) or not isinstance(outputs, list):
            return {}  # Native folder/template forms are for direct Python use.
        mask = params.get("vector_mask")
        mask = tuple(mask) if mask else None
        dependencies = [*inputs, *([mask[1]] if mask else [])]
        if not all(os.path.isfile(filename) for filename in dependencies):
            return {}  # A producing step may still need to create these inputs.
        fingerprints = tuple((os.stat(p).st_mtime_ns, os.stat(p).st_size) for p in dependencies)
        return {"output_images": list(_numbered_outputs(tuple(inputs), tuple(outputs), mask, fingerprints))}

    def can_reuse_outputs(self, *, params):
        # Running restores the flat return list needed by downstream aggregate
        # steps. Native resume still reuses each already completed numbered file.
        return not params.get("seperate_multipolygons", False)

    executor_parameter = "executor"
    tile_threads_parameter = "tile_threads"
    target = "vhrharmonize.plugins.mask_rasters:mask_rasters"
    scope = "aggregate"
    input_dependency_paths = frozenset(["input_images"])
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset(["output_images"])
    output_dependency_paths = output_path_resolution_paths
    output_parent_creation_paths = output_path_resolution_paths
    output_collision_check_paths = output_path_resolution_paths
    output_validation_paths = output_path_resolution_paths
    output_invalid_removal_paths = output_path_resolution_paths
    output_temporary_cleanup_paths = output_path_resolution_paths
    output_hpc_staging_paths = output_path_resolution_paths
    output_hpc_download_paths = output_path_resolution_paths
    output_reuse_paths = output_path_resolution_paths
    output_overview_calculation_paths = frozenset(["output_images"])
