"""local_block_adjustment: thin adapter for the installed SpectralMatch function."""

from .base import FunctionPlugin
from vhrharmonize.io.progress import call_with_progress, reports_progress


def _upstream():
    from spectralmatch import local_block_adjustment

    return local_block_adjustment


@reports_progress
def local_block_adjustment(*args, **kwargs):
    """Call SpectralMatch's local_block_adjustment with its native parameters."""

    for name in (
        "dask_scheduler",
        "load_block_maps",
        "number_of_blocks",
        "override_bounds_canvas_coords",
        "save_block_maps",
        "vector_mask",
        "window_scales",
    ):
        if isinstance(kwargs.get(name), list):
            kwargs[name] = tuple(kwargs[name])
    return call_with_progress(_upstream(), *args, **kwargs)


local_block_adjustment.__parameter_sources__ = (_upstream,)


class LocalBlockAdjustment(FunctionPlugin):
    target = "vhrharmonize.plugins.local_block_adjustment:local_block_adjustment"
    scope = "aggregate"
    input_dependency_paths = frozenset(["input_images"])
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset(["output_images"])
    output_dependency_paths = output_path_resolution_paths
    output_target_paths = output_path_resolution_paths
    output_parent_creation_paths = output_path_resolution_paths
    output_collision_check_paths = output_path_resolution_paths
    output_validation_paths = output_path_resolution_paths
    output_invalid_removal_paths = output_path_resolution_paths
    output_temporary_cleanup_paths = output_path_resolution_paths
    output_hpc_staging_paths = output_path_resolution_paths
    output_hpc_download_paths = output_path_resolution_paths
    output_reuse_paths = output_path_resolution_paths
    output_context_checkpoint_paths = frozenset(["output_images"])
    output_overview_calculation_paths = frozenset(["output_images"])
