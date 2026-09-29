"""compute_overviews: thin adapter for the installed SpectralMatch function."""

from .base import FunctionPlugin
from vhrharmonize.io.progress import call_with_progress, reports_progress


def _upstream():
    from spectralmatch import compute_overviews

    return compute_overviews


@reports_progress
def compute_overviews(*args, **kwargs):
    """Call SpectralMatch's compute_overviews with its native parameters."""

    for name in ("dask_scheduler", "window_scales"):
        if isinstance(kwargs.get(name), list):
            kwargs[name] = tuple(kwargs[name])
    return call_with_progress(_upstream(), *args, **kwargs)


compute_overviews.__parameter_sources__ = (_upstream,)


class ComputeOverviews(FunctionPlugin):
    target = "vhrharmonize.plugins.compute_overviews:compute_overviews"
    scope = "aggregate"
    input_dependency_paths = frozenset(["input_images_paths"])
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset(["output_image_paths"])
    output_dependency_paths = output_path_resolution_paths
    output_target_paths = output_path_resolution_paths
    output_parent_creation_paths = output_path_resolution_paths
    output_collision_check_paths = output_path_resolution_paths
    output_validation_paths = output_path_resolution_paths
    output_invalid_removal_paths = output_path_resolution_paths
    output_temporary_cleanup_paths = output_path_resolution_paths
    output_hpc_staging_paths = output_path_resolution_paths
    output_hpc_download_paths = output_path_resolution_paths
    output_context_checkpoint_paths = frozenset(["output_image_paths"])
