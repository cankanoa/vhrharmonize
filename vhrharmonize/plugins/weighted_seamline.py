"""weighted_seamline: thin adapter for the installed SpectralMatch function."""

from .base import FunctionPlugin
from vhrharmonize.io.progress import call_with_progress, reports_progress


def _upstream():
    from spectralmatch import weighted_seamline

    return weighted_seamline


@reports_progress
def weighted_seamline(*args, **kwargs):
    """Call SpectralMatch's weighted_seamline with its native parameters."""

    return call_with_progress(_upstream(), *args, **kwargs)


weighted_seamline.__parameter_sources__ = (_upstream,)


class WeightedSeamline(FunctionPlugin):
    target = "vhrharmonize.plugins.weighted_seamline:weighted_seamline"
    scope = "aggregate"
    input_dependency_paths = frozenset(["input_images", "input_polygons"])
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset(["output_mask"])
    output_dependency_paths = output_path_resolution_paths
    output_parent_creation_paths = output_path_resolution_paths
    output_collision_check_paths = output_path_resolution_paths
    output_validation_paths = output_path_resolution_paths
    output_invalid_removal_paths = output_path_resolution_paths
    output_temporary_cleanup_paths = output_path_resolution_paths
    output_hpc_staging_paths = output_path_resolution_paths
    output_hpc_download_paths = output_path_resolution_paths
    output_reuse_paths = output_path_resolution_paths
