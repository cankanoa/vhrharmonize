"""global_regression: thin adapter for the installed SpectralMatch function."""

from .base import FunctionPlugin
from vhrharmonize.io.progress import call_with_progress, reports_progress


def _upstream():
    from spectralmatch import global_regression

    return global_regression


@reports_progress
def global_regression(*args, **kwargs):
    """Call SpectralMatch's global_regression with its native parameters."""

    for name in ("dask_scheduler", "specify_model_images", "vector_mask", "window_scales"):
        if isinstance(kwargs.get(name), list):
            kwargs[name] = tuple(kwargs[name])
    return call_with_progress(_upstream(), *args, **kwargs)


global_regression.__parameter_sources__ = (_upstream,)


class GlobalRegression(FunctionPlugin):
    target = "vhrharmonize.plugins.global_regression:global_regression"
    scope = "aggregate"
    input_dependency_paths = frozenset(["input_images", "load_adjustments", "pif_load_ties"])
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset(["output_images", "save_adjustments"])
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
