"""align_rasters: thin adapter for the installed SpectralMatch function."""

from .base import FunctionPlugin


def _upstream():
    from spectralmatch import align_rasters

    return align_rasters


def align_rasters(*args, **kwargs):
    """Call SpectralMatch's align_rasters with its native parameters."""

    for name in ("dask_scheduler", "resolution"):
        if isinstance(kwargs.get(name), list):
            kwargs[name] = tuple(kwargs[name])
    return _upstream()(*args, **kwargs)


align_rasters.__parameter_sources__ = (_upstream,)


class AlignRasters(FunctionPlugin):
    target = "vhrharmonize.plugins.align_rasters:align_rasters"
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
