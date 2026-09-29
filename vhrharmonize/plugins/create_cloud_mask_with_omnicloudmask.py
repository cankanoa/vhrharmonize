"""create_cloud_mask_with_omnicloudmask: thin adapter for the installed SpectralMatch function."""

from .base import FunctionPlugin
from vhrharmonize.io.progress import call_with_progress, reports_progress


def _upstream():
    from spectralmatch import create_cloud_mask_with_omnicloudmask

    return create_cloud_mask_with_omnicloudmask


@reports_progress
def create_cloud_mask_with_omnicloudmask(*args, **kwargs):
    """Call SpectralMatch's create_cloud_mask_with_omnicloudmask with its native parameters."""

    for name in ("dask_scheduler",):
        if isinstance(kwargs.get(name), list):
            kwargs[name] = tuple(kwargs[name])
    return call_with_progress(_upstream(), *args, **kwargs)


create_cloud_mask_with_omnicloudmask.__parameter_sources__ = (_upstream,)


class CreateCloudMaskWithOmnicloudmask(FunctionPlugin):
    target = "vhrharmonize.plugins.create_cloud_mask_with_omnicloudmask:create_cloud_mask_with_omnicloudmask"
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
