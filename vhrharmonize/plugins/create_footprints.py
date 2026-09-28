"""create_footprints: thin adapter for the installed SpectralMatch function."""

from .base import FunctionPlugin


def _upstream():
    from spectralmatch import create_footprints

    return create_footprints


def create_footprints(*args, **kwargs):
    """Call SpectralMatch's create_footprints with its native parameters."""

    for name in ("dask_scheduler",):
        if isinstance(kwargs.get(name), list):
            kwargs[name] = tuple(kwargs[name])
    return _upstream()(*args, **kwargs)


create_footprints.__parameter_sources__ = (_upstream,)


class CreateFootprints(FunctionPlugin):
    target = "vhrharmonize.plugins.create_footprints:create_footprints"
    scope = "aggregate"
    input_dependency_paths = frozenset(["input_images", "metadata_csv"])
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset(["output_polygons"])
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
    output_context_checkpoint_paths = frozenset(["output_polygons"])
