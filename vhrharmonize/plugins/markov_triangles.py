"""markov_triangles: thin adapter for the installed SpectralMatch function."""

from .base import FunctionPlugin


def _upstream():
    from spectralmatch import markov_triangles

    return markov_triangles


def markov_triangles(*args, **kwargs):
    """Call SpectralMatch's markov_triangles with its native parameters."""

    for name in ("dask_scheduler",):
        if isinstance(kwargs.get(name), list):
            kwargs[name] = tuple(kwargs[name])
    return _upstream()(*args, **kwargs)


markov_triangles.__parameter_sources__ = (_upstream,)


class MarkovTriangles(FunctionPlugin):
    target = "vhrharmonize.plugins.markov_triangles:markov_triangles"
    scope = "aggregate"
    input_dependency_paths = frozenset(["input_images", "input_polygons", "foreground_path"])
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset(["output_mask"])
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
    output_context_checkpoint_paths = frozenset(["output_mask"])
