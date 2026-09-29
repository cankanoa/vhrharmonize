"""merge_vectors: thin adapter for the installed SpectralMatch function."""

from .base import FunctionPlugin
from vhrharmonize.io.progress import call_with_progress, reports_progress


def _upstream():
    from spectralmatch import merge_vectors

    return merge_vectors


@reports_progress
def merge_vectors(*args, **kwargs):
    """Call SpectralMatch's merge_vectors with its native parameters."""

    for name in ("create_name_attribute",):
        if isinstance(kwargs.get(name), list):
            kwargs[name] = tuple(kwargs[name])
    return call_with_progress(_upstream(), *args, **kwargs)


merge_vectors.__parameter_sources__ = (_upstream,)


class MergeVectors(FunctionPlugin):
    target = "vhrharmonize.plugins.merge_vectors:merge_vectors"
    scope = "aggregate"
    input_dependency_paths = frozenset(["input_vectors"])
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset(["merged_vector_path"])
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
    output_context_checkpoint_paths = frozenset(["merged_vector_path"])
