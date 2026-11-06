"""create_footprints: thin adapter for the installed SpectralMatch function."""

from .spectralmatch_base import SpectralMatchPlugin
from vhrharmonize.io.progress import call_with_progress, reports_progress


def _upstream():
    from spectralmatch import create_footprints

    return create_footprints


@reports_progress
def create_footprints(*args, **kwargs):
    """Call SpectralMatch's create_footprints with its native parameters."""

    for name in ("dask_scheduler",):
        if isinstance(kwargs.get(name), list):
            kwargs[name] = tuple(kwargs[name])
    return call_with_progress(_upstream(), *args, **kwargs)


create_footprints.__parameter_sources__ = (_upstream,)


class CreateFootprints(SpectralMatchPlugin):
    executor_parameter = "executor"
    record_filter_parameters = frozenset({"outline_path", "image_rank_function"})

    def retained_record_indices(self, params, returned):
        """Keep workflow scene order while omitting images with no output coverage."""
        import geopandas as gpd

        retained = set(gpd.read_file(returned, layer=params.get("output_layer"))["image_path"])
        return [index for index, path in enumerate(params["input_images"]) if path in retained]

    target = "vhrharmonize.plugins.create_footprints:create_footprints"
    scope = "aggregate"
    input_dependency_paths = frozenset(["input_images", "metadata_csv", "outline_path"])
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset(["output_polygons"])
    output_dependency_paths = output_path_resolution_paths
    output_parent_creation_paths = output_path_resolution_paths
    output_collision_check_paths = output_path_resolution_paths
    output_validation_paths = output_path_resolution_paths
    output_invalid_removal_paths = output_path_resolution_paths
    output_temporary_cleanup_paths = output_path_resolution_paths
    output_hpc_staging_paths = output_path_resolution_paths
    output_hpc_download_paths = output_path_resolution_paths
    # SpectralMatch checks streaming completeness before applying native resume.
    output_reuse_paths = frozenset()
