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
    removes_scenes = True

    def scene_ids_to_remove(self, *, params, returned, scene_ids):
        """Remove images absent from the written footprints after this step runs."""
        import os
        import geopandas as gpd

        def identity(filename):
            return os.path.normcase(os.path.realpath(filename))

        frame = gpd.read_file(params["output_polygons"], layer=params.get("output_layer", "footprints"))
        retained = {identity(filename) for filename in frame["image_path"]}
        paths = params["input_images"]
        if not isinstance(paths, list) or len(paths) != len(scene_ids):
            raise ValueError("Scene removal requires one explicit input image per scene")
        return [scene_id for scene_id, filename in zip(scene_ids, paths)
                if identity(filename) not in retained]

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
