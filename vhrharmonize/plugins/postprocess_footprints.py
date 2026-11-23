"""postprocess_footprints: thin adapter for the installed SpectralMatch function."""

from .spectralmatch_base import SpectralMatchPlugin
from vhrharmonize.io.progress import call_with_progress, reports_progress


def _upstream():
    from spectralmatch import postprocess_footprints

    return postprocess_footprints


@reports_progress
def postprocess_footprints(*args, **kwargs):
    """Call SpectralMatch's postprocess_footprints with its native parameters."""

    for name in ("dask_scheduler",):
        if isinstance(kwargs.get(name), list):
            kwargs[name] = tuple(kwargs[name])
    return call_with_progress(_upstream(), *args, **kwargs)


postprocess_footprints.__parameter_sources__ = (_upstream,)


class PostprocessFootprints(SpectralMatchPlugin):
    executor_parameter = "executor"
    removes_scenes = True

    def can_reuse_outputs(self, *, params):
        """Do not reuse an unfinished streaming footprint file."""
        import os

        output = params.get("output_polygons")
        return isinstance(output, str) and not os.path.exists(output + ".incomplete")

    def scene_ids_to_remove(self, *, params, returned, scene_ids):
        """Select scenes from fresh or reused output, without reopening input polygons.

        The workflow's footprint metadata must preserve ``scene_id``. Raster
        filenames may have changed since polygonization, so they are not used
        to infer scene identity.
        """
        import geopandas as gpd

        frame = gpd.read_file(params["output_polygons"], layer=params.get("output_layer", "footprints"))
        frame = frame.loc[frame.geometry.notna() & ~frame.geometry.is_empty]
        if frame.empty:
            raise ValueError("No coverage remains after postprocess_footprints; all scene footprints were removed.")
        if "scene_id" not in frame:
            raise ValueError(
                "postprocess_footprints scene removal requires a scene_id attribute; "
                "include scene_id in create_footprints metadata_records and regenerate the footprints."
            )
        identifiers = frame["scene_id"]
        if identifiers.isna().any() or identifiers.astype(str).str.strip().eq("").any():
            raise ValueError("Postprocessed footprints contain missing scene_id values.")
        retained = set(identifiers.astype(str))
        unknown = retained - set(scene_ids)
        if unknown:
            raise ValueError(
                "Postprocessed footprints contain scene IDs outside the current workflow: "
                + ", ".join(sorted(unknown))
                + ". Regenerate the footprints for the current inputs."
            )
        return [scene_id for scene_id in scene_ids if scene_id not in retained]

    target = "vhrharmonize.plugins.postprocess_footprints:postprocess_footprints"
    scope = "aggregate"
    input_dependency_paths = frozenset(["input_polygons"])
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
    output_reuse_paths = output_path_resolution_paths
