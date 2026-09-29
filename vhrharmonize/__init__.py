from importlib import import_module
from typing import Any

_EXPORTS = {
    "joint_coregistration": (".plugins.joint_coregistration", "joint_coregistration"),
    "global_regression": (".plugins.global_regression", "global_regression"),
    "local_block_adjustment": (".plugins.local_block_adjustment", "local_block_adjustment"),
    "align_rasters": (".plugins.align_rasters", "align_rasters"),
    "create_footprints": (".plugins.create_footprints", "create_footprints"),
    "postprocess_footprints": (".plugins.postprocess_footprints", "postprocess_footprints"),
    "voronoi_center_seamline": (".plugins.voronoi_center_seamline", "voronoi_center_seamline"),
    "weighted_seamline": (".plugins.weighted_seamline", "weighted_seamline"),
    "markov_triangles": (".plugins.markov_triangles", "markov_triangles"),
    "mask_rasters": (".plugins.mask_rasters", "mask_rasters"),
    "merge_rasters": (".plugins.merge_rasters", "merge_rasters"),
    "merge_vectors": (".plugins.merge_vectors", "merge_vectors"),
    "band_math": (".plugins.band_math", "band_math"),
    "create_cloud_mask_with_omnicloudmask": (".plugins.create_cloud_mask_with_omnicloudmask", "create_cloud_mask_with_omnicloudmask"),
    "process_raster_values_to_vector_polygons": (".plugins.process_raster_values_to_vector_polygons", "process_raster_values_to_vector_polygons"),
    "compute_overviews": (".plugins.compute_overviews", "compute_overviews"),
    "search_paths": (".plugins.search_paths", "search_paths"),
    "create_paths": (".plugins.create_paths", "create_paths"),
    "match_paths": (".plugins.match_paths", "match_paths"),

    "load_workflow": (".workflow.api", "load_workflow"),
    "run_workflow": (".workflow.api", "run_workflow"),
    "run_plugin": (".workflow.api", "run_plugin"),
    "ProgressSnapshot": (".progress", "ProgressSnapshot"),
    "ProgressCallback": (".progress", "ProgressCallback"),
    "TimingEvent": (".progress", "TimingEvent"),
    "TimingCallback": (".progress", "TimingCallback"),
    "StatisticsRecorder": (".statistics", "StatisticsRecorder"),
    "summarize_statistics": (".statistics", "summarize_statistics"),
    "read_progress_snapshot": (".progress", "read_progress_snapshot"),
    "render_progress": (".progress", "render_progress"),
    "get_slurm_progress": (".slurm", "get_slurm_progress"),
    "prepare_slurm_plan": (".slurm", "prepare_slurm_plan"),
    "upload_slurm_files": (".slurm", "upload_slurm_files"),
    "start_slurm_job": (".slurm", "start_slurm_job"),
    "update_status_slurm_file": (".slurm", "update_status_slurm_file"),
    "stop_slurm_job": (".slurm", "stop_slurm_job"),
    "close_hpc_connection": (".slurm", "close_hpc_connection"),
    "download_slurm_outputs": (".slurm", "download_slurm_outputs"),
    "materialize_geometry": (".io.metadata", "materialize_geometry"),
    "read_metadata": (".io.metadata", "read_metadata"),
    "Workflow": (".workflow.engine", "Workflow"),
    "get_image_percentile_value": (".io.geospatial", "get_image_percentile_value"),
    "qgis_gcps_to_geojson": (
        ".plugins.orthorectification",
        "qgis_gcps_to_geojson",
    ),
    "gcp_refined_rpc_orthorectification": (
        ".plugins.orthorectification",
        "gcp_refined_rpc_orthorectification",
    ),
    "pansharpen_image": (".plugins.pansharpen", "pansharpen_image"),
    "run_flaash": (".plugins.atmospheric_correction", "run_flaash"),
    "run_py6s": (".plugins.atmospheric_correction", "run_py6s"),
    "cloudmask_raster": (".plugins.cloud_mask", "cloudmask_raster"),
    "apply_binary_cloud_mask_to_image": (
        ".plugins.cloud_mask",
        "apply_binary_cloud_mask_to_image",
    ),
    "resolve_output_resolution_for_crs": (
        ".plugins.orthorectification",
        "resolve_output_resolution_for_crs",
    ),
    "align_image_pair": (".plugins.alignment", "align_image_pair"),
    "AlignmentResult": (".plugins.alignment", "AlignmentResult"),
}

__all__ = list(_EXPORTS.keys())


def __getattr__(name: str) -> Any:
    """Lazily resolve a public package export.
    Args:
        name: Export name requested from the package namespace.
    Returns:
        The resolved export object.
    """
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = _EXPORTS[name]
    value = getattr(import_module(module_name, __name__), attr_name)
    globals()[name] = value
    return value
