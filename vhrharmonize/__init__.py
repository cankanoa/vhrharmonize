from importlib import import_module
from typing import Any

_EXPORTS = {
    "find_files": (".providers.worldview", "find_files"),
    "shp_to_gpkg": (".io.geospatial", "shp_to_gpkg"),
    "get_image_percentile_value": (".io.geospatial", "get_image_percentile_value"),
    "qgis_gcps_to_geojson": (
        ".preprocess.orthorectification",
        "qgis_gcps_to_geojson",
    ),
    "gcp_refined_rpc_orthorectification": (
        ".preprocess.orthorectification",
        "gcp_refined_rpc_orthorectification",
    ),
    "pansharpen_image": (".preprocess.pansharpening", "pansharpen_image"),
    "run_flaash": (".preprocess.atmospheric_correction", "run_flaash"),
    "run_py6s": (".preprocess.atmospheric_correction", "run_py6s"),
    "create_cloud_mask_with_omnicloudmask": (
        ".preprocess.cloudmasking",
        "create_cloud_mask_with_omnicloudmask",
    ),
    "cloudmask_raster": (".preprocess.cloudmasking", "cloudmask_raster"),
    "apply_binary_cloud_mask_to_image": (
        ".preprocess.cloudmasking",
        "apply_binary_cloud_mask_to_image",
    ),
    "resolve_output_resolution_for_crs": (
        ".preprocess.orthorectification",
        "resolve_output_resolution_for_crs",
    ),
    "radiometric_normalization": (
        ".preprocess.radiometric_normalization",
        "radiometric_normalization",
    ),
    "align_image_pair": (".preprocess.alignment", "align_image_pair"),
    "AlignmentResult": (".preprocess.alignment", "AlignmentResult"),
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
