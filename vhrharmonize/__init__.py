"""Public package API with lazy-loaded exports."""

from importlib import import_module

_EXPORTS = {
    "find_files": (".providers.worldview.files", "find_files"),
    "find_roots": (".providers.worldview.files", "find_roots"),
    "get_metadata_from_files": (".providers.worldview.files", "get_metadata_from_files"),
    "shp_to_gpkg": (".io.geospatial", "shp_to_gpkg"),
    "get_image_largest_value": (".io.geospatial", "get_image_largest_value"),
    "get_image_percentile_value": (".io.geospatial", "get_image_percentile_value"),
    "wsl_to_windows_path": (".io.geospatial", "wsl_to_windows_path"),
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
    "create_cloud_mask_with_omnicloudmask": (
        ".preprocess.cloudmasking",
        "create_cloud_mask_with_omnicloudmask",
    ),
    "apply_binary_cloud_mask_to_image": (
        ".preprocess.cloudmasking",
        "apply_binary_cloud_mask_to_image",
    ),
    "tile_image": (".preprocess.tiling", "tile_image"),
    "run_elastix_registration": (
        ".preprocess.registration",
        "run_elastix_registration",
    ),
    "estimate_elastix_transform": (
        ".preprocess.registration",
        "estimate_elastix_transform",
    ),
    "apply_elastix_transform": (
        ".preprocess.registration",
        "apply_elastix_transform",
    ),
}

__all__ = list(_EXPORTS.keys())


def __getattr__(name):
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = _EXPORTS[name]
    value = getattr(import_module(module_name, __name__), attr_name)
    globals()[name] = value
    return value
