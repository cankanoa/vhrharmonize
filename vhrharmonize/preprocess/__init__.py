"""Preprocessing step modules."""

from importlib import import_module
from typing import Any

_EXPORTS = {
    "log": (".helpers", "log"),
    "atmospheric_correction": (".atmospheric_correction", "atmospheric_correction"),
    "run_py6s": (".atmospheric_correction", "run_py6s"),
    "run_flaash": (".atmospheric_correction", "run_flaash"),
    "gcp_refined_rpc_orthorectification": (
        ".orthorectification",
        "gcp_refined_rpc_orthorectification",
    ),
    "resolve_output_resolution_for_crs": (
        ".orthorectification",
        "resolve_output_resolution_for_crs",
    ),
    "pansharpen_image": (".pansharpening", "pansharpen_image"),
    "cloudmask_raster": (".cloudmasking", "cloudmask_raster"),
    "create_cloud_mask_with_omnicloudmask": (
        ".cloudmasking",
        "create_cloud_mask_with_omnicloudmask",
    ),
    "apply_binary_cloud_mask_to_image": (
        ".cloudmasking",
        "apply_binary_cloud_mask_to_image",
    ),
    "radiometric_normalization": (".radiometric_normalization", "radiometric_normalization"),
    "write_seamline_metadata_gpkg": (".seamline_metadata", "write_seamline_metadata_gpkg"),
    "align_image_pair": (".alignment", "align_image_pair"),
    "AlignmentResult": (".alignment", "AlignmentResult"),
}

__all__ = list(_EXPORTS.keys())


def __getattr__(name: str) -> Any:
    """Lazily resolve a public preprocessing export.
    Args:
        name: Export name requested from the preprocessing namespace.
    Returns:
        The resolved export object.
    """
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = _EXPORTS[name]
    value = getattr(import_module(module_name, __name__), attr_name)
    globals()[name] = value
    return value
