"""Preprocessing step modules."""

from importlib import import_module

_EXPORTS = {
    "atmospheric_correction": (".atmospheric_correction", "atmospheric_correction"),
    "gcp_refined_rpc_orthorectification": (
        ".orthorectification",
        "gcp_refined_rpc_orthorectification",
    ),
    "pansharpen_image": (".pansharpening", "pansharpen_image"),
    "create_cloud_mask_with_omnicloudmask": (
        ".cloudmasking",
        "create_cloud_mask_with_omnicloudmask",
    ),
    "apply_binary_cloud_mask_to_image": (
        ".cloudmasking",
        "apply_binary_cloud_mask_to_image",
    ),
    "apply_rrn": (".radiometric_normalization", "apply_rrn"),
    "apply_arn": (".radiometric_normalization", "apply_arn"),
    "tile_image": (".tiling", "tile_image"),
    "estimate_elastix_transform": (".registration", "estimate_elastix_transform"),
    "apply_elastix_transform": (".registration", "apply_elastix_transform"),
    "run_elastix_registration": (".registration", "run_elastix_registration"),
}

__all__ = list(_EXPORTS.keys())


def __getattr__(name):
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = _EXPORTS[name]
    value = getattr(import_module(module_name, __name__), attr_name)
    globals()[name] = value
    return value
