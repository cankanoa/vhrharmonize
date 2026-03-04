"""Geospatial IO and utility helpers."""

from importlib import import_module

_EXPORTS = {
    "convert_dn_to_radiance": (".geospatial", "convert_dn_to_radiance"),
    "change_last_folder_and_add_suffix": (".geospatial", "change_last_folder_and_add_suffix"),
    "wsl_to_windows_path": (".geospatial", "wsl_to_windows_path"),
    "windows_to_wsl_path": (".geospatial", "windows_to_wsl_path"),
    "shp_to_gpkg": (".geospatial", "shp_to_gpkg"),
    "convert_dat_to_tif": (".geospatial", "convert_dat_to_tif"),
    "copy_tif_file": (".geospatial", "copy_tif_file"),
    "spectral_scale_image": (".geospatial", "spectral_scale_image"),
    "stretch_spectral_values": (".geospatial", "stretch_spectral_values"),
    "get_image_largest_value": (".geospatial", "get_image_largest_value"),
    "get_image_percentile_value": (".geospatial", "get_image_percentile_value"),
    "replace_band_continuous_values_in_largest_segment": (
        ".geospatial",
        "replace_band_continuous_values_in_largest_segment",
    ),
    "scale_gcps_geojson": (".geospatial", "scale_gcps_geojson"),
    "translate_gcp_image_to_origin": (".geospatial", "translate_gcp_image_to_origin"),
}

__all__ = list(_EXPORTS.keys())


def __getattr__(name):
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = _EXPORTS[name]
    value = getattr(import_module(module_name, __name__), attr_name)
    globals()[name] = value
    return value
