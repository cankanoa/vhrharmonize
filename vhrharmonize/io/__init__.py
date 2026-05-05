"""Geospatial IO and utility helpers."""

from importlib import import_module
from typing import Any

_EXPORTS = {
    "calculate_raster_overviews": (".geospatial", "calculate_raster_overviews"),
    "shp_to_gpkg": (".geospatial", "shp_to_gpkg"),
    "get_image_percentile_value": (".geospatial", "get_image_percentile_value"),
    "StepOutputPlan": (".workflow_utils", "StepOutputPlan"),
    "build_output_path_from_input": (".workflow_utils", "build_output_path_from_input"),
    "plan_step_outputs": (".workflow_utils", "plan_step_outputs"),
    "resolve_output_dir": (".workflow_utils", "resolve_output_dir"),
    "resolve_temp_dir": (".workflow_utils", "resolve_temp_dir"),
    "resolve_relative_to_input": (".workflow_utils", "resolve_relative_to_input"),
}

__all__ = list(_EXPORTS.keys())


def __getattr__(name: str) -> Any:
    """Lazily resolve a public IO export.
    Args:
        name: Export name requested from the IO namespace.
    Returns:
        The resolved export object.
    """
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = _EXPORTS[name]
    value = getattr(import_module(module_name, __name__), attr_name)
    globals()[name] = value
    return value
