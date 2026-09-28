"""Download elevation data: Python implementation and workflow plugin."""

from __future__ import annotations

import os

from typing import Optional

import requests

from .base import FunctionPlugin
from vhrharmonize.io.logging import _log, _logged_operation


DEFAULT_OPENTOPOGRAPHY_GLOBALDEM_ENDPOINT = "https://portal.opentopography.org/API/globaldem"
DEFAULT_OPENTOPOGRAPHY_DEMTYPE = "SRTMGL1_Ellip"


@_logged_operation("fetch_dem", inputs=(), outputs=("output_tif_path",))
def download_opentopography_dem_for_bbox(
    *,
    min_lon: float,
    min_lat: float,
    max_lon: float,
    max_lat: float,
    output_tif_path: str,
    api_key: Optional[str] = None,
    demtype: str = DEFAULT_OPENTOPOGRAPHY_DEMTYPE,
    endpoint: str = DEFAULT_OPENTOPOGRAPHY_GLOBALDEM_ENDPOINT,
    timeout_s: float = 120.0,
    log_to_console: bool = False,
    scene_basename: str | None = None,
) -> str:
    """Download an OpenTopography DEM subset for a WGS84 bbox."""

    if min_lon > max_lon or min_lat > max_lat:
        raise ValueError("Invalid bbox bounds.")
    api_key = api_key or os.environ.get("OPENTOPOGRAPHY_API_KEY") or os.environ.get("OT_API_KEY")
    if not api_key:
        raise ValueError(
            "OpenTopography API key is required. Set api_key or OPENTOPOGRAPHY_API_KEY."
        )

    os.makedirs(os.path.dirname(output_tif_path) or ".", exist_ok=True)
    params = {
        "demtype": demtype,
        "south": f"{min_lat:.8f}",
        "north": f"{max_lat:.8f}",
        "west": f"{min_lon:.8f}",
        "east": f"{max_lon:.8f}",
        "outputFormat": "GTiff",
        "API_Key": api_key,
    }
    _log(
        f"Downloading {demtype} DEM to {os.path.basename(output_tif_path)}",
        enabled=log_to_console,
        step="dem",
        scene_basename=scene_basename,
    )
    response = requests.get(endpoint, params=params, timeout=timeout_s)
    response.raise_for_status()

    content_type = (response.headers.get("content-type") or "").lower()
    if "application/json" in content_type:
        raise ValueError(
            f"OpenTopography DEM request returned JSON instead of GeoTIFF: {response.text}"
        )
    if "text/html" in content_type:
        raise ValueError("OpenTopography DEM request returned HTML instead of GeoTIFF.")

    with open(output_tif_path, "wb") as handle:
        handle.write(response.content)
    _log(
        f"Wrote DEM {os.path.basename(output_tif_path)}",
        enabled=log_to_console,
        step="dem",
        scene_basename=scene_basename,
    )
    return output_tif_path


class FetchDEM(FunctionPlugin):
    output_path_resolution_paths = frozenset({"output_tif_path"})
    output_dependency_paths = output_path_resolution_paths
    output_target_paths = output_path_resolution_paths
    output_parent_creation_paths = output_path_resolution_paths
    output_collision_check_paths = output_path_resolution_paths
    output_reuse_paths = output_path_resolution_paths
    output_validation_paths = output_path_resolution_paths
    output_invalid_removal_paths = output_path_resolution_paths
    output_overview_calculation_paths = frozenset({"output_tif_path"})
    output_temporary_cleanup_paths = output_path_resolution_paths
    output_hpc_staging_paths = output_path_resolution_paths
    output_hpc_download_paths = output_path_resolution_paths
    output_context_checkpoint_paths = frozenset({"output_tif_path"})
    target = "vhrharmonize.plugins.fetch_dem:download_opentopography_dem_for_bbox"


__all__ = [
    "FetchDEM",
    "download_opentopography_dem_for_bbox",
    "DEFAULT_OPENTOPOGRAPHY_DEMTYPE",
    "DEFAULT_OPENTOPOGRAPHY_GLOBALDEM_ENDPOINT",
]
