"""Download elevation data: Python implementation and workflow plugin."""

from __future__ import annotations

import math
import os
from typing import Optional

import rasterio
import requests
from rasterio.transform import GCPTransformer, RPCTransformer
from rasterio.warp import transform, transform_bounds

from .base import FunctionPlugin
from vhrharmonize.io.logging import _log, _logged_operation


DEFAULT_OPENTOPOGRAPHY_GLOBALDEM_ENDPOINT = "https://portal.opentopography.org/API/globaldem"
DEFAULT_OPENTOPOGRAPHY_DEMTYPE = "SRTMGL1_Ellip"


def _image_bounds(image_path: str) -> tuple[float, float, float, float]:
    """Read WGS84 coverage from affine, RPC, or GCP georeferencing."""
    with rasterio.open(image_path) as source:
        if source.crs and not source.transform.is_identity:
            return transform_bounds(source.crs, "EPSG:4326", *source.bounds, densify_pts=21)
        gcps, gcp_crs = source.gcps
        rpcs = source.rpcs
        if rpcs:
            transformer = RPCTransformer(rpcs, RPC_HEIGHT=rpcs.height_off)
            crs = "EPSG:4326"
        elif gcps and gcp_crs:
            transformer = GCPTransformer(gcps)
            crs = gcp_crs
        elif source.crs:
            return transform_bounds(source.crs, "EPSG:4326", *source.bounds, densify_pts=21)
        else:
            raise ValueError(
                f"Cannot determine image bounds without CRS, RPCs or GCPs: {image_path}"
            )
        # Sample all edges, not just corners, for nonlinear sensor models.
        rows, cols = [], []
        for index in range(23):
            fraction = index / 22
            rows.extend([0, source.height, fraction * source.height, fraction * source.height])
            cols.extend([fraction * source.width, fraction * source.width, 0, source.width])
        with transformer:
            xs, ys = transformer.xy(rows, cols, offset="ul")
        lons, lats = transform(crs, "EPSG:4326", xs, ys)
        return min(lons), min(lats), max(lons), max(lats)


@_logged_operation("fetch_dem", inputs=("image_paths",), outputs=("output_tif_path",))
def download_opentopography_dem_for_bbox(
    *,
    min_lon: float | None = None,
    min_lat: float | None = None,
    max_lon: float | None = None,
    max_lat: float | None = None,
    image_paths: str | list[str] | None = None,
    output_tif_path: str,
    api_key: Optional[str] = None,
    demtype: str = DEFAULT_OPENTOPOGRAPHY_DEMTYPE,
    endpoint: str = DEFAULT_OPENTOPOGRAPHY_GLOBALDEM_ENDPOINT,
    timeout_s: float = 120.0,
    log_to_console: bool = False,
    scene_basename: str | None = None,
) -> str:
    """Download one DEM covering an explicit WGS84 bbox or the union of images.

    Supply image_paths as one raster path or a nonempty list; each image may use
    a different CRS. Alternatively supply all four bbox coordinates. The two
    modes are mutually exclusive. RPC footprints use the model's reference
    height, so their extent is approximate before a DEM is available.
    """

    bounds = (min_lon, min_lat, max_lon, max_lat)
    if image_paths is not None:
        if any(value is not None for value in bounds):
            raise ValueError("Supply image_paths or bbox coordinates, not both.")
        paths = (
            [image_paths]
            if isinstance(image_paths, (str, os.PathLike))
            else list(image_paths)
        )
        if not paths:
            raise ValueError("image_paths must contain at least one image.")
        image_bounds = [_image_bounds(path) for path in paths]
        if any(
            not all(math.isfinite(v) for v in b) or b[0] >= b[2] or b[1] >= b[3]
            for b in image_bounds
        ):
            raise ValueError(
                "Images must have finite, nonempty bounds that do not cross the antimeridian."
            )
        min_lon = min(b[0] for b in image_bounds)
        min_lat = min(b[1] for b in image_bounds)
        max_lon = max(b[2] for b in image_bounds)
        max_lat = max(b[3] for b in image_bounds)
    elif any(value is None for value in bounds):
        raise ValueError("Supply image_paths or all four bbox coordinates.")
    if not all(math.isfinite(v) for v in (min_lon, min_lat, max_lon, max_lat)) or not (
        -180 <= min_lon < max_lon <= 180 and -90 <= min_lat < max_lat <= 90
    ):
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
    input_dependency_paths = frozenset({"image_paths"})
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
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
