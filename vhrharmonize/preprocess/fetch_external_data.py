"""Shared external atmosphere data fetchers."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import date, datetime, timedelta, timezone
import os
from typing import Dict, List, Optional, Tuple

import requests

from vhrharmonize.preprocess.helpers import log


# ---------------------------------------------------------------------------
# DEM
# ---------------------------------------------------------------------------


DEFAULT_OPENTOPOGRAPHY_GLOBALDEM_ENDPOINT = "https://portal.opentopography.org/API/globaldem"
DEFAULT_OPENTOPOGRAPHY_DEMTYPE = "SRTMGL1_Ellip"


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
            "OpenTopography API key is required for dem_file_path='online'. "
            "Set --dem-online-api-key or OPENTOPOGRAPHY_API_KEY."
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
    log(
        f"Downloading {demtype} DEM to {os.path.basename(output_tif_path)}",
        enabled=log_to_console,
        step="dem",
        scene_basename=scene_basename,
    )
    response = requests.get(endpoint, params=params, timeout=timeout_s)
    response.raise_for_status()

    content_type = (response.headers.get("content-type") or "").lower()
    if "application/json" in content_type:
        raise ValueError(f"OpenTopography DEM request returned JSON instead of GeoTIFF: {response.text}")
    if "text/html" in content_type:
        raise ValueError("OpenTopography DEM request returned HTML instead of GeoTIFF.")

    with open(output_tif_path, "wb") as handle:
        handle.write(response.content)
    log(
        f"Wrote DEM {os.path.basename(output_tif_path)}",
        enabled=log_to_console,
        step="dem",
        scene_basename=scene_basename,
    )
    return output_tif_path


# ---------------------------------------------------------------------------
# NASA POWER
# ---------------------------------------------------------------------------


@dataclass
class AtmosphereEstimate:
    """Scene-level atmosphere estimate from NASA APIs."""

    aot550: Optional[float]
    water_vapor: Optional[float]
    ozone_cm_atm: Optional[float]
    source: str
    sample_count: int
    date_used: str


def _linspace(min_v: float, max_v: float, n: int) -> List[float]:
    if n <= 1:
        return [(min_v + max_v) / 2.0]
    step = (max_v - min_v) / float(n - 1)
    return [min_v + (i * step) for i in range(n)]


def _grid_points_from_bbox_wgs84(
    min_lon: float,
    min_lat: float,
    max_lon: float,
    max_lat: float,
    grid_size: int,
) -> List[Tuple[float, float]]:
    if grid_size < 1:
        raise ValueError("grid_size must be >= 1")
    lons = _linspace(min_lon, max_lon, grid_size)
    lats = _linspace(min_lat, max_lat, grid_size)
    return [(lon, lat) for lat in lats for lon in lons]


def _parse_power_value(value) -> Optional[float]:
    if value is None:
        return None
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    if v <= -900:
        return None
    return v


def _fetch_power_daily_point(
    day_utc: date,
    lon: float,
    lat: float,
    *,
    endpoint: str = "https://power.larc.nasa.gov/api/temporal/daily/point",
    timeout_s: float = 30.0,
) -> Dict[str, Optional[float]]:
    ymd = day_utc.strftime("%Y%m%d")
    params = {
        "parameters": "AOD_55,PW,TO3",
        "community": "RE",
        "longitude": f"{lon:.8f}",
        "latitude": f"{lat:.8f}",
        "start": ymd,
        "end": ymd,
        "format": "JSON",
    }
    resp = requests.get(endpoint, params=params, timeout=timeout_s)
    resp.raise_for_status()
    payload = resp.json()
    param = payload.get("properties", {}).get("parameter", {})

    aod = _parse_power_value((param.get("AOD_55", {}) or {}).get(ymd))
    pw_cm = _parse_power_value((param.get("PW", {}) or {}).get(ymd))
    to3_du = _parse_power_value((param.get("TO3", {}) or {}).get(ymd))

    return {
        "aot550": aod,
        "water_vapor": pw_cm,
        "ozone_cm_atm": (to3_du / 1000.0) if to3_du is not None else None,
    }


def fetch_power_atmosphere_for_bbox(
    *,
    day_utc: date,
    min_lon: float,
    min_lat: float,
    max_lon: float,
    max_lat: float,
    grid_size: int = 3,
    search_days: int = 1,
    endpoint: str = "https://power.larc.nasa.gov/api/temporal/daily/point",
    timeout_s: float = 30.0,
    log_to_console: bool = False,
    scene_basename: str | None = None,
) -> AtmosphereEstimate:
    """Fetch AOT550, water vapor, and ozone for a bbox via NASA POWER daily API."""
    if min_lon > max_lon or min_lat > max_lat:
        raise ValueError("Invalid bbox bounds.")
    if grid_size < 1:
        raise ValueError("grid_size must be >= 1")
    if search_days < 0:
        raise ValueError("search_days must be >= 0")

    points = _grid_points_from_bbox_wgs84(min_lon, min_lat, max_lon, max_lat, grid_size)

    offsets = [0]
    for i in range(1, search_days + 1):
        offsets.extend([i, -i])

    for delta in offsets:
        query_day = day_utc + timedelta(days=delta)
        aod_vals: List[float] = []
        wv_vals: List[float] = []
        oz_vals: List[float] = []
        for lon, lat in points:
            item = _fetch_power_daily_point(
                query_day,
                lon,
                lat,
                endpoint=endpoint,
                timeout_s=timeout_s,
            )
            if item["aot550"] is not None:
                aod_vals.append(float(item["aot550"]))
            if item["water_vapor"] is not None:
                wv_vals.append(float(item["water_vapor"]))
            if item["ozone_cm_atm"] is not None:
                oz_vals.append(float(item["ozone_cm_atm"]))

        if aod_vals or wv_vals or oz_vals:
            log(
                f"Using NASA POWER values from {query_day.isoformat()}",
                enabled=log_to_console,
                step="fetch_atmosphere",
                scene_basename=scene_basename,
            )
            return AtmosphereEstimate(
                aot550=(sum(aod_vals) / len(aod_vals)) if aod_vals else None,
                water_vapor=(sum(wv_vals) / len(wv_vals)) if wv_vals else None,
                ozone_cm_atm=(sum(oz_vals) / len(oz_vals)) if oz_vals else None,
                source="nasa_power_daily_point_grid",
                sample_count=len(points),
                date_used=query_day.isoformat(),
            )

    log(
        "No external atmosphere values found",
        enabled=log_to_console,
        step="fetch_atmosphere",
        scene_basename=scene_basename,
    )
    return AtmosphereEstimate(
        aot550=None,
        water_vapor=None,
        ozone_cm_atm=None,
        source="nasa_power_daily_point_grid",
        sample_count=len(points),
        date_used=day_utc.isoformat(),
    )


# ---------------------------------------------------------------------------
# MODIS / GEE
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ModisWaterVaporEstimate:
    """Provider-neutral MODIS atmosphere lookup result."""

    status: str
    selected_collection: Optional[str]
    modis_raw_value: Optional[float]
    modis_scaled_g_cm2: Optional[float]
    water_vapor_preset: Optional[float]
    aod_raw: Optional[float]
    aod_scaled: Optional[float]
    default_visibility: Optional[float]
    modtran_atm: Optional[str]
    image_time_utc: Optional[str]
    abs_time_diff_hours: Optional[float]
    error: Optional[str] = None

    def to_dict(self) -> Dict[str, Optional[float] | Optional[str]]:
        return asdict(self)


def _load_env_file(env_file_path: Optional[str]) -> Dict[str, str]:
    values: Dict[str, str] = {}
    if not env_file_path or not os.path.exists(env_file_path):
        return values
    with open(env_file_path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip().strip('"').strip("'")
            if key:
                values[key] = value
    return values


def _resolve_ee_project(cli_value: Optional[str], env_values: Dict[str, str]) -> Optional[str]:
    if cli_value:
        return cli_value
    for key in ("GEE_PROJECT", "EE_PROJECT", "GOOGLE_EARTH_ENGINE_PROJECT"):
        if env_values.get(key):
            return env_values[key]
    return None


def init_ee_client(
    *,
    ee_project: Optional[str] = None,
    authenticate: bool = False,
    env_file: Optional[str] = None,
):
    """Initialize and return an Earth Engine client."""
    import ee

    env_values = _load_env_file(env_file)
    resolved_project = _resolve_ee_project(ee_project, env_values)
    if authenticate:
        ee.Authenticate()
    if resolved_project:
        ee.Initialize(project=resolved_project)
    else:
        ee.Initialize()
    return ee


def _fetch_collection_value(
    ee,
    collection_id: str,
    band_name: str,
    *,
    min_lon: float,
    min_lat: float,
    max_lon: float,
    max_lat: float,
    scene_dt_utc: datetime,
    hours_window: int,
    reduce_scale_m: int,
):
    start = (scene_dt_utc - timedelta(hours=hours_window)).strftime("%Y-%m-%dT%H:%M:%S")
    end = (scene_dt_utc + timedelta(hours=hours_window)).strftime("%Y-%m-%dT%H:%M:%S")
    ee_geom = ee.Geometry.Rectangle([min_lon, min_lat, max_lon, max_lat])
    scene_ms = int(scene_dt_utc.timestamp() * 1000)

    col = ee.ImageCollection(collection_id).filterDate(start, end).filterBounds(ee_geom)
    size = col.size().getInfo()
    if size == 0:
        return None

    def _add_time_diff(img):
        diff = ee.Number(img.get("system:time_start")).subtract(scene_ms).abs()
        return img.set("abs_time_diff", diff)

    best = col.map(_add_time_diff).sort("abs_time_diff").first()
    best_info = best.getInfo()
    props = best_info.get("properties", {})
    bands = [band.get("id") for band in best_info.get("bands", [])]
    if band_name not in bands:
        return {
            "collection": collection_id,
            "band_found": False,
            "bands": bands,
            "raw_value": None,
            "image_time_utc": None,
            "abs_time_diff_hours": None,
        }

    stats = best.select(band_name).reduceRegion(
        reducer=ee.Reducer.median(),
        geometry=ee_geom,
        scale=reduce_scale_m,
        bestEffort=True,
        maxPixels=1_000_000_000,
    ).getInfo()
    raw_value = stats.get(band_name) if isinstance(stats, dict) else None
    time_start_ms = props.get("system:time_start")
    abs_diff_ms = abs(int(time_start_ms) - scene_ms) if time_start_ms is not None else None
    image_time_utc = (
        datetime.fromtimestamp(time_start_ms / 1000, tz=timezone.utc).isoformat()
        if time_start_ms is not None
        else None
    )
    return {
        "collection": collection_id,
        "band_found": True,
        "bands": bands,
        "raw_value": raw_value,
        "image_time_utc": image_time_utc,
        "abs_time_diff_hours": (abs_diff_ms / 3_600_000) if abs_diff_ms is not None else None,
    }


def _fetch_first_available_band_value(
    ee,
    collection_id: str,
    band_candidates: List[str],
    *,
    min_lon: float,
    min_lat: float,
    max_lon: float,
    max_lat: float,
    scene_dt_utc: datetime,
    hours_window: int,
    reduce_scale_m: int,
):
    start = (scene_dt_utc - timedelta(hours=hours_window)).strftime("%Y-%m-%dT%H:%M:%S")
    end = (scene_dt_utc + timedelta(hours=hours_window)).strftime("%Y-%m-%dT%H:%M:%S")
    ee_geom = ee.Geometry.Rectangle([min_lon, min_lat, max_lon, max_lat])
    scene_ms = int(scene_dt_utc.timestamp() * 1000)

    col = ee.ImageCollection(collection_id).filterDate(start, end).filterBounds(ee_geom)
    size = col.size().getInfo()
    if size == 0:
        return None

    def _add_time_diff(img):
        diff = ee.Number(img.get("system:time_start")).subtract(scene_ms).abs()
        return img.set("abs_time_diff", diff)

    best = col.map(_add_time_diff).sort("abs_time_diff").first()
    best_info = best.getInfo()
    props = best_info.get("properties", {})
    bands = [band.get("id") for band in best_info.get("bands", [])]

    selected_band = None
    for band in band_candidates:
        if band in bands:
            selected_band = band
            break
    if selected_band is None:
        return {
            "collection": collection_id,
            "band_found": False,
            "selected_band": None,
            "bands": bands,
            "raw_value": None,
            "image_time_utc": None,
            "abs_time_diff_hours": None,
        }

    stats = best.select(selected_band).reduceRegion(
        reducer=ee.Reducer.median(),
        geometry=ee_geom,
        scale=reduce_scale_m,
        bestEffort=True,
        maxPixels=1_000_000_000,
    ).getInfo()
    raw_value = stats.get(selected_band) if isinstance(stats, dict) else None
    time_start_ms = props.get("system:time_start")
    abs_diff_ms = abs(int(time_start_ms) - scene_ms) if time_start_ms is not None else None
    image_time_utc = (
        datetime.fromtimestamp(time_start_ms / 1000, tz=timezone.utc).isoformat()
        if time_start_ms is not None
        else None
    )
    return {
        "collection": collection_id,
        "band_found": True,
        "selected_band": selected_band,
        "bands": bands,
        "raw_value": raw_value,
        "image_time_utc": image_time_utc,
        "abs_time_diff_hours": (abs_diff_ms / 3_600_000) if abs_diff_ms is not None else None,
    }


def _choose_best_candidate(candidates: List[Dict]) -> Optional[Dict]:
    valid = [
        candidate
        for candidate in candidates
        if candidate
        and candidate.get("band_found")
        and candidate.get("raw_value") is not None
        and candidate.get("abs_time_diff_hours") is not None
    ]
    if not valid:
        return None
    return sorted(valid, key=lambda item: item["abs_time_diff_hours"])[0]


def _choose_modtran_atm_by_lat_month(lat_deg: float, month: int) -> str:
    if abs(lat_deg) <= 23.5:
        return "Tropical Atmosphere"
    if month in (12, 1, 2):
        return "Mid-Latitude Winter"
    return "Mid-Latitude Summer"


def _visibility_from_aod(aod: float) -> float:
    if aod < 0.10:
        return 30.0
    if aod < 0.20:
        return 20.0
    if aod < 0.35:
        return 12.0
    if aod < 0.60:
        return 8.0
    return 5.0


def fetch_modis_water_vapor_for_bbox(
    *,
    scene_datetime_utc: datetime,
    min_lon: float,
    min_lat: float,
    max_lon: float,
    max_lat: float,
    ee=None,
    ee_project: Optional[str] = None,
    authenticate: bool = False,
    env_file: Optional[str] = None,
    hours_window: int = 24,
    terra_collection: str = "MODIS/061/MOD08_D3",
    aqua_collection: str = "MODIS/061/MYD08_D3",
    band_name: str = "Atmospheric_Water_Vapor_Mean",
    aod_band_candidates: Optional[List[str]] = None,
    aod_scale_factor: float = 0.001,
    modis_scale_factor: float = 0.001,
    modtran_baseline_water_vapor: float = 2.92,
    reduce_scale_m: int = 1000,
    log_to_console: bool = False,
    scene_basename: str | None = None,
) -> ModisWaterVaporEstimate:
    """Fetch MODIS water vapor and companion AOD-derived FLAASH inputs for a bbox."""
    log(
        "Querying MODIS atmosphere data",
        enabled=log_to_console,
        step="fetch_atmosphere",
        scene_basename=scene_basename,
    )
    if ee is None:
        ee = init_ee_client(ee_project=ee_project, authenticate=authenticate, env_file=env_file)
    if scene_datetime_utc.tzinfo is None:
        scene_datetime_utc = scene_datetime_utc.replace(tzinfo=timezone.utc)
    else:
        scene_datetime_utc = scene_datetime_utc.astimezone(timezone.utc)

    aod_band_candidates = aod_band_candidates or [
        "Aerosol_Optical_Depth_Land_Ocean_Mean",
        "Aerosol_Optical_Depth_Land_Ocean_Mean_Mean",
        "Aerosol_Optical_Depth_Land_Ocean",
    ]

    centroid_lat = (min_lat + max_lat) / 2.0
    terra = _fetch_collection_value(
        ee,
        terra_collection,
        band_name,
        min_lon=min_lon,
        min_lat=min_lat,
        max_lon=max_lon,
        max_lat=max_lat,
        scene_dt_utc=scene_datetime_utc,
        hours_window=hours_window,
        reduce_scale_m=reduce_scale_m,
    )
    aqua = _fetch_collection_value(
        ee,
        aqua_collection,
        band_name,
        min_lon=min_lon,
        min_lat=min_lat,
        max_lon=max_lon,
        max_lat=max_lat,
        scene_dt_utc=scene_datetime_utc,
        hours_window=hours_window,
        reduce_scale_m=reduce_scale_m,
    )
    best = _choose_best_candidate([terra, aqua])
    if best is None:
        log(
            "No MODIS match found",
            enabled=log_to_console,
            step="fetch_atmosphere",
            scene_basename=scene_basename,
        )
        return ModisWaterVaporEstimate(
            status="no_modis_match",
            selected_collection=None,
            modis_raw_value=None,
            modis_scaled_g_cm2=None,
            water_vapor_preset=None,
            aod_raw=None,
            aod_scaled=None,
            default_visibility=None,
            modtran_atm=None,
            image_time_utc=None,
            abs_time_diff_hours=None,
        )

    raw_value = float(best["raw_value"])
    scaled = raw_value * modis_scale_factor
    preset = scaled / modtran_baseline_water_vapor

    terra_aod = _fetch_first_available_band_value(
        ee,
        terra_collection,
        aod_band_candidates,
        min_lon=min_lon,
        min_lat=min_lat,
        max_lon=max_lon,
        max_lat=max_lat,
        scene_dt_utc=scene_datetime_utc,
        hours_window=hours_window,
        reduce_scale_m=reduce_scale_m,
    )
    aqua_aod = _fetch_first_available_band_value(
        ee,
        aqua_collection,
        aod_band_candidates,
        min_lon=min_lon,
        min_lat=min_lat,
        max_lon=max_lon,
        max_lat=max_lat,
        scene_dt_utc=scene_datetime_utc,
        hours_window=hours_window,
        reduce_scale_m=reduce_scale_m,
    )
    best_aod = _choose_best_candidate([terra_aod, aqua_aod])
    aod_raw = float(best_aod["raw_value"]) if best_aod and best_aod.get("raw_value") is not None else None
    aod_scaled = aod_raw * aod_scale_factor if aod_raw is not None else None
    default_visibility = _visibility_from_aod(aod_scaled) if aod_scaled is not None else None
    modtran_atm = _choose_modtran_atm_by_lat_month(centroid_lat, scene_datetime_utc.month)
    result = ModisWaterVaporEstimate(
        status="ok",
        selected_collection=best["collection"],
        modis_raw_value=raw_value,
        modis_scaled_g_cm2=scaled,
        water_vapor_preset=preset,
        aod_raw=aod_raw,
        aod_scaled=aod_scaled,
        default_visibility=default_visibility,
        modtran_atm=modtran_atm,
        image_time_utc=best["image_time_utc"],
        abs_time_diff_hours=best["abs_time_diff_hours"],
    )
    log(
        "Fetched MODIS atmosphere values",
        enabled=log_to_console,
        step="fetch_atmosphere",
        scene_basename=scene_basename,
    )
    return result


__all__ = [
    "AtmosphereEstimate",
    "DEFAULT_OPENTOPOGRAPHY_DEMTYPE",
    "DEFAULT_OPENTOPOGRAPHY_GLOBALDEM_ENDPOINT",
    "ModisWaterVaporEstimate",
    "download_opentopography_dem_for_bbox",
    "fetch_modis_water_vapor_for_bbox",
    "fetch_power_atmosphere_for_bbox",
    "init_ee_client",
]
