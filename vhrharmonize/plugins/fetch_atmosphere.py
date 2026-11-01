"""Shared external atmosphere data fetchers."""

from __future__ import annotations

import os
import json
import math
from pathlib import Path
from dataclasses import asdict, dataclass, field
from datetime import date, datetime, timedelta, timezone
from typing import Any, Dict, List, Optional, Tuple

import requests

from .base import FunctionPlugin
from vhrharmonize.io.logging import _log, _logged_operation
from vhrharmonize.io.progress import progress


_ATMOSPHERE_PARAMETERS = ("aot550", "water_vapor", "ozone_cm_atm")


class _AtmosphereDiagnostics:
    """Keep the same warning records in console logs and retrieval JSON."""

    def __init__(self, source, log_to_console, scene_basename):
        self.source = source
        self.log_to_console = log_to_console
        self.scene_basename = scene_basename
        self.warnings = []
        self.lookup_attempts = []

    def warn(self, code, message, **details):
        self.warnings.append(dict(code=code, message=message, source=self.source, **details))
        _log("WARNING | " + message, enabled=self.log_to_console,
             step="fetch_atmosphere", scene_basename=self.scene_basename)

    def metadata(self):
        return dict(warnings=self.warnings, lookup_attempts=self.lookup_attempts)

    def call(self, function, *args, lookup, **kwargs):
        try:
            return function(*args, **kwargs)
        except Exception as error:
            description = f"{type(error).__name__}: {error}"
            self.lookup_attempts.append(dict(**lookup, status="error", error=description))
            self.warn("lookup_error", f"{self.source} lookup failed: {description}. No automatic source switch or default substitution was performed.", **lookup)
            # Preserve the original exception and its behavior; the adapter saves
            # these diagnostics before re-raising it.
            error.atmosphere_diagnostics = self.metadata()
            raise


def _validate_fallback_values(values):
    values = {} if values is None else dict(values)
    unknown = set(values) - set(_ATMOSPHERE_PARAMETERS)
    if unknown:
        raise ValueError(f"Unsupported atmospheric fallback parameters: {sorted(unknown)}")
    result = {}
    for name, value in values.items():
        if isinstance(value, bool):
            raise ValueError(f"fallback_values.{name} must be a finite nonnegative number.")
        value = float(value)
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"fallback_values.{name} must be a finite nonnegative number.")
        result[name] = value
    return result


def _apply_atmosphere_defaults(values, fallbacks, diagnostics):
    """Substitute only explicitly configured defaults, retaining an audit trail."""
    values = dict(values)
    used = {}
    units = {"aot550": "", "water_vapor": " g/cm^2", "ozone_cm_atm": " cm-atm"}
    for name in _ATMOSPHERE_PARAMETERS:
        if values[name] is not None:
            continue
        if name in fallbacks:
            values[name] = used[name] = fallbacks[name]
            diagnostics.warn("default_used",
                f"{diagnostics.source} supplied no {name}; using configured default {fallbacks[name]}{units[name]}.",
                parameter=name, value=fallbacks[name])
        else:
            diagnostics.warn("parameter_unavailable",
                f"{diagnostics.source} supplied no {name}; no fallback is configured and the value remains null.",
                parameter=name)
    return values, used


@dataclass
class AtmosphereEstimate:
    """Scene-level atmosphere estimate from NASA APIs."""

    aot550: Optional[float]
    water_vapor: Optional[float]
    ozone_cm_atm: Optional[float]
    source: str
    sample_count: int
    date_used: str
    date_requested: Optional[str] = None
    warnings: List[Dict[str, Any]] = field(default_factory=list)
    lookup_attempts: List[Dict[str, Any]] = field(default_factory=list)
    defaults_used: Dict[str, float] = field(default_factory=dict)


def _linspace(min_v: float, max_v: float, n: int) -> List[float]:
    """Return evenly spaced values.
    Args:
        min_v: Minimum value.
        max_v: Maximum value.
        n: Number of values to generate.
    Returns:
        Evenly spaced values across the range.
    """
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
    """Build a regular point grid for a WGS84 bbox.
    Args:
        min_lon: Minimum longitude.
        min_lat: Minimum latitude.
        max_lon: Maximum longitude.
        max_lat: Maximum latitude.
        grid_size: Number of samples per side.
    Returns:
        Grid point longitude and latitude pairs.
    """
    if grid_size < 1:
        raise ValueError("grid_size must be >= 1")
    lons = _linspace(min_lon, max_lon, grid_size)
    lats = _linspace(min_lat, max_lat, grid_size)
    return [(lon, lat) for lat in lats for lon in lons]


def _parse_power_value(value: object) -> Optional[float]:
    """Parse a NASA POWER scalar value.
    Args:
        value: Raw POWER value.
    Returns:
        Parsed float value or None.
    """
    if value is None:
        return None
    v = float(value)
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
    """Fetch POWER atmosphere values for one point and day.
    Args:
        day_utc: Query day in UTC.
        lon: Longitude in WGS84.
        lat: Latitude in WGS84.
        endpoint: NASA POWER endpoint.
        timeout_s: Request timeout in seconds.
    Returns:
        POWER atmosphere values for the point.
    """
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


@_logged_operation("fetch_atmosphere", inputs=(), outputs=(), worker_progress=True)
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
    fallback_values: dict[str, float] | None = None,
    log_to_console: bool = False,
    scene_basename: str | None = None,
) -> AtmosphereEstimate:
    """Fetch AOT550, water vapor, and ozone with logged/persisted diagnostics.

    Args:
        day_utc: Requested acquisition date in UTC.
        min_lon: Bounding-box west longitude.
        min_lat: Bounding-box south latitude.
        max_lon: Bounding-box east longitude.
        max_lat: Bounding-box north latitude.
        grid_size: Number of sample points per bounding-box side.
        search_days: Search this many days before and after acquisition.
        endpoint: NASA POWER daily point API URL.
        timeout_s: Timeout per point request in seconds.
        fallback_values: Explicit defaults for missing aot550, water_vapor (g/cm^2), and ozone_cm_atm (cm-atm). None leaves missing values null.
        log_to_console: Print lookup and fallback warnings.
        scene_basename: Optional scene label in console logs.
    """
    if min_lon > max_lon or min_lat > max_lat:
        raise ValueError("Invalid bbox bounds.")
    if grid_size < 1:
        raise ValueError("grid_size must be >= 1")
    if search_days < 0:
        raise ValueError("search_days must be >= 0")
    fallbacks = _validate_fallback_values(fallback_values)
    diagnostics = _AtmosphereDiagnostics("nasa_power", log_to_console, scene_basename)
    points = _grid_points_from_bbox_wgs84(min_lon, min_lat, max_lon, max_lat, grid_size)
    offsets = [0]
    for i in range(1, search_days + 1):
        offsets.extend([i, -i])

    values = dict.fromkeys(_ATMOSPHERE_PARAMETERS)
    date_used = day_utc
    for delta in offsets:
        query_day = day_utc + timedelta(days=delta)
        if delta:
            diagnostics.warn("alternate_date_lookup",
                f"No atmospheric values found on earlier NASA POWER lookup dates; trying {query_day.isoformat()} instead of acquisition date {day_utc.isoformat()}.",
                requested_date=day_utc.isoformat(), query_date=query_day.isoformat(), offset_days=delta)
        samples = {name: [] for name in _ATMOSPHERE_PARAMETERS}
        for lon, lat in progress(points, desc=f"Fetching atmosphere for {query_day}", unit="points"):
            lookup = dict(query_date=query_day.isoformat(), longitude=lon, latitude=lat)
            item = diagnostics.call(_fetch_power_daily_point, query_day, lon, lat,
                endpoint=endpoint, timeout_s=timeout_s, lookup=lookup)
            missing = [name for name in _ATMOSPHERE_PARAMETERS if item[name] is None]
            diagnostics.lookup_attempts.append(dict(**lookup, values=dict(item), missing_parameters=missing,
                status="missing" if len(missing) == 3 else "partial" if missing else "ok"))
            if missing:
                diagnostics.warn("lookup_missing_values",
                    f"NASA POWER lookup on {query_day.isoformat()} at ({lon:g}, {lat:g}) returned no {', '.join(missing)}.",
                    **lookup, parameters=missing)
            for name in _ATMOSPHERE_PARAMETERS:
                if item[name] is not None:
                    samples[name].append(float(item[name]))
        if any(samples.values()):
            values = {name: sum(items) / len(items) if items else None for name, items in samples.items()}
            date_used = query_day
            if delta:
                diagnostics.warn("alternate_date_used",
                    f"Using NASA POWER values from {query_day.isoformat()} instead of acquisition date {day_utc.isoformat()}.",
                    requested_date=day_utc.isoformat(), date_used=query_day.isoformat(),
                    parameters=[name for name, value in values.items() if value is not None])
            _log(f"Using NASA POWER values from {query_day.isoformat()}", enabled=log_to_console,
                 step="fetch_atmosphere", scene_basename=scene_basename)
            break
    else:
        diagnostics.warn("no_atmosphere_match", "No external atmosphere values found from NASA POWER on any searched date.",
            dates_attempted=[(day_utc + timedelta(days=delta)).isoformat() for delta in offsets])
    values, used = _apply_atmosphere_defaults(values, fallbacks, diagnostics)
    return AtmosphereEstimate(**values, source="nasa_power_daily_point_grid", sample_count=len(points),
        date_used=date_used.isoformat(), date_requested=day_utc.isoformat(), defaults_used=used,
        **diagnostics.metadata())


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
    aot550: Optional[float] = None
    water_vapor: Optional[float] = None
    ozone_cm_atm: Optional[float] = None
    scene_datetime_utc: Optional[str] = None
    warnings: List[Dict[str, Any]] = field(default_factory=list)
    lookup_attempts: List[Dict[str, Any]] = field(default_factory=list)
    defaults_used: Dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert the MODIS estimate to a dictionary.
        Args:
            self: MODIS estimate instance.
        Returns:
            Dictionary representation of the estimate.
        """
        return asdict(self)


def _load_env_file(env_file_path: Optional[str]) -> Dict[str, str]:
    """Load key-value pairs from an env file.
    Args:
        env_file_path: Optional env file path.
    Returns:
        Loaded environment-style key-value pairs.
    """
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
    """Resolve the Earth Engine project id.
    Args:
        cli_value: CLI-provided Earth Engine project id.
        env_values: Loaded env-file values.
    Returns:
        Resolved Earth Engine project id or None.
    """
    if cli_value:
        return cli_value
    for key in ("GEE_PROJECT", "EE_PROJECT", "GOOGLE_EARTH_ENGINE_PROJECT"):
        if env_values.get(key):
            return env_values[key]
    return None


def _init_ee_client(
    *,
    ee_project: Optional[str] = None,
    authenticate: bool = False,
    env_file: Optional[str] = None,
) -> Any:
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
    ee: Any,
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
) -> Optional[Dict[str, object]]:
    """Fetch a MODIS collection value near a scene time.
    Args:
        ee: Initialized Earth Engine module or client.
        collection_id: Earth Engine collection id.
        band_name: Band name to query.
        min_lon: Minimum longitude.
        min_lat: Minimum latitude.
        max_lon: Maximum longitude.
        max_lat: Maximum latitude.
        scene_dt_utc: Scene datetime in UTC.
        hours_window: Search window in hours.
        reduce_scale_m: Reduce-region scale in meters.
    Returns:
        Best collection match summary or None.
    """
    start = (scene_dt_utc - timedelta(hours=hours_window)).strftime("%Y-%m-%dT%H:%M:%S")
    end = (scene_dt_utc + timedelta(hours=hours_window)).strftime("%Y-%m-%dT%H:%M:%S")
    ee_geom = ee.Geometry.Rectangle([min_lon, min_lat, max_lon, max_lat])
    scene_ms = int(scene_dt_utc.timestamp() * 1000)

    col = ee.ImageCollection(collection_id).filterDate(start, end).filterBounds(ee_geom)
    size = col.size().getInfo()
    if size == 0:
        return None

    def _add_time_diff(img: Any) -> Any:
        """Attach absolute time difference to an image.
        Args:
            img: Earth Engine image.
        Returns:
            Image with an abs_time_diff property.
        """
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

    stats = (
        best.select(band_name)
        .reduceRegion(
            reducer=ee.Reducer.median(),
            geometry=ee_geom,
            scale=reduce_scale_m,
            bestEffort=True,
            maxPixels=1_000_000_000,
        )
        .getInfo()
    )
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
    ee: Any,
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
) -> Optional[Dict[str, object]]:
    """Fetch the first available candidate band near a scene time.
    Args:
        ee: Initialized Earth Engine module or client.
        collection_id: Earth Engine collection id.
        band_candidates: Candidate band names in priority order.
        min_lon: Minimum longitude.
        min_lat: Minimum latitude.
        max_lon: Maximum longitude.
        max_lat: Maximum latitude.
        scene_dt_utc: Scene datetime in UTC.
        hours_window: Search window in hours.
        reduce_scale_m: Reduce-region scale in meters.
    Returns:
        Best candidate-band summary or None.
    """
    start = (scene_dt_utc - timedelta(hours=hours_window)).strftime("%Y-%m-%dT%H:%M:%S")
    end = (scene_dt_utc + timedelta(hours=hours_window)).strftime("%Y-%m-%dT%H:%M:%S")
    ee_geom = ee.Geometry.Rectangle([min_lon, min_lat, max_lon, max_lat])
    scene_ms = int(scene_dt_utc.timestamp() * 1000)

    col = ee.ImageCollection(collection_id).filterDate(start, end).filterBounds(ee_geom)
    size = col.size().getInfo()
    if size == 0:
        return None

    def _add_time_diff(img: Any) -> Any:
        """Attach absolute time difference to an image.
        Args:
            img: Earth Engine image.
        Returns:
            Image with an abs_time_diff property.
        """
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

    stats = (
        best.select(selected_band)
        .reduceRegion(
            reducer=ee.Reducer.median(),
            geometry=ee_geom,
            scale=reduce_scale_m,
            bestEffort=True,
            maxPixels=1_000_000_000,
        )
        .getInfo()
    )
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
    """Choose the closest valid MODIS candidate.
    Args:
        candidates: Candidate result dictionaries.
    Returns:
        Best candidate or None.
    """
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
    """Choose a MODTRAN atmosphere profile.
    Args:
        lat_deg: Scene centroid latitude in degrees.
        month: Scene month number.
    Returns:
        Suggested MODTRAN atmosphere profile name.
    """
    if abs(lat_deg) <= 23.5:
        return "Tropical Atmosphere"
    if month in (12, 1, 2):
        return "Mid-Latitude Winter"
    return "Mid-Latitude Summer"


def _visibility_from_aod(aod: float) -> float:
    """Estimate visibility from AOD.
    Args:
        aod: Aerosol optical depth.
    Returns:
        Approximate visibility in kilometers.
    """
    if aod < 0.10:
        return 30.0
    if aod < 0.20:
        return 20.0
    if aod < 0.35:
        return 12.0
    if aod < 0.60:
        return 8.0
    return 5.0


def _record_modis_lookup(diagnostics, candidate, details):
    """Describe missing images, missing bands, and unusable regional samples."""
    requested = details["bands_requested"]
    if candidate is None:
        diagnostics.lookup_attempts.append(dict(**details, status="no_images"))
        diagnostics.warn("lookup_no_images",
            f"MODIS lookup found no images in {details['collection']} for {details['parameter']} within +/-{details['hours_window']} hours of {details['scene_datetime_utc']}.", **details)
        return
    selected = candidate.get("selected_band") or (requested[0] if candidate.get("band_found") else None)
    if "bands" in candidate:
        for band in requested:
            if band in candidate["bands"]:
                break
            diagnostics.warn("band_not_found", f"MODIS band {band} not found in {details['collection']}.",
                **details, band=band)
    if not candidate.get("band_found"):
        status = "band_not_found"
    elif candidate.get("raw_value") is None:
        status = "no_valid_pixels"
    elif candidate.get("abs_time_diff_hours") is None:
        status = "missing_timestamp"
    else:
        status = "ok"
    diagnostics.lookup_attempts.append(dict(**details, status=status, selected_band=selected,
        raw_value=candidate.get("raw_value"), image_time_utc=candidate.get("image_time_utc"),
        abs_time_diff_hours=candidate.get("abs_time_diff_hours")))
    if status != "ok":
        diagnostics.warn("lookup_missing_values",
            f"MODIS lookup in {details['collection']} supplied no usable {details['parameter']} ({status}).", **details, reason=status)
    if status == "ok" and selected in requested and requested.index(selected) > 0:
        diagnostics.warn("alternate_band_used",
            f"Using alternate MODIS band {selected} in {details['collection']} for {details['parameter']}.",
            **details, selected_band=selected)


def _record_modis_selection(diagnostics, best, terra, aqua, parameter, scene_datetime):
    if best is aqua and _choose_best_candidate([terra]) is None:
        diagnostics.warn("alternate_collection_used",
            f"Using MODIS collection {best['collection']} for {parameter} because the Terra lookup was unavailable.",
            parameter=parameter, collection=best["collection"])
    timestamp = best.get("image_time_utc")
    if timestamp:
        acquired = datetime.fromisoformat(timestamp.replace("Z", "+00:00"))
        if acquired.tzinfo is None:
            acquired = acquired.replace(tzinfo=timezone.utc)
        used_date = acquired.astimezone(timezone.utc).date()
        if used_date != scene_datetime.date():
            diagnostics.warn("alternate_date_used",
                f"Using MODIS {parameter} from {used_date.isoformat()} instead of acquisition date {scene_datetime.date().isoformat()}.",
                parameter=parameter, collection=best["collection"], image_time_utc=timestamp,
                requested_date=scene_datetime.date().isoformat(), date_used=used_date.isoformat(),
                abs_time_diff_hours=best.get("abs_time_diff_hours"))


@_logged_operation("fetch_atmosphere", inputs=(), outputs=(), worker_progress=True)
def fetch_modis_water_vapor_for_bbox(
    *,
    scene_datetime_utc: datetime,
    min_lon: float,
    min_lat: float,
    max_lon: float,
    max_lat: float,
    ee: Any = None,
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
    fallback_values: dict[str, float] | None = None,
    log_to_console: bool = False,
    scene_basename: str | None = None,
) -> ModisWaterVaporEstimate:
    """Fetch MODIS water vapor and AOD-derived inputs for a bbox.
    Args:
        scene_datetime_utc: Scene acquisition datetime in UTC.
        min_lon: Minimum longitude.
        min_lat: Minimum latitude.
        max_lon: Maximum longitude.
        max_lat: Maximum latitude.
        ee: Optional initialized Earth Engine module or client.
        ee_project: Optional Earth Engine project id.
        authenticate: Whether to run Earth Engine authentication.
        env_file: Optional env file containing Earth Engine settings.
        hours_window: Search window in hours around the scene time.
        terra_collection: Terra collection id.
        aqua_collection: Aqua collection id.
        band_name: Water vapor band name.
        aod_band_candidates: Candidate AOD band names.
        aod_scale_factor: Scale factor for AOD values.
        modis_scale_factor: Scale factor for MODIS water vapor values.
        modtran_baseline_water_vapor: Baseline water vapor used to derive the FLAASH preset.
        reduce_scale_m: Reduce-region scale in meters.
        fallback_values: Explicit defaults for missing aot550, water_vapor (g/cm^2), and ozone_cm_atm (cm-atm). None leaves missing values null.
        log_to_console: Whether to emit console logs.
        scene_basename: Optional scene basename for log prefixes.
    Returns:
        MODIS atmosphere lookup result.
    """

    fallbacks = _validate_fallback_values(fallback_values)
    diagnostics = _AtmosphereDiagnostics("modis_gee", log_to_console, scene_basename)
    _log("Querying MODIS atmosphere data", enabled=log_to_console,
         step="fetch_atmosphere", scene_basename=scene_basename)
    if scene_datetime_utc.tzinfo is None:
        scene_datetime_utc = scene_datetime_utc.replace(tzinfo=timezone.utc)
    else:
        scene_datetime_utc = scene_datetime_utc.astimezone(timezone.utc)
    if ee is None:
        ee = diagnostics.call(_init_ee_client, ee_project=ee_project,
            authenticate=authenticate, env_file=env_file, lookup={"operation": "initialize_earth_engine"})
    aod_band_candidates = aod_band_candidates or [
        "Aerosol_Optical_Depth_Land_Ocean_Mean",
        "Aerosol_Optical_Depth_Land_Ocean_Mean_Mean",
        "Aerosol_Optical_Depth_Land_Ocean",
    ]
    bounds = dict(min_lon=min_lon, min_lat=min_lat, max_lon=max_lon, max_lat=max_lat)

    def lookup(collection, parameter, bands, function):
        requested = [bands] if isinstance(bands, str) else bands
        details = dict(collection=collection, parameter=parameter, bands_requested=requested,
            scene_datetime_utc=scene_datetime_utc.isoformat(), hours_window=hours_window, **bounds)
        candidate = diagnostics.call(function, ee, collection, bands,
            scene_dt_utc=scene_datetime_utc, hours_window=hours_window,
            reduce_scale_m=reduce_scale_m, lookup=details, **bounds)
        _record_modis_lookup(diagnostics, candidate, details)
        return candidate

    terra = lookup(terra_collection, "water_vapor", band_name, _fetch_collection_value)
    aqua = lookup(aqua_collection, "water_vapor", band_name, _fetch_collection_value)
    best = _choose_best_candidate([terra, aqua])
    if best is None:
        diagnostics.warn("no_atmosphere_match", "No MODIS water-vapor match found in either collection within the configured time window.")
        diagnostics.warn("lookup_skipped", "MODIS AOD lookup was not attempted because no water-vapor match was found.", parameter="aot550")
        values, used = _apply_atmosphere_defaults(dict.fromkeys(_ATMOSPHERE_PARAMETERS), fallbacks, diagnostics)
        return ModisWaterVaporEstimate(status="no_modis_match", selected_collection=None,
            modis_raw_value=None, modis_scaled_g_cm2=None, water_vapor_preset=None,
            aod_raw=None, aod_scaled=None, default_visibility=None, modtran_atm=None,
            image_time_utc=None, abs_time_diff_hours=None, **values, defaults_used=used,
            scene_datetime_utc=scene_datetime_utc.isoformat(), **diagnostics.metadata())

    _record_modis_selection(diagnostics, best, terra, aqua, "water_vapor", scene_datetime_utc)
    raw_value = float(best["raw_value"])
    scaled = raw_value * modis_scale_factor
    preset = scaled / modtran_baseline_water_vapor
    terra_aod = lookup(terra_collection, "aot550", aod_band_candidates, _fetch_first_available_band_value)
    aqua_aod = lookup(aqua_collection, "aot550", aod_band_candidates, _fetch_first_available_band_value)
    best_aod = _choose_best_candidate([terra_aod, aqua_aod])
    if best_aod is not None:
        _record_modis_selection(diagnostics, best_aod, terra_aod, aqua_aod, "aot550", scene_datetime_utc)
    aod_raw = float(best_aod["raw_value"]) if best_aod is not None else None
    aod_scaled = aod_raw * aod_scale_factor if aod_raw is not None else None
    default_visibility = _visibility_from_aod(aod_scaled) if aod_scaled is not None else None
    modtran_atm = _choose_modtran_atm_by_lat_month((min_lat + max_lat) / 2.0, scene_datetime_utc.month)
    values, used = _apply_atmosphere_defaults(
        dict(aot550=aod_scaled, water_vapor=scaled, ozone_cm_atm=None), fallbacks, diagnostics)
    result = ModisWaterVaporEstimate(status="ok", selected_collection=best["collection"],
        modis_raw_value=raw_value, modis_scaled_g_cm2=scaled, water_vapor_preset=preset,
        aod_raw=aod_raw, aod_scaled=aod_scaled, default_visibility=default_visibility,
        modtran_atm=modtran_atm, image_time_utc=best["image_time_utc"],
        abs_time_diff_hours=best["abs_time_diff_hours"], **values, defaults_used=used,
        scene_datetime_utc=scene_datetime_utc.isoformat(), **diagnostics.metadata())
    _log("Fetched MODIS atmosphere values", enabled=log_to_console,
         step="fetch_atmosphere", scene_basename=scene_basename)
    return result


__all__ = [
    "FetchAtmosphere",
    "AtmosphereEstimate",
    "ModisWaterVaporEstimate",
    "fetch_modis_water_vapor_for_bbox",
    "fetch_power_atmosphere_for_bbox",
]


class FetchAtmosphere(FunctionPlugin):
    input_dependency_paths = frozenset({"env_file"})
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset({"output_path"})
    output_dependency_paths = output_path_resolution_paths
    output_parent_creation_paths = output_path_resolution_paths
    output_collision_check_paths = output_path_resolution_paths
    output_reuse_paths = output_path_resolution_paths
    output_validation_paths = output_path_resolution_paths
    output_invalid_removal_paths = output_path_resolution_paths
    output_temporary_cleanup_paths = output_path_resolution_paths
    output_hpc_staging_paths = output_path_resolution_paths
    output_hpc_download_paths = output_path_resolution_paths

    def run(self, *, params, shared):
        from vhrharmonize.io.metadata import write_json

        params = dict(params)
        source = params.pop("source", shared.get("source", "nasa_power"))
        if source == "nasa_power":
            self.target = "vhrharmonize.plugins.fetch_atmosphere:fetch_power_atmosphere_for_bbox"
            params["day_utc"] = date.fromisoformat(
                str(params.get("day_utc", shared.get("day_utc")))
            )
        elif source == "modis_gee":
            self.target = "vhrharmonize.plugins.fetch_atmosphere:fetch_modis_water_vapor_for_bbox"
            value = params.get("scene_datetime_utc", shared.get("scene_datetime_utc"))
            params["scene_datetime_utc"] = datetime.fromisoformat(str(value).replace("Z", "+00:00"))
        else:
            raise ValueError(f"Unsupported atmosphere source: {source}")
        output_path = params.pop("output_path", None)
        if not output_path:
            raise ValueError("fetch_atmosphere requires param:output_path")
        try:
            result = super().run(params=params, shared=shared)
        except Exception as error:
            diagnostics = getattr(error, "atmosphere_diagnostics", None)
            if diagnostics is not None:
                write_json(output_path, dict(status="error", source=source,
                    error=f"{type(error).__name__}: {error}", **diagnostics))
            raise
        write_json(output_path, result)
        return result
