"""NASA atmosphere parameter fetchers for Py6S configuration."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, timedelta
from typing import Dict, List, Optional, Tuple

import requests


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
    except Exception:
        return None
    # POWER uses negative sentinels (for example -999) for missing values.
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

    aod_map = param.get("AOD_55", {}) or {}
    pw_map = param.get("PW", {}) or {}
    to3_map = param.get("TO3", {}) or {}

    aod = _parse_power_value(aod_map.get(ymd))
    pw_cm = _parse_power_value(pw_map.get(ymd))
    to3_du = _parse_power_value(to3_map.get(ymd))

    ozone_cm_atm = (to3_du / 1000.0) if to3_du is not None else None
    # PW in POWER is precipitable water in cm, numerically equivalent to g/cm^2.
    water_vapor_g_cm2 = pw_cm

    return {
        "aot550": aod,
        "water_vapor": water_vapor_g_cm2,
        "ozone_cm_atm": ozone_cm_atm,
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
) -> AtmosphereEstimate:
    """Fetch AOT550, water vapor, and ozone for a bbox via NASA POWER daily API.

    The bbox is sampled by an NxN point grid and aggregated by mean over valid points.
    If no valid values exist for the scene date, dates in +/- ``search_days`` are tried.
    """
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
            return AtmosphereEstimate(
                aot550=(sum(aod_vals) / len(aod_vals)) if aod_vals else None,
                water_vapor=(sum(wv_vals) / len(wv_vals)) if wv_vals else None,
                ozone_cm_atm=(sum(oz_vals) / len(oz_vals)) if oz_vals else None,
                source="nasa_power_daily_point_grid",
                sample_count=len(points),
                date_used=query_day.isoformat(),
            )

    return AtmosphereEstimate(
        aot550=None,
        water_vapor=None,
        ozone_cm_atm=None,
        source="nasa_power_daily_point_grid",
        sample_count=len(points),
        date_used=day_utc.isoformat(),
    )


__all__ = [
    "AtmosphereEstimate",
    "fetch_power_atmosphere_for_bbox",
]
