"""Provider-agnostic standardized metadata models."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime
from typing import Any, Dict, Iterator, List, Mapping, Optional, Tuple

from vhrharmonize.providers.worldview import WorldViewMetadata, parse_worldview_basename

_WV_BAND_ORDER = [
    "BAND_C",
    "BAND_B",
    "BAND_G",
    "BAND_Y",
    "BAND_R",
    "BAND_RE",
    "BAND_N",
    "BAND_N2",
]

_WV_RADIOMETRIC_GAIN_OFFSET_2018V0: Dict[str, Dict[str, Tuple[float, float]]] = {
    "WV02": {
        "BAND_C": (0.977, -5.552),
        "BAND_B": (0.987, -4.443),
        "BAND_G": (0.987, -2.523),
        "BAND_Y": (0.987, -4.929),
        "BAND_R": (0.987, -3.914),
        "BAND_RE": (0.987, -5.675),
        "BAND_N": (0.987, -1.546),
        "BAND_N2": (0.987, -3.564),
    },
    "WV03": {
        "BAND_C": (0.905, -8.604),
        "BAND_B": (0.940, -5.809),
        "BAND_G": (0.938, -4.996),
        "BAND_Y": (0.962, -3.649),
        "BAND_R": (0.964, -3.021),
        "BAND_RE": (1.000, -4.521),
        "BAND_N": (0.961, -5.522),
        "BAND_N2": (0.978, -2.992),
    },
}


def _walk_values(data: Any) -> Iterator[tuple[str, Any]]:
    """Yield nested mapping values.
    Args:
        data: Nested mapping or list structure.
    Returns:
        Iterator of key-value pairs from nested mappings.
    """
    if isinstance(data, dict):
        for key, value in data.items():
            yield key, value
            yield from _walk_values(value)
    elif isinstance(data, list):
        for value in data:
            yield from _walk_values(value)


def _collect_band_groups(data: Mapping[str, Any]) -> List[Tuple[str, Mapping[str, Any]]]:
    """Collect nested WorldView band groups.
    Args:
        data: Nested metadata mapping.
    Returns:
        Ordered band-group tuples.
    """
    groups: List[Tuple[str, Mapping[str, Any]]] = []
    for key, value in _walk_values(data):
        if isinstance(value, dict) and key.upper().startswith("BAND_"):
            groups.append((key, value))
    return groups


@dataclass(frozen=True)
class StandardizedMetadata:
    """Provider-neutral metadata fields used by preprocessing workflows."""

    provider: str
    photo_basename: Optional[str] = None
    acquisition_datetime_utc: Optional[datetime] = None
    sensor_id: Optional[str] = None
    cloud_cover: Optional[float] = None
    sun_azimuth: Optional[float] = None
    sun_elevation: Optional[float] = None
    sun_zenith: Optional[float] = None
    view_azimuth: Optional[float] = None
    off_nadir_angle: Optional[float] = None
    view_zenith: Optional[float] = None
    line_of_sight_azimuth: Optional[float] = None
    line_of_sight_zenith: Optional[float] = None
    row_resolution: Optional[float] = None
    column_resolution: Optional[float] = None
    product_resolution: Optional[float] = None
    band_order: Tuple[str, ...] = ()
    dn_to_radiance_factors: Tuple[float, ...] = ()
    dn_to_radiance_gains: Tuple[float, ...] = ()
    dn_to_radiance_offsets: Tuple[float, ...] = ()
    radiometric_adjustment_version: Optional[str] = None
    provider_field_map: Mapping[str, str] = field(default_factory=dict)
    source_metadata: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_worldview_metadata(cls: type["StandardizedMetadata"], worldview_metadata: WorldViewMetadata) -> "StandardizedMetadata":
        """Build standardized metadata from WorldView metadata.
        Args:
            cls: Dataclass type being constructed.
            worldview_metadata: Parsed WorldView metadata object.
        Returns:
            Standardized metadata instance.
        """
        acquisition_dt = worldview_metadata.find_first_datetime(
            "firstLineTime", "TLCTime", "earliestAcqTime", "latestAcqTime"
        )
        sun_azimuth = worldview_metadata.find_first_number("meanSunAz")
        sun_elevation = worldview_metadata.find_first_number("meanSunEl")
        off_nadir = worldview_metadata.find_first_number("meanOffNadirViewAngle")
        sat_azimuth = worldview_metadata.find_first_number("meanSatAz")
        sensor_id = worldview_metadata.find_first("satId")
        cloud_cover_value = worldview_metadata.find_first_number("cloudCover")
        cloud_cover = None if cloud_cover_value is None else cloud_cover_value * 100.0

        band_factors: Dict[str, float] = {}
        discovered_band_order: List[str] = []
        for band_name, band_group in _collect_band_groups(worldview_metadata.raw_metadata):
            abs_cal = None
            effective_bandwidth = None
            if isinstance(band_group, dict):
                abs_cal = band_group.get("absCalFactor")
                effective_bandwidth = band_group.get("effectiveBandwidth")
            try:
                abs_cal = float(abs_cal) if abs_cal is not None else None
                effective_bandwidth = float(effective_bandwidth) if effective_bandwidth is not None else None
            except (TypeError, ValueError):
                abs_cal = None
                effective_bandwidth = None
            if abs_cal is None or effective_bandwidth in (None, 0.0):
                continue
            discovered_band_order.append(band_name)
            band_factors[band_name] = abs_cal / effective_bandwidth

        ordered_bands = [band for band in _WV_BAND_ORDER if band in band_factors]
        if not ordered_bands:
            ordered_bands = discovered_band_order

        gains: List[float] = []
        offsets: List[float] = []
        version: Optional[str] = None
        sensor_key = str(sensor_id).upper() if sensor_id is not None else None
        gain_offset_table = _WV_RADIOMETRIC_GAIN_OFFSET_2018V0.get(sensor_key or "")
        if gain_offset_table and all(band in gain_offset_table for band in ordered_bands):
            version = "2018v0"
            gains = [gain_offset_table[band][0] for band in ordered_bands]
            offsets = [gain_offset_table[band][1] for band in ordered_bands]

        return cls(
            provider="worldview",
            photo_basename=worldview_metadata.photo_basename,
            acquisition_datetime_utc=acquisition_dt,
            sensor_id=str(sensor_id) if sensor_id is not None else None,
            cloud_cover=cloud_cover,
            sun_azimuth=sun_azimuth,
            sun_elevation=sun_elevation,
            sun_zenith=(90.0 - sun_elevation) if sun_elevation is not None else None,
            view_azimuth=sat_azimuth,
            off_nadir_angle=off_nadir,
            view_zenith=off_nadir,
            line_of_sight_azimuth=sat_azimuth,
            line_of_sight_zenith=(180.0 - off_nadir) if off_nadir is not None else None,
            row_resolution=worldview_metadata.find_first_number("meanProductRowGSD"),
            column_resolution=worldview_metadata.find_first_number("meanProductColGSD"),
            product_resolution=worldview_metadata.find_first_number("meanProductGSD"),
            band_order=tuple(ordered_bands),
            dn_to_radiance_factors=tuple(band_factors[band] for band in ordered_bands),
            dn_to_radiance_gains=tuple(gains),
            dn_to_radiance_offsets=tuple(offsets),
            radiometric_adjustment_version=version,
            provider_field_map={
                "satId": "sensor_id",
                "firstLineTime": "acquisition_datetime_utc",
                "TLCTime": "acquisition_datetime_utc",
                "earliestAcqTime": "acquisition_datetime_utc",
                "latestAcqTime": "acquisition_datetime_utc",
                "cloudCover": "cloud_cover",
                "meanSunAz": "sun_azimuth",
                "meanSunEl": "sun_elevation",
                "meanSatAz": "view_azimuth",
                "meanOffNadirViewAngle": "off_nadir_angle",
                "meanProductRowGSD": "row_resolution",
                "meanProductColGSD": "column_resolution",
                "meanProductGSD": "product_resolution",
            },
            source_metadata=worldview_metadata.raw_metadata,
        )

    def resolve_scene_datetime(self) -> datetime:
        """Resolve a scene acquisition datetime.
        Args:
            self: Standardized metadata instance.
        Returns:
            Resolved acquisition datetime in UTC.
        """
        if self.acquisition_datetime_utc is not None:
            return self.acquisition_datetime_utc
        if not self.photo_basename:
            raise ValueError("No acquisition datetime or photo basename available.")
        parts = parse_worldview_basename(self.photo_basename)
        if parts is None:
            raise ValueError(f"Could not parse acquisition datetime from basename: {self.photo_basename}")
        return parts.acquisition_datetime_utc

    def to_dict(self) -> Dict[str, Any]:
        """Convert standardized metadata to a dictionary.
        Args:
            self: Standardized metadata instance.
        Returns:
            Dictionary representation of the metadata.
        """
        payload = asdict(self)
        if self.acquisition_datetime_utc is not None:
            payload["acquisition_datetime_utc"] = self.acquisition_datetime_utc.isoformat()
        return payload


__all__ = ["StandardizedMetadata"]
