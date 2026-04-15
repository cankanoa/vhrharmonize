"""WorldView IMD metadata parsing."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Dict, Iterator, List, Mapping, Optional, Tuple


def _split_imd_statements(text: str) -> Iterator[str]:
    buffer: List[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        buffer.append(line)
        if ";" not in line and not line.startswith(("BEGIN_GROUP", "END_GROUP", "BEGIN_OBJECT", "END_OBJECT")):
            continue
        statement = " ".join(buffer).strip()
        buffer.clear()
        if statement.endswith(";"):
            statement = statement[:-1].strip()
        if statement:
            yield statement
    if buffer:
        statement = " ".join(buffer).strip()
        if statement.endswith(";"):
            statement = statement[:-1].strip()
        if statement:
            yield statement


def _split_top_level_csv(raw: str) -> List[str]:
    parts: List[str] = []
    buffer: List[str] = []
    depth = 0
    in_quotes = False
    quote_char = ""
    for char in raw:
        if char in {'"', "'"}:
            if in_quotes and char == quote_char:
                in_quotes = False
                quote_char = ""
            elif not in_quotes:
                in_quotes = True
                quote_char = char
        elif not in_quotes:
            if char in "([{":
                depth += 1
            elif char in ")]}":
                depth = max(0, depth - 1)
            elif char == "," and depth == 0:
                parts.append("".join(buffer).strip())
                buffer.clear()
                continue
        buffer.append(char)
    if buffer:
        parts.append("".join(buffer).strip())
    return [part for part in parts if part]


def _parse_scalar(raw_value: str) -> Any:
    value = raw_value.strip()
    if not value:
        return ""
    if (value.startswith('"') and value.endswith('"')) or (value.startswith("'") and value.endswith("'")):
        return value[1:-1]
    upper = value.upper()
    if upper == "TRUE":
        return True
    if upper == "FALSE":
        return False
    if upper in {"NULL", "NONE"}:
        return None
    if value.startswith(("(", "[", "{")) and value.endswith((")", "]", "}")):
        inner = value[1:-1].strip()
        if not inner:
            return []
        return [_parse_scalar(part) for part in _split_top_level_csv(inner)]
    if re.fullmatch(r"[-+]?\d+", value):
        try:
            return int(value)
        except ValueError:
            return value
    if re.fullmatch(r"[-+]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", value):
        try:
            return float(value)
        except ValueError:
            return value
    return value


def _append_group(container: Dict[str, Any], key: str, value: Dict[str, Any]) -> None:
    existing = container.get(key)
    if existing is None:
        container[key] = value
    elif isinstance(existing, list):
        existing.append(value)
    else:
        container[key] = [existing, value]


def parse_worldview_imd_text(imd_text: str) -> Dict[str, Any]:
    """Parse an IMD file into nested Python objects."""
    root: Dict[str, Any] = {}
    stack: List[Dict[str, Any]] = [root]

    for statement in _split_imd_statements(imd_text):
        if "=" not in statement:
            continue
        key, raw_value = statement.split("=", 1)
        key = key.strip()
        raw_value = raw_value.strip()

        if key in {"BEGIN_GROUP", "BEGIN_OBJECT"}:
            group_name = str(_parse_scalar(raw_value))
            group_payload: Dict[str, Any] = {}
            _append_group(stack[-1], group_name, group_payload)
            stack.append(group_payload)
            continue
        if key in {"END_GROUP", "END_OBJECT"}:
            if len(stack) > 1:
                stack.pop()
            continue
        stack[-1][key] = _parse_scalar(raw_value)

    return root


def parse_worldview_imd_file(imd_file: str) -> Dict[str, Any]:
    """Parse a WorldView IMD file into nested Python objects."""
    with open(imd_file, "r", encoding="utf-8") as handle:
        return parse_worldview_imd_text(handle.read())


def _walk_values(data: Any) -> Iterator[Tuple[str, Any]]:
    if isinstance(data, dict):
        for key, value in data.items():
            yield key, value
            yield from _walk_values(value)
    elif isinstance(data, list):
        for value in data:
            yield from _walk_values(value)


@dataclass(frozen=True)
class WorldViewMetadata:
    """Provider-specific parsed IMD metadata."""

    imd_file: str
    photo_basename: str
    raw_metadata: Mapping[str, Any]

    @classmethod
    def from_imd_file(cls, imd_file: str, *, photo_basename: Optional[str] = None) -> "WorldViewMetadata":
        basename = photo_basename or os.path.splitext(os.path.basename(imd_file))[0]
        return cls(
            imd_file=imd_file,
            photo_basename=basename,
            raw_metadata=parse_worldview_imd_file(imd_file),
        )

    def find_first(self, *keys: str) -> Any:
        for wanted_key in keys:
            for current_key, current_value in _walk_values(self.raw_metadata):
                if current_key == wanted_key:
                    return current_value
        return None

    def find_first_number(self, *keys: str) -> Optional[float]:
        value = self.find_first(*keys)
        if value is None:
            return None
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    def find_first_datetime(self, *keys: str) -> Optional[datetime]:
        raw_value = self.find_first(*keys)
        if raw_value in (None, ""):
            return None
        try:
            parsed = datetime.fromisoformat(str(raw_value).replace("Z", "+00:00"))
        except ValueError:
            return None
        if parsed.tzinfo is None:
            parsed = parsed.replace(tzinfo=timezone.utc)
        return parsed.astimezone(timezone.utc)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "imd_file": self.imd_file,
            "photo_basename": self.photo_basename,
            "raw_metadata": dict(self.raw_metadata),
        }


def load_worldview_metadata(imd_file: str, *, photo_basename: Optional[str] = None) -> WorldViewMetadata:
    """Load provider-specific WorldView metadata from an IMD file."""
    return WorldViewMetadata.from_imd_file(imd_file, photo_basename=photo_basename)


__all__ = [
    "WorldViewMetadata",
    "load_worldview_metadata",
    "parse_worldview_imd_file",
    "parse_worldview_imd_text",
]
