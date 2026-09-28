"""Decode metadata file formats without interpreting sensor fields."""
from __future__ import annotations
import json
import os
import re
from pathlib import Path
from typing import Any, Dict, Iterator, List
import xml.etree.ElementTree as ET
import yaml

def _split_imd_statements(text: str) -> Iterator[str]:
    """Split assignments while preserving quoted delimiters and multiline arrays."""
    buffer = []
    quote = None
    depth = 0
    escaped = False
    for char in text:
        if quote:
            buffer.append(char)
            if char == quote and not escaped:
                quote = None
            escaped = char == "\\" and not escaped
            continue
        if char in {"'", '"'}:
            quote = char
        elif char in "([{":
            depth += 1
        elif char in ")]}":
            depth -= 1
        statement = "".join(buffer).strip()
        if char == ";" and depth == 0 or char == "\n" and depth == 0 and statement.startswith(("BEGIN_GROUP", "END_GROUP", "BEGIN_OBJECT", "END_OBJECT")):
            if statement:
                yield statement
            buffer = []
        else:
            buffer.append(char)
    if quote or depth:
        raise ValueError("Unterminated metadata string or array")
    statement = "".join(buffer).strip()
    if statement:
        yield statement


def _split_top_level_csv(raw: str) -> List[str]:
    """Split a top-level CSV-like string.
    Args:
        raw: Raw comma-delimited string.
    Returns:
        Top-level comma-delimited tokens.
    """
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
    """Parse a scalar IMD value.
    Args:
        raw_value: Raw IMD value text.
    Returns:
        Parsed Python scalar or nested list value.
    """
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
    """Append a parsed group into a container.
    Args:
        container: Target parsed metadata container.
        key: Group key to append.
        value: Parsed group payload.
    Returns:
        None.
    """
    existing = container.get(key)
    if existing is None:
        container[key] = value
    elif isinstance(existing, list):
        existing.append(value)
    else:
        container[key] = [existing, value]


def parse_imd_text(imd_text: str) -> Dict[str, Any]:
    """Parse an IMD file into nested Python objects."""
    root: Dict[str, Any] = {}
    stack: List[Dict[str, Any]] = [root]
    groups = []

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
            groups.append(group_name)
            continue
        if key in {"END_GROUP", "END_OBJECT"}:
            if not groups or groups.pop() != str(_parse_scalar(raw_value)):
                raise ValueError(f"Mismatched metadata group: {raw_value}")
            stack.pop()
            continue
        stack[-1][key] = _parse_scalar(raw_value)

    if groups:
        raise ValueError(f"Unclosed metadata groups: {groups}")
    return root


def parse_imd_file(imd_file: str) -> Dict[str, Any]:
    """Parse a IMD file into nested Python objects."""
    with open(imd_file, "r", encoding="utf-8") as handle:
        return parse_imd_text(handle.read())



def read_metadata(path: str) -> dict:
    suffix = Path(path).suffix.lower()
    if suffix == ".imd":
        result = parse_imd_file(path)
    elif suffix in {".yaml", ".yml"}:
        result = yaml.safe_load(Path(path).read_text())
    elif suffix == ".xml":
        def decode(element):
            if not len(element):
                return _parse_scalar(element.text or "")
            result = dict(element.attrib)
            for child in element:
                key = child.tag.rsplit("}", 1)[-1]
                value = decode(child)
                if key in result:
                    if not isinstance(result[key], list):
                        result[key] = [result[key]]
                    result[key].append(value)
                else:
                    result[key] = value
            return result
        root = ET.parse(path).getroot()
        result = {root.tag.rsplit("}", 1)[-1]: decode(root)}
    elif suffix in {".json", ".geojson"}:
        result = json.loads(Path(path).read_text())
    else:
        raise ValueError(f"Unsupported metadata format: {path}")
    if not isinstance(result, dict):
        raise ValueError(f"Metadata must decode to an object: {path}")
    # Enforce the contract: metadata passed between plugins is always JSON.
    return json.loads(json.dumps(result, default=str, allow_nan=False))


def materialize_geometry(geometry: dict, *, epsg: int = 4326):
    """Materialize GeoJSON geometry in WGS84, optionally reprojecting it."""
    from shapely.geometry import shape
    from shapely.ops import transform
    from pyproj import Transformer
    result = shape(geometry)
    if result.is_empty or not result.is_valid or result.area == 0:
        raise ValueError("Metadata geometry must be a valid nonempty polygon")
    if epsg != 4326:
        result = transform(Transformer.from_crs(4326, epsg, always_xy=True).transform, result)
    return result


def write_json(filename, value):
    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    temporary = f"{filename}.writing-{os.getpid()}"
    try:
        Path(temporary).write_text(json.dumps(value, indent=2, allow_nan=False, default=str) + "\n")
        os.replace(temporary, filename)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)

