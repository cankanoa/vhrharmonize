"""WorldView provider discovery and IMD metadata parsing."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Tuple


@dataclass(frozen=True)
class WorldViewFilenameParts:
    """Parsed components of a WorldView bundle filename."""

    basename: str
    acquisition_token: str
    product_code: str
    catalog_id: str
    sequence_id: str
    scene_id: str

    @property
    def image_role(self) -> Optional[str]:
        """Return the parsed image role.
        Args:
            self: Parsed filename parts instance.
        Returns:
            WorldView image role string or None.
        """
        product = self.product_code.upper()
        if product.startswith("M"):
            return "mul"
        if product.startswith("P"):
            return "pan"
        return None

    @property
    def acquisition_datetime_utc(self) -> datetime:
        """Return the acquisition datetime in UTC.
        Args:
            self: Parsed filename parts instance.
        Returns:
            Parsed acquisition datetime.
        """
        month_map = {
            "JAN": 1,
            "FEB": 2,
            "MAR": 3,
            "APR": 4,
            "MAY": 5,
            "JUN": 6,
            "JUL": 7,
            "AUG": 8,
            "SEP": 9,
            "OCT": 10,
            "NOV": 11,
            "DEC": 12,
        }
        match = re.match(r"^(\d{2})([A-Z]{3})(\d{2})(\d{2})(\d{2})(\d{2})$", self.acquisition_token)
        if not match:
            raise ValueError(f"Could not parse WorldView acquisition token: {self.acquisition_token}")
        return datetime(
            year=2000 + int(match.group(3)),
            month=month_map[match.group(2)],
            day=int(match.group(1)),
            hour=int(match.group(4)),
            minute=int(match.group(5)),
            second=int(match.group(6)),
            tzinfo=timezone.utc,
        )


@dataclass
class WorldViewImage:
    """WorldView image files and metadata for a single basename."""

    filename_parts: WorldViewFilenameParts
    root_folder_path: str
    tif_file: str
    imd_file: Optional[str]
    shp_file: Optional[str]
    til_file: Optional[str]
    worldview_metadata: Optional[object] = None
    standardized_metadata: Optional[object] = None
    step_file_paths: Dict[str, str] = field(default_factory=dict)

    @property
    def basename(self) -> str:
        """Return the image basename.
        Args:
            self: WorldView image instance.
        Returns:
            Parsed image basename.
        """
        return self.filename_parts.basename

    @property
    def image_role(self) -> Optional[str]:
        """Return the image role.
        Args:
            self: WorldView image instance.
        Returns:
            Parsed image role string or None.
        """
        return self.filename_parts.image_role


@dataclass
class WorldViewScene:
    """A single WorldView acquisition bundle with one MUL image and one PAN image."""

    scene_id: str
    catalog_id: str
    root_folder_path: str
    mul_image: Optional[WorldViewImage] = None
    pan_image: Optional[WorldViewImage] = None
    step_outputs: Dict[str, List[str]] = field(default_factory=dict)
    metadata_report_path: Optional[str] = None

    @property
    def primary_image(self) -> Optional[WorldViewImage]:
        """Return the preferred primary image.
        Args:
            self: WorldView scene instance.
        Returns:
            Primary scene image or None.
        """
        return self.mul_image or self.pan_image

    @property
    def primary_basename(self) -> Optional[str]:
        """Return the primary scene basename.
        Args:
            self: WorldView scene instance.
        Returns:
            Primary image basename or None.
        """
        image = self.primary_image
        return image.basename if image else None

    def get_image(self, role: str) -> Optional[WorldViewImage]:
        """Return an image for a requested role.
        Args:
            role: Requested image role.
        Returns:
            Matching scene image or None.
        """
        role_norm = role.strip().lower()
        if role_norm == "mul":
            return self.mul_image
        if role_norm == "pan":
            return self.pan_image
        return None

    def set_image(self, image: WorldViewImage) -> None:
        """Assign an image to the scene.
        Args:
            image: WorldView image to attach to the scene.
        Returns:
            None.
        """
        if image.image_role == "mul":
            if self.mul_image is not None and self.mul_image.basename != image.basename:
                raise ValueError(
                    f"Scene {self.scene_id}_{self.catalog_id} has more than one multispectral image: "
                    f"{self.mul_image.basename}, {image.basename}"
                )
            self.mul_image = image
            return
        if image.image_role == "pan":
            if self.pan_image is not None and self.pan_image.basename != image.basename:
                raise ValueError(
                    f"Scene {self.scene_id}_{self.catalog_id} has more than one panchromatic image: "
                    f"{self.pan_image.basename}, {image.basename}"
                )
            self.pan_image = image
            return
        raise ValueError(f"Unsupported WorldView image role for scene assignment: {image.image_role}")

    def iter_images(self) -> List[WorldViewImage]:
        """Return available scene images.
        Args:
            self: WorldView scene instance.
        Returns:
            List of non-null scene images.
        """
        return [image for image in (self.mul_image, self.pan_image) if image is not None]

    def to_dict(self) -> Dict[str, object]:
        """Convert a scene to a flat dictionary.
        Args:
            self: WorldView scene instance.
        Returns:
            Flattened scene dictionary.
        """
        mul_image = self.mul_image
        pan_image = self.pan_image
        payload: Dict[str, object] = {
            "scene_id": self.scene_id,
            "catalog_id": self.catalog_id,
            "root_folder_path": self.root_folder_path,
            "step_outputs": dict(self.step_outputs),
        }
        if mul_image is not None:
            payload.update(
                {
                    "mul_photo_basename": mul_image.basename,
                    "mul_imd_file": mul_image.imd_file,
                    "mul_tif_file": mul_image.tif_file,
                    "mul_shp_file": mul_image.shp_file,
                    "mul_til_file": mul_image.til_file,
                }
            )
        if pan_image is not None:
            payload.update(
                {
                    "pan_photo_basename": pan_image.basename,
                    "pan_imd_file": pan_image.imd_file,
                    "pan_tif_file": pan_image.tif_file,
                    "pan_shp_file": pan_image.shp_file,
                    "pan_til_file": pan_image.til_file,
                }
            )
        return payload


def parse_worldview_basename(name: str) -> Optional[WorldViewFilenameParts]:
    """Parse a WorldView basename using `-` and `_` separators."""
    basename = os.path.splitext(os.path.basename(name))[0]
    parts = re.split(r"[-_]", basename)
    if len(parts) < 5:
        return None
    acquisition_token, product_code, catalog_id, sequence_id, scene_id = parts[:5]
    if not re.fullmatch(r"\d{12}", catalog_id):
        return None
    if not re.fullmatch(r"P\d{3}", scene_id, flags=re.IGNORECASE):
        return None
    return WorldViewFilenameParts(
        basename=basename,
        acquisition_token=acquisition_token,
        product_code=product_code,
        catalog_id=catalog_id,
        sequence_id=sequence_id,
        scene_id=scene_id.upper(),
    )


def _scene_root_for_path(path: str) -> str:
    """Resolve the canonical scene root for a file path.
    Args:
        path: Candidate scene file path.
    Returns:
        Scene root directory path.
    """
    directory = os.path.dirname(path)
    base = os.path.basename(directory).upper()
    if base.endswith("_MUL") or base.endswith("_PAN"):
        return os.path.dirname(directory)
    return directory


def _find_worldview_shp_file(
    image_directory: str,
    scene_root: str,
    basename: str,
) -> Optional[str]:
    """Find a scene shapefile for a basename.
    Args:
        image_directory: Directory containing the image bundle.
        scene_root: Root scene directory.
        basename: Image basename to match.
    Returns:
        Matching shapefile path or None.
    """
    for candidate in sorted(os.listdir(image_directory)):
        candidate_path = os.path.join(image_directory, candidate)
        if not os.path.isfile(candidate_path):
            continue
        if os.path.splitext(candidate)[1].lower() != ".shp":
            continue
        if basename in os.path.splitext(candidate)[0]:
            return candidate_path

    gis_files_directory = None
    for candidate in sorted(os.listdir(scene_root)):
        candidate_path = os.path.join(scene_root, candidate)
        if os.path.isdir(candidate_path) and candidate.upper() == "GIS_FILES":
            gis_files_directory = candidate_path
            break

    if gis_files_directory is None:
        return None

    for candidate in sorted(os.listdir(gis_files_directory)):
        candidate_path = os.path.join(gis_files_directory, candidate)
        if not os.path.isfile(candidate_path):
            continue
        if os.path.splitext(candidate)[1].lower() != ".shp":
            continue
        if basename in os.path.splitext(candidate)[0]:
            return candidate_path

    return None


def _companion_files_for_stem(directory: str, scene_root: str, basename: str) -> Dict[str, Optional[str]]:
    """Collect companion files for an image basename.
    Args:
        directory: Directory containing image bundle files.
        scene_root: Root scene directory.
        basename: Image basename to match.
    Returns:
        Mapping of discovered companion file paths.
    """
    companions = {"imd_file": None, "tif_file": None, "shp_file": None, "til_file": None}
    for candidate in os.listdir(directory):
        candidate_path = os.path.join(directory, candidate)
        if not os.path.isfile(candidate_path):
            continue
        stem, extension = os.path.splitext(candidate)
        if stem != basename:
            continue
        ext = extension.lower()
        if ext == ".imd":
            companions["imd_file"] = candidate_path
        elif ext == ".tif":
            companions["tif_file"] = candidate_path
        elif ext == ".til":
            companions["til_file"] = candidate_path
    companions["shp_file"] = _find_worldview_shp_file(directory, scene_root, basename)
    return companions


def discover_worldview_scene_tree_from_tif_files(
    tif_files: Iterable[str],
    filter_basenames: Optional[Iterable[str]] = None,
) -> Dict[str, Dict[str, WorldViewScene]]:
    """Discover WorldView scenes from tif files grouped by scene id then catalog id."""
    allowed = {value.strip() for value in (filter_basenames or []) if value and value.strip()}
    scene_tree: Dict[str, Dict[str, WorldViewScene]] = {}

    for tif_file in sorted({os.path.abspath(path) for path in tif_files}):
        filename_parts = parse_worldview_basename(tif_file)
        if filename_parts is None:
            continue
        if allowed and filename_parts.basename not in allowed:
            continue

        scene_root = _scene_root_for_path(tif_file)
        companions = _companion_files_for_stem(
            os.path.dirname(tif_file),
            scene_root,
            filename_parts.basename,
        )
        scene_bucket = scene_tree.setdefault(filename_parts.scene_id, {})
        scene = scene_bucket.setdefault(
            filename_parts.catalog_id,
            WorldViewScene(
                scene_id=filename_parts.scene_id,
                catalog_id=filename_parts.catalog_id,
                root_folder_path=scene_root,
            ),
        )
        image = WorldViewImage(
            filename_parts=filename_parts,
            root_folder_path=scene_root,
            tif_file=tif_file,
            imd_file=companions["imd_file"],
            shp_file=companions["shp_file"],
            til_file=companions["til_file"],
        )
        scene.set_image(image)

    return scene_tree


def iter_worldview_scenes(scene_tree: Mapping[str, Mapping[str, WorldViewScene]]) -> List[WorldViewScene]:
    """Flatten a scene tree into a stable list of scenes."""
    scenes: List[WorldViewScene] = []
    for scene_id in sorted(scene_tree):
        for catalog_id in sorted(scene_tree[scene_id]):
            scenes.append(scene_tree[scene_id][catalog_id])
    return scenes


def _split_imd_statements(text: str) -> Iterator[str]:
    """Split IMD text into logical statements.
    Args:
        text: Raw IMD text.
    Returns:
        Iterator of normalized IMD statements.
    """
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


@dataclass(frozen=True)
class WorldViewMetadata:
    """Provider-specific parsed IMD metadata."""

    imd_file: str
    photo_basename: str
    raw_metadata: Mapping[str, Any]

    @classmethod
    def from_imd_file(cls: type["WorldViewMetadata"], imd_file: str, *, photo_basename: Optional[str] = None) -> "WorldViewMetadata":
        """Build metadata from an IMD file.
        Args:
            cls: Dataclass type being constructed.
            imd_file: IMD file path.
            photo_basename: Optional basename override.
        Returns:
            Parsed WorldView metadata instance.
        """
        basename = photo_basename or os.path.splitext(os.path.basename(imd_file))[0]
        return cls(
            imd_file=imd_file,
            photo_basename=basename,
            raw_metadata=parse_worldview_imd_file(imd_file),
        )

    def find_first(self, *keys: str) -> Any:
        """Find the first matching raw metadata value.
        Args:
            keys: Metadata keys to search for in order.
        Returns:
            First matching value or None.
        """
        for wanted_key in keys:
            for current_key, current_value in _walk_values(self.raw_metadata):
                if current_key == wanted_key:
                    return current_value
        return None

    def find_first_number(self, *keys: str) -> Optional[float]:
        """Find the first matching numeric metadata value.
        Args:
            keys: Metadata keys to search for in order.
        Returns:
            First matching numeric value or None.
        """
        value = self.find_first(*keys)
        if value is None:
            return None
        return float(value)

    def find_first_datetime(self, *keys: str) -> Optional[datetime]:
        """Find the first matching datetime metadata value.
        Args:
            keys: Metadata keys to search for in order.
        Returns:
            First matching UTC datetime or None.
        """
        raw_value = self.find_first(*keys)
        if raw_value in (None, ""):
            return None
        parsed = datetime.fromisoformat(str(raw_value).replace("Z", "+00:00"))
        if parsed.tzinfo is None:
            parsed = parsed.replace(tzinfo=timezone.utc)
        return parsed.astimezone(timezone.utc)

    def to_dict(self) -> Dict[str, Any]:
        """Convert parsed metadata to a dictionary.
        Args:
            self: WorldView metadata instance.
        Returns:
            Dictionary representation of the metadata.
        """
        return {
            "imd_file": self.imd_file,
            "photo_basename": self.photo_basename,
            "raw_metadata": dict(self.raw_metadata),
        }


def load_worldview_metadata(imd_file: str, *, photo_basename: Optional[str] = None) -> WorldViewMetadata:
    """Load provider-specific WorldView metadata from an IMD file."""
    return WorldViewMetadata.from_imd_file(imd_file, photo_basename=photo_basename)


def enrich_worldview_scenes_with_metadata(scenes: Iterable[WorldViewScene]) -> List[WorldViewScene]:
    """Attach provider-specific and standardized metadata to discovered scenes."""
    from vhrharmonize.providers.standardized import StandardizedMetadata

    enriched: List[WorldViewScene] = []
    for scene in scenes:
        for image in scene.iter_images():
            if image.imd_file is None:
                continue
            worldview_metadata = load_worldview_metadata(image.imd_file, photo_basename=image.basename)
            image.worldview_metadata = worldview_metadata
            image.standardized_metadata = StandardizedMetadata.from_worldview_metadata(worldview_metadata)
        enriched.append(scene)
    return enriched


def load_worldview_scenes_from_tif_files(
    tif_files: Iterable[str],
    filter_basenames: Optional[Iterable[str]] = None,
) -> List[WorldViewScene]:
    """Discover WorldView scenes from tif files and attach parsed metadata."""
    scene_tree = discover_worldview_scene_tree_from_tif_files(
        tif_files,
        filter_basenames=filter_basenames,
    )
    return enrich_worldview_scenes_with_metadata(iter_worldview_scenes(scene_tree))


def find_files(root_folder_path: str, filter_basenames: Optional[Iterable[str]] = None) -> Dict[str, Dict[str, object]]:
    """Discover WorldView scenes under a folder root.
    Args:
        root_folder_path: Root folder to scan for tif files.
        filter_basenames: Optional basenames to keep.
    Returns:
        Flattened scene dictionary keyed by scene and catalog id.
    """
    tif_files: List[str] = []
    for root, _, files in os.walk(root_folder_path):
        for file_name in files:
            if os.path.splitext(file_name)[1].lower() != ".tif":
                continue
            tif_files.append(os.path.join(root, file_name))
    return {
        f"{scene.scene_id}_{scene.catalog_id}": scene.to_dict()
        for scene in load_worldview_scenes_from_tif_files(tif_files, filter_basenames=filter_basenames)
    }


__all__ = [
    "WorldViewFilenameParts",
    "WorldViewImage",
    "WorldViewMetadata",
    "WorldViewScene",
    "discover_worldview_scene_tree_from_tif_files",
    "enrich_worldview_scenes_with_metadata",
    "find_files",
    "iter_worldview_scenes",
    "load_worldview_metadata",
    "load_worldview_scenes_from_tif_files",
    "parse_worldview_basename",
    "parse_worldview_imd_file",
    "parse_worldview_imd_text",
]
