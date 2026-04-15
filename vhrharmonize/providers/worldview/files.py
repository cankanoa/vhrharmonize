"""WorldView file discovery, grouping, and scene models."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Dict, Iterable, List, Mapping, Optional


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
        product = self.product_code.upper()
        if product.startswith("M"):
            return "mul"
        if product.startswith("P"):
            return "pan"
        return None

    @property
    def acquisition_datetime_utc(self) -> datetime:
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
        return self.filename_parts.basename

    @property
    def image_role(self) -> Optional[str]:
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
        return self.mul_image or self.pan_image

    @property
    def primary_basename(self) -> Optional[str]:
        image = self.primary_image
        return image.basename if image else None

    def get_image(self, role: str) -> Optional[WorldViewImage]:
        role_norm = role.strip().lower()
        if role_norm == "mul":
            return self.mul_image
        if role_norm == "pan":
            return self.pan_image
        return None

    def set_image(self, image: WorldViewImage) -> None:
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
        return [image for image in (self.mul_image, self.pan_image) if image is not None]

    def to_dict(self) -> Dict[str, object]:
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
    directory = os.path.dirname(path)
    base = os.path.basename(directory).upper()
    if base.endswith("_MUL") or base.endswith("_PAN"):
        return os.path.dirname(directory)
    return directory


def _companion_files_for_stem(directory: str, basename: str) -> Dict[str, Optional[str]]:
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
        elif ext == ".shp":
            companions["shp_file"] = candidate_path
        elif ext == ".til":
            companions["til_file"] = candidate_path
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

        companions = _companion_files_for_stem(os.path.dirname(tif_file), filename_parts.basename)
        scene_root = _scene_root_for_path(tif_file)
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


def enrich_worldview_scenes_with_metadata(scenes: Iterable[WorldViewScene]) -> List[WorldViewScene]:
    """Attach provider-specific and standardized metadata to discovered scenes."""
    from vhrharmonize.providers.standardized import StandardizedMetadata
    from vhrharmonize.providers.worldview.metadata import load_worldview_metadata

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


def find_files(root_folder_path, filter_basenames=None):
    """Backward-compatible flattened discovery wrapper for a folder root."""
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
    "WorldViewScene",
    "discover_worldview_scene_tree_from_tif_files",
    "enrich_worldview_scenes_with_metadata",
    "find_files",
    "iter_worldview_scenes",
    "load_worldview_scenes_from_tif_files",
    "parse_worldview_basename",
]
