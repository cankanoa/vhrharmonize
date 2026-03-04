"""WorldView metadata wrappers."""

from typing import Dict, Tuple

from .files import get_metadata_from_files


def load_scene_metadata(root_file_path: str, imd_file: str, photo_basename: str) -> Tuple[Dict, Dict, Dict]:
    """Load scene/photo overrides and parsed IMD metadata."""
    return get_metadata_from_files(root_file_path, imd_file, photo_basename)
