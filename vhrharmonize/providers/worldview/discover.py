"""WorldView scene discovery wrappers."""

from typing import Dict, List

from .files import find_files, find_roots


def discover_scenes(input_dir: str, filter_basenames: List[str] | None = None) -> List[Dict]:
    """Discover WorldView scenes from Root files and associated imagery."""
    scenes: List[Dict] = []
    for root_folder_path, root_file_path in find_roots(input_dir):
        grouped = find_files(root_folder_path, root_file_path, filter_basenames or [])
        scenes.extend(grouped.values())
    return scenes
