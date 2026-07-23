"""Shared CLI config loading helpers."""

from __future__ import annotations

from typing import Any, Dict

import yaml


def load_yaml_config(path: str) -> Dict[str, Any]:
    """Load and normalize a YAML config file.
    Args:
        path: YAML config file path.
    Returns:
        Loaded config mapping with hyphenated keys normalized to underscores.
    """
    with open(path, "r", encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle) or {}
    if not isinstance(loaded, dict):
        raise ValueError("Config YAML root must be a mapping/dictionary.")

    normalized: Dict[str, Any] = {}
    for key, value in loaded.items():
        normalized[str(key).replace("-", "_")] = value
    return normalized


__all__ = ["load_yaml_config"]
