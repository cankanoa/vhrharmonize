"""Shared CLI config loading helpers."""

from __future__ import annotations

from typing import Any, Dict


def load_yaml_config(path: str) -> Dict[str, Any]:
    """Load a YAML config file and normalize hyphenated keys to underscores."""
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError(
            "PyYAML is required for --config-yaml. Install it with `pip install pyyaml`."
        ) from exc

    with open(path, "r", encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle) or {}
    if not isinstance(loaded, dict):
        raise ValueError("Config YAML root must be a mapping/dictionary.")

    normalized: Dict[str, Any] = {}
    for key, value in loaded.items():
        normalized[str(key).replace("-", "_")] = value
    return normalized


__all__ = ["load_yaml_config"]
