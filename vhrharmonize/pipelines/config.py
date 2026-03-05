"""Pipeline config loading utilities."""

from typing import Any, Dict


def load_yaml_config(path: str) -> Dict[str, Any]:
    """Load a YAML config file and normalize hyphenated keys to underscores."""
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError("PyYAML is required to load pipeline config files.") from exc

    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    if not isinstance(cfg, dict):
        raise ValueError("YAML root must be a mapping")
    return {k.replace("-", "_"): v for k, v in cfg.items()}
