"""Persistent editable configuration files for the QGIS wrapper."""

from __future__ import annotations

import re
import shutil
from pathlib import Path


class ConfigManager:
    """Keep immutable defaults separate from user-edited plugin configs."""

    WORLDVIEW = "example.worldview.yml"
    HPC = "example.hpc.yml"
    SLURM = "example.slurm.sbatch"
    FILES = (WORLDVIEW, HPC, SLURM)

    def __init__(self, plugin_dir: Path):
        self.plugin_dir = Path(plugin_dir)
        self.default_dir = self.plugin_dir / "default_configs"
        self.config_dir = self.plugin_dir / "configs"

    def ensure_user_configs(self) -> None:
        """Create each editable config once, preserving later user changes."""
        self.config_dir.mkdir(parents=True, exist_ok=True)
        for filename in self.FILES:
            target = self.path(filename)
            if not target.exists():
                shutil.copy2(self.default_path(filename), target)

    def default_path(self, filename: str) -> Path:
        return self.default_dir / filename

    def path(self, filename: str) -> Path:
        return self.config_dir / filename

    def read(self, filename: str) -> str:
        return self.path(filename).read_text(encoding="utf-8")

    def write(self, filename: str, text: str) -> None:
        self.path(filename).write_text(text, encoding="utf-8")

    def reset(self, filename: str) -> str:
        """Replace an editable config with its immutable build-time default."""
        shutil.copy2(self.default_path(filename), self.path(filename))
        return self.read(filename)

    def staged_hpc_path(self) -> Path:
        """Resolve the configured staged HPC path without requiring valid YAML."""
        text = self.read(self.HPC)
        run_id = self._yaml_scalar(text, "run_id") or "1"
        template = self._yaml_scalar(text, "staged_hpc_file") or "configs/{run_id}.staged.hpc.yml"
        rendered = template.replace("{run_id}", run_id)
        path = Path(rendered).expanduser()
        return path if path.is_absolute() else self.plugin_dir / path

    @staticmethod
    def _yaml_scalar(text: str, key: str) -> str | None:
        match = re.search(rf"^\s*{re.escape(key)}\s*:\s*(.*?)\s*(?:#.*)?$", text, re.MULTILINE)
        if not match:
            return None
        value = match.group(1).strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in {'\"', "'"}:
            value = value[1:-1]
        return value or None
