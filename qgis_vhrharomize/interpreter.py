"""Locate the Python executable which matches the running QGIS process."""

from __future__ import annotations

import platform
import sys
from pathlib import Path

from qgis.core import Qgis, QgsMessageLog


def _log(message: str) -> None:
    QgsMessageLog.logMessage(message, "VHRHarmonize", level=Qgis.MessageLevel.Info)


def python_command() -> str:
    """Return the Python command associated with the running QGIS process."""
    if (Path(sys.prefix) / "conda-meta").exists():
        _log(f"Using the Conda environment's Python at {sys.executable}")
        return sys.executable

    if platform.system() == "Windows":
        for filename in ("python.exe", "python3.exe"):
            path = Path(sys.prefix) / filename
            if path.is_file():
                _log(f"Using Windows Python at {path}")
                return str(path)
        _log(f"Using Windows fallback at {sys.executable}")
        return sys.executable

    if platform.system() == "Darwin":
        base_paths = (Path(sys.prefix), Path(sys.prefix) / "bin", Path(sys.executable).parent)
        for base_path in base_paths:
            for filename in ("python", "python3"):
                path = base_path / filename
                if path.is_file():
                    _log(f"Using macOS Python at {path}")
                    return str(path)
        _log(f"Using macOS fallback at {sys.executable}")
        return sys.executable

    _log(f"Using Python at {sys.executable}")
    return sys.executable
