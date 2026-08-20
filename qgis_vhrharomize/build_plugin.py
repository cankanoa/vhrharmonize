"""Build the VHRHarmonize QGIS plugin from repository sources."""

import json
import re
import shutil
import tomllib
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
PLUGIN_DIR = ROOT / "qgis_vhrharomize"
ZIP_PATH = ROOT / "qgis_vhrharomize.zip"


def generate_optional_dependency_commands(
    input_toml_path: Path,
    output_json_path: Path,
) -> None:
    """Write install templates for the base and optional dependencies."""
    with input_toml_path.open("rb") as handle:
        pyproject = tomllib.load(handle)
    project = pyproject["project"]
    base_dependencies = list(project.get("dependencies", []))
    optional_dependencies = project.get("optional-dependencies", {})
    ordered_names = sorted(optional_dependencies, key=lambda name: name != "defaults")
    groups = []
    for name in ordered_names:
        dependencies = optional_dependencies[name]
        groups.append(
            {
                "name": name,
                "dependencies": list(dependencies),
                "command": [
                    "python",
                    "-u",
                    "-m",
                    "pip",
                    "install",
                    *dependencies,
                    "--target",
                    "{target}",
                    "--upgrade",
                ],
            }
        )
    output_json_path.write_text(
        json.dumps(
            {
                "version": 1,
                "base": {
                    "dependencies": base_dependencies,
                    "command": [
                        "python",
                        "-u",
                        "-m",
                        "pip",
                        "install",
                        *base_dependencies,
                        "--target",
                        "{target}",
                        "--upgrade",
                    ],
                },
                "groups": groups,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def copy_vhrharmonize_package() -> None:
    target = PLUGIN_DIR / "vhrharmonize"
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(
        ROOT / "vhrharmonize",
        target,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", ".DS_Store"),
    )


def copy_default_configs() -> None:
    target = PLUGIN_DIR / "default_configs"
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(ROOT / "configs", target, ignore=shutil.ignore_patterns("__pycache__", ".DS_Store"))
    runtime_configs = PLUGIN_DIR / "configs"
    if runtime_configs.exists():
        shutil.rmtree(runtime_configs)


def read_metadata_version() -> str:
    metadata = (PLUGIN_DIR / "metadata.txt").read_text(encoding="utf-8")
    match = re.search(r"(?m)^version=(.+)$", metadata)
    if match is None:
        raise ValueError("QGIS metadata.txt does not contain a version")
    return match.group(1).strip()


def validate_metadata_version() -> None:
    with (ROOT / "pyproject.toml").open("rb") as handle:
        package_version = tomllib.load(handle)["project"]["version"]
    metadata_version = read_metadata_version()
    if metadata_version != package_version:
        raise ValueError(
            "QGIS metadata version "
            f"{metadata_version} does not match package version {package_version}"
        )


def build_zip() -> None:
    if ZIP_PATH.exists():
        ZIP_PATH.unlink()
    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(PLUGIN_DIR.rglob("*")):
            relative = path.relative_to(ROOT)
            if path.is_dir() or "__pycache__" in path.parts or path.suffix in {".pyc", ".pyo"}:
                continue
            archive.write(path, relative.as_posix())


def build() -> Path:
    validate_metadata_version()
    copy_vhrharmonize_package()
    copy_default_configs()
    shutil.copy2(ROOT / "LICENSE", PLUGIN_DIR / "LICENSE")
    legacy_requirements = PLUGIN_DIR / "requirements.txt"
    if legacy_requirements.exists():
        legacy_requirements.unlink()
    generate_optional_dependency_commands(
        ROOT / "pyproject.toml",
        PLUGIN_DIR / "optional_dependency_commands.json",
    )
    build_zip()
    return ZIP_PATH


if __name__ == "__main__":
    output = build()
    print(f"Built {output}")
