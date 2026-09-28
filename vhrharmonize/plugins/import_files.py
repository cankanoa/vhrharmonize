"""Discover files and build plain scene dictionaries from explicit metadata rules."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

from wcmatch import glob
from vhrharmonize.io.metadata import read_metadata
from vhrharmonize.workflow.values import path, resolve
from .base import FunctionPlugin

FLAGS = glob.GLOBSTAR | glob.BRACE | glob.EXTGLOB | glob.GLOBTILDE | glob.NEGATE
_RESERVED_FIELDS = {"scene_id", "file_path", "source_paths", "temp_dir", "output_dir"}


def _metadata_rules(create_metadata_json):
    if create_metadata_json is None:
        return {}
    if not isinstance(create_metadata_json, dict):
        raise ValueError("create_metadata_json must be a name-to-rule mapping")
    for name, rule in create_metadata_json.items():
        if not isinstance(name, str) or not name:
            raise ValueError("create_metadata_json field names must be nonempty strings")
        if name in _RESERVED_FIELDS:
            raise ValueError(f"create_metadata_json cannot replace reserved field {name!r}")
        if not isinstance(rule, dict) or set(rule) not in ({"path"}, {"to_json"}):
            raise ValueError(
                f"create_metadata_json.{name} must contain exactly one of path or to_json"
            )
    return create_metadata_json


def _create_metadata(item, rules):
    context = {"var": item}
    base_dir = str(Path(item["file_path"]).parent)
    for name, rule in rules.items():
        kind, template = next(iter(rule.items()))
        value = path(resolve(template, context, returned=item), base_dir=base_dir)
        if kind == "to_json":
            if not isinstance(value, str):
                raise ValueError(f"create_metadata_json.{name}.to_json must resolve to one path")
            item[name] = read_metadata(value)
            item["source_paths"].append(value)
        else:
            matches = sorted({p for p in glob.glob(value, flags=FLAGS) if os.path.isfile(p)})
            item[name] = matches
            item["source_paths"].extend(matches)
    item["source_paths"] = list(dict.fromkeys(item["source_paths"]))
    return item


def _scene_id(value):
    if not isinstance(value, str) or not value.strip():
        raise ValueError("scene_id must resolve to a nonempty string")
    return value


def import_file(
    file_path: str, *, scene_id: str = "var:file_path", create_metadata_json: dict | None = None
) -> dict:
    """Return one file's configured scene ID, path, metadata and protected source paths.

    Args:
        file_path: Existing file to import.
        scene_id: Nonempty string or per-file reference/expression; default var:file_path.
            Evaluated after metadata rules, before YAML scene mappings.
        create_metadata_json: Named rules with either path (companion glob, returning
            a sorted path list) or to_json (metadata file, returning its decoded data).
            Relative paths use this file's parent. Expressions read var.file_path
            and fields produced by earlier rules. No JSON file is written.
    """
    _scene_id(scene_id)
    rules = _metadata_rules(create_metadata_json)
    filename = Path(file_path).expanduser().resolve()
    if not filename.is_file():
        raise FileNotFoundError(filename)
    item = _create_metadata({"file_path": str(filename), "source_paths": [str(filename)]}, rules)
    item["scene_id"] = _scene_id(resolve(scene_id, {"var": item}, returned=item))
    return item


def import_files(
    search_glob,
    *,
    scene_id: str = "var:file_path",
    create_metadata_json: dict | None = None,
    temp_dir: str = "sys",
    output_dir: str = "./output",
    temp_dir_scope: str = "const",
    output_dir_scope: str = "var",
    where=True,
) -> dict:
    """Import matching files and publish directory roots and plain scene objects.

    Relative paths resolve from each discovered file. Constant directory roots
    use the first file's directory (cwd when there are no matches). Use absolute
    roots for a shared project directory. Each directory has its own scope:
    temp_dir_scope defaults to const and output_dir_scope defaults to var.
    The search glob itself is relative to the Python working directory.

    Per-file metadata and filter expressions may be passed as literal:expr:... in
    YAML. They read var.file_path and imported fields, before YAML var assignments.
    Constants used as arguments are resolved by the caller; the function receives
    no workflow context. Decoded metadata stays in memory; core's
    output_metadata_path setting saves the final processing context separately.

    Args:
        search_glob: Pattern or list of patterns; relative to the Python working directory.
        scene_id: Nonempty string or per-file reference/expression; default var:file_path.
            Resolved after metadata rules and returned in each scene as scene_id.
            Later imports merge records with the same ID; IDs must be unique within this import.
        create_metadata_json: Mapping of field names to {path: glob} or {to_json: file}.
            Fields are added directly to each scene. Rules run in declaration order.
        temp_dir: Temporary directory path; sys creates a system temporary directory.
        output_dir: Persistent output directory, relative to the matched file.
        temp_dir_scope: const (default) shares one temporary root; var uses per-file roots.
        output_dir_scope: var (default) uses per-file output roots; const shares one root.
        where: Boolean or per-file expression controlling which scenes are imported.
    """
    _scene_id(scene_id)
    if not search_glob:
        raise ValueError("import_files.param:search_glob is required")
    for name, scope in (("temp_dir_scope", temp_dir_scope), ("output_dir_scope", output_dir_scope)):
        if scope not in {"const", "var"}:
            raise ValueError(f"{name} must be const or var")
    rules = _metadata_rules(create_metadata_json)
    patterns = search_glob if isinstance(search_glob, list) else [search_glob]
    files = sorted(
        {
            os.path.abspath(p)
            for pattern in patterns
            for p in glob.glob(os.path.expanduser(pattern), flags=FLAGS)
            if os.path.isfile(p)
        }
    )

    def directories(context, base_dir, scope):
        result = {}
        for name, template, selected_scope in (
            ("temp_dir", temp_dir, temp_dir_scope),
            ("output_dir", output_dir, output_dir_scope),
        ):
            if selected_scope != scope:
                continue
            value = resolve(template, context)
            result[name] = (
                tempfile.mkdtemp(prefix="vhr-")
                if name == "temp_dir" and value == "sys"
                else path(value, base_dir=base_dir)
            )
        return result

    base_dir = str(Path(files[0]).parent) if files else os.getcwd()
    published = directories({}, base_dir, "const")
    scenes = []
    seen_ids = set()
    for filename in files:
        item = {"file_path": filename, "source_paths": [filename]}
        context = {"var": item}
        item.update(directories(context, str(Path(filename).parent), "var"))
        _create_metadata(item, rules)
        item["scene_id"] = _scene_id(resolve(scene_id, context, returned=item))
        if resolve(where, context, returned=item):
            if item["scene_id"] in seen_ids:
                raise ValueError(
                    f"Scene identifiers must be unique within an import: {item['scene_id']!r}"
                )
            seen_ids.add(item["scene_id"])
            scenes.append(item)
    return {"const": published, "scenes": scenes}


class ImportFiles(FunctionPlugin):
    """Add imported scenes and missing fields while retaining earlier workflow state."""

    target = "vhrharmonize.plugins.import_files:import_files"
    scene_records_return = "scenes"
    scene_records_mode = "merge"
    constant_values_return = "const"
    scene_id_return = "scene_id"  # The function computes this from its scene_id argument.
    source_file_protection_paths_return = "source_paths"
    temporary_directory_context_paths = ("var.temp_dir", "const.temp_dir")
    output_directory_context_paths = ("var.output_dir", "const.output_dir")
