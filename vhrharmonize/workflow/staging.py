"""Materialize declared plugin file arguments for remote Slurm execution."""

from copy import deepcopy
import hashlib
import os
from pathlib import Path
from vhrharmonize.plugins.base import INPUT_PATH_FEATURES, OUTPUT_PATH_FEATURES
from .config import step_settings, shared_blocks
from .engine import Workflow, _within, _required_arguments
from .values import contains_pending, lookup, evaluate_settings, resolve, Deferred, Pending, path
from .paths import directory_bindings


def _hash(value):
    return hashlib.sha256(value.encode()).hexdigest()[:12]


def _literal(value):
    """Preserve strings that might otherwise be interpreted as typed references."""
    if isinstance(value, str):
        return "literal:" + value
    if isinstance(value, list):
        return [_literal(v) for v in value]
    if isinstance(value, dict):
        return {k: _literal(v) for k, v in value.items()}
    return value


def stage_workflow(config, *, config_dir, remote_output_dir, remote_temp_dir, remote_reference_dir):
    workflow = Workflow(config, config_dir=config_dir).plan()
    if workflow.barrier_index is not None:
        raise ValueError(
            "HPC staging requires scenes and file paths known during planning; "
            "its required upstream processing must finish before these scenes can be staged"
        )
    staged = deepcopy(workflow.config)
    shared_names = [name for name, settings in staged.items() if settings.get("plugin") == "shared"]
    shared_name = shared_names[0] if shared_names else "shared"
    while shared_name in staged and shared_name not in shared_names:
        shared_name = "_" + shared_name
    merged_shared = {key: value for block in shared_blocks(staged) for key, value in block.items()}
    staged = {name: block for name, block in staged.items() if name not in shared_names}
    staged = {shared_name: {**merged_shared, "plugin": "shared", "core:run": True}, **staged}
    records = deepcopy(workflow.initial_records)
    path_map, downloads = {}, {}
    root_map = {}
    contexts = [
        workflow.initial_context,
        *[record["context"] for record in records],
        *[node.context for node in workflow.nodes],
    ]
    for context in contexts:
        for role, selectors in workflow.directory_locations.items():
            remote = remote_temp_dir if role == "temp_dir" else remote_output_dir
            for selector, local in directory_bindings(context, {role: selectors}).items():
                root_map[local] = (
                    os.path.join(remote, "scenes", _hash(local))
                    if selector.startswith("var.")
                    else remote
                )
    for node in workflow.nodes:
        output_staging = set(node.paths("output_hpc_staging_paths"))
        output_downloads = set(node.paths("output_hpc_download_paths"))
        for filename in sorted(output_staging | output_downloads):
            directories = node.directories
            temporary = _within(filename, directories["temp_dir"])
            roots = directories["temp_dir"] if temporary else directories["output_dir"]
            local_root = next(
                (root for root in sorted(roots, key=len, reverse=True) if _within(filename, root)),
                None,
            )
            remote_root = root_map.get(
                local_root, remote_temp_dir if temporary else remote_output_dir
            )
            relative = (
                os.path.relpath(filename, local_root)
                if local_root
                else f"{_hash(str(Path(filename).parent))}/{Path(filename).name}"
            )
            path_map[filename] = os.path.join(remote_root, relative)
            if not temporary and filename in output_downloads:
                downloads[filename] = path_map[filename]
        checkpoint_outputs = node.paths("output_context_checkpoint_paths")
        anchor = checkpoint_outputs[0] if checkpoint_outputs else None
        if node.checkpoint and anchor in path_map:
            path_map[node.checkpoint] = path_map[anchor] + ".context.json"
            if (
                node.dynamic_names
                and anchor in output_downloads
                and not _within(node.checkpoint, node.directories["temp_dir"])
            ):
                downloads[node.checkpoint] = path_map[node.checkpoint]
    required = set()

    def constants(node):
        bindings = workflow.constant_steps.get(node.step_index)
        return bindings.planned if bindings else None

    for node in workflow.nodes:
        if node.status == "loaded":
            required.update(node.demanded_paths & set(node.paths("output_hpc_staging_paths")))
            if (
                node.dynamic_names
                and node.checkpoint
                and os.path.isfile(node.checkpoint)
                and node.paths("output_context_checkpoint_paths")[0]
                in node.paths("output_hpc_staging_paths")
            ):
                required.add(node.checkpoint)
        elif node.status == "processing":
            required.update(
                p
                for p in [*node.paths("input_hpc_staging_paths"), *node.requirements]
                if p not in path_map
            )
        params = evaluate_settings(
            node.step["settings"],
            node.pre_context,
            records=node.collection_snapshot if node.record is None else None,
            planning=True,
            constants=constants(node),
        )[0]
        for name in node.file_features["input_hpc_staging_paths"] - params.keys():
            if name in workflow.shared:
                try:
                    params[name] = resolve(workflow.shared[name], node.context)
                except Deferred:
                    params[name] = Pending(name)
        unresolved = [
            name
            for name in node.file_features["input_hpc_staging_paths"]
            if contains_pending(params.get(name))
        ]
        if unresolved:
            raise ValueError(
                f"HPC file arguments must resolve during planning: {sorted(unresolved)}"
            )
    for filename in sorted(required):
        path_map.setdefault(
            filename,
            os.path.join(
                remote_reference_dir,
                "inputs",
                _hash(str(Path(filename).parent)),
                Path(filename).name,
            ),
        )
    missing = [filename for filename in required if not os.path.exists(filename)]
    if missing:
        raise FileNotFoundError(f"Required workflow inputs are missing: {missing}")
    uploads = {filename: path_map[filename] for filename in sorted(required)}

    def rewrite(value):
        if isinstance(value, str):
            if value in path_map:
                return path_map[value]
            if value in root_map:
                return root_map[value]
            return value
        if isinstance(value, list):
            return [rewrite(v) for v in value]
        if isinstance(value, dict):
            return {k: rewrite(v) for k, v in value.items()}
        return value

    for record in records:
        record["context"] = rewrite(record["context"])
        record["source_paths"] = rewrite(record["source_paths"])

    def materialize(node, key, value, *, file_value=False, transfer=True):
        settings = step_settings(staged, node.step)
        value = rewrite(value) if transfer else value
        if node.record is None:
            settings[key] = _literal(value)
        else:
            name = f"staged_{node.step_index}_{key.replace(':','_').replace('.','_')}"
            records[node.record]["context"]["var"][name] = value
            settings[key] = "var:" + name
            if file_value:
                records[node.record].setdefault("path_fields", []).append(name)

    # Scene-independent constants stay workflow-wide, including for empty scenes.
    # Never replace their right-hand sides with per-scene staged var: references.
    for step_index, bindings in workflow.constant_steps.items():
        settings = step_settings(staged, workflow.steps[step_index])
        for name, value in bindings.planned.items():
            if not contains_pending(value):
                settings[name.replace(".", ":", 1)] = _literal(rewrite(value))

    for node in workflow.nodes:
        # Freeze pre-call assignments at their local planning values. This also
        # preserves suffix updates and file aliases after moving to another host.
        for name, value in node.updates.items():
            if not contains_pending(value):
                value = lookup(node.context, name)
                materialize(node, name.replace(".", ":", 1), value)
        for name, value in node.file_arguments(*INPUT_PATH_FEATURES, *OUTPUT_PATH_FEATURES).items():
            transfer = name in (
                node.file_features["input_hpc_staging_paths"]
                | node.file_features["output_hpc_staging_paths"]
                | node.file_features["output_hpc_download_paths"]
            )
            materialize(node, "param:" + name, value, file_value=transfer, transfer=transfer)
        if node.requirements:
            params = evaluate_settings(
                node.step["settings"],
                node.pre_context,
                records=node.collection_snapshot if node.record is None else None,
                planning=True,
                constants=constants(node),
            )[0]
            normalized = _required_arguments(
                params,
                node.step,
                node.context,
                node.requirements,
                records=node.collection_snapshot if node.record is None else None,
            )
            for name, value in normalized.items():
                if name not in node.file_arguments(
                    *INPUT_PATH_FEATURES, *OUTPUT_PATH_FEATURES
                ) and not contains_pending(value):
                    if rewrite(value) != params[name]:
                        materialize(node, "param:" + name, value)
            materialize(node, "core:requires", node.requirements, file_value=True)
    for record in records:
        record["file_paths"] = list(set(path_map.values()) | {remote_output_dir, remote_temp_dir})
        # IDs need not be paths. Rebase only IDs that are actual staged file paths.
        record["id"] = path_map.get(record["id"], record["id"])
    if workflow.controls["output_metadata_path"] is not None:
        metadata_targets = {}
        final_contexts = []
        for node in workflow.final_nodes():
            if node.record is not None or not records:
                final_contexts.append((node.record, node.context, node.collection_snapshot))
            else:
                final_contexts.extend(
                    (index, {**context, "const": node.context["const"]}, node.collection_snapshot)
                    for index, context in enumerate(node.collection_snapshot)
                )
        if not workflow.nodes:
            final_contexts = [
                (index, record["context"], None) for index, record in enumerate(workflow.records)
            ]
        for record_index, context, collection in final_contexts:
            local = path(
                resolve(workflow.controls["output_metadata_path"], context, records=collection),
                base_dir=config_dir,
            )
            root = next(
                (root for root in sorted(root_map, key=len, reverse=True) if _within(local, root)),
                None,
            )
            remote = (
                os.path.join(root_map[root], os.path.relpath(local, root))
                if root
                else os.path.join(remote_output_dir, "metadata", _hash(local), Path(local).name)
            )
            downloads[local] = remote
            path_map[local] = remote
            if not workflow.controls["delete_final_json_first"] and os.path.isfile(local):
                uploads[local] = remote
            if record_index is None:
                metadata_targets[None] = remote
            else:
                records[record_index]["context"]["var"]["staged_final_json_path"] = remote
                records[record_index].setdefault("path_fields", []).append("staged_final_json_path")
        staged.setdefault(shared_name, {})["core:output_metadata_path"] = (
            _literal(metadata_targets[None])
            if None in metadata_targets
            else "var:staged_final_json_path"
        )
    staged = rewrite(staged)

    # Preserve a single workflow-wide constant scope, even with no imported scenes.
    shared = {
        k: v for k, v in staged.get(shared_name, {}).items() if not k.startswith(("const:", "var:"))
    }
    shared.update(
        {"const:" + k: _literal(rewrite(v)) for k, v in workflow.initial_context["const"].items()}
    )
    staged = {shared_name: shared, **{k: v for k, v in staged.items() if k != shared_name}}
    # Replace locally evaluated discovery with an ordinary registered snapshot
    # plugin. Remote execution neither repeats globs nor imports local documents.
    for index in workflow.preflight_steps:
        step_settings(staged, workflow.steps[index])["core:run"] = False
    restore_name = "restore_scenes"
    while restore_name in staged:
        restore_name = "_" + restore_name
    staged = {
        shared_name: staged.pop(shared_name),
        **(
            {
                restore_name: {
                    "plugin": "restore_scenes",
                    "core:run": True,
                    "param:records": _literal(records),
                    "param:directory_locations": workflow.directory_locations,
                }
            }
            if "var" in workflow.initial_context
            else {}
        ),
        **staged,
    }
    return staged, uploads, downloads
