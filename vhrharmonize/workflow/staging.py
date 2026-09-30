"""Stage required files by replacing declared roots in an ordinary workflow recipe."""

from copy import deepcopy
import os
from pathlib import Path
import tempfile

from .context_io import CONTEXT_CONTROLS, SAVE_CONTROLS, validate_snapshot
from .engine import Workflow, _paths, _within
from .paths import directory_bindings, validate_path_mappings
from .registry import load_plugin
from .staging_files import context_path_value, remote_directory, selected_values, rewrite_file_assignments
from .values import Deferred, lookup, path, remap_paths, contains_pending
from vhrharmonize.io.metadata import write_json


def stage_workflow(config, *, config_dir, remote_work_dir=None, path_mappings=None,
                   remote_output_dir=None, remote_temp_dir=None, remote_reference_dir=None,
                   context_staging_dir=None, upload_groups=None):
    """Map directory roots or flatten selected files, and transfer dependencies.

    Modern staging requires mapping coverage for workflow data and context files.
    Legacy remote directory arguments remain accepted for existing callers.
    """
    modern = remote_work_dir is not None
    if modern:
        remote_output_dir = os.path.join(remote_work_dir, "products")
        remote_temp_dir = os.path.join(remote_work_dir, "files")
        remote_reference_dir = os.path.join(remote_work_dir, "inputs")
    elif not all((remote_output_dir, remote_temp_dir, remote_reference_dir)):
        raise ValueError("HPC staging requires remote_work_dir")
    mappings = validate_path_mappings({} if path_mappings is None else path_mappings)
    workflow = Workflow(config, config_dir=config_dir).plan()
    if workflow.barrier_index is not None:
        raise ValueError("HPC staging requires known scenes and paths; use run_to_step_before_prepare to produce them first")
    staged = deepcopy(workflow.config)
    context_groups = [[workflow.initial_context, *[record["context"] for record in workflow.initial_records]]]
    for index in dict.fromkeys(node.step_index for node in workflow.nodes):
        nodes = [node for node in workflow.nodes if node.step_index == index]
        context_groups.extend((
            [context for node in nodes for context in [node.pre_context, *node.collection_snapshot]],
            [context for node in nodes for context in [node.context, *node.collection_result]],
        ))
    contexts = [context for group in context_groups for context in group]
    roots, root_selectors, file_bindings, files = {}, {}, {}, set()
    bindings = {}
    declared_paths = set(workflow.discovery_sources)
    for node in workflow.nodes:
        declared_paths.update(node.paths("input_hpc_staging_paths", "output_hpc_staging_paths"))
        declared_paths.update(node.requirements)
    directory_selectors = {selector.replace(".", ":", 1)
                           for selectors in workflow.directory_locations.values() for selector in selectors}

    def add_root(local, remote):
        local = path(local, base_dir=config_dir)
        if local in roots and roots[local] != remote:
            raise ValueError(f"Conflicting HPC path mappings for {local}")
        roots[local] = remote.rstrip("/") or "/"

    for selector, template in mappings.items():
        selected = selected_values(selector, context_groups)
        bindings[selector] = selected
        file_binding = False
        for context, value in selected:
            remote = remote_directory(template, context, selector)
            values = value if isinstance(value, list) else [value]
            file_binding |= isinstance(value, list)  # Preserve lists, including empty companions.
            for local in values:
                local = path(local, base_dir=config_dir)
                is_file = os.path.isfile(local) or (
                    not os.path.isdir(local) and selector not in directory_selectors and local in declared_paths)
                destination = os.path.join(remote, os.path.basename(local)) if is_file else remote
                add_root(local, destination)
                root_selectors.setdefault(local, selector)
                if is_file:
                    files.add(local)
                    file_binding = True
        if file_binding:
            file_bindings[selector] = selected

    if not modern:
        # Compatibility for the earlier three-directory API, with no path hashes.
        roots.setdefault(os.path.abspath(config_dir), remote_reference_dir)
        for context in contexts:
            for role, selectors in workflow.directory_locations.items():
                for _, local in directory_bindings(context, {role: selectors}).items():
                    roots[local] = remote_temp_dir if role == "temp_dir" else remote_output_dir

    def mapped(filename, owner="workflow"):
        local = path(filename, base_dir=config_dir)
        root = next((root for root in sorted(roots, key=len, reverse=True) if _within(local, root)), None)
        if root is None:
            raise ValueError(f"{owner}: required path has no HPC root mapping or explicit file mapping: {local}")
        return roots[root] if local == root else os.path.join(roots[root], os.path.relpath(local, root))

    needed, downloads, destinations = set(), {}, {}
    sources, products = {}, {}
    generated = {filename for node in workflow.nodes if node.status == "processing" for filename in node.paths("output_hpc_staging_paths")}

    def register(filename, owner):
        remote = mapped(filename, owner)
        previous = destinations.get(remote)
        if previous is not None and previous != filename:
            raise ValueError(f"HPC path collision: {previous} and {filename} both map to {remote}")
        destinations[remote] = filename
        return remote

    for node in workflow.nodes:
        if not node.needed:
            continue
        owner = node.step["name"]
        for feature, labels in (("input_hpc_staging_paths", sources), ("output_hpc_staging_paths", products)):
            for name, value in node.file_arguments(feature).items():
                reference = node.step["settings"].get("param:" + name, "param:" + name)
                if not isinstance(reference, str) or not reference.startswith(("var:", "const:")):
                    reference = "param:" + name
                for filename in _paths(value):
                    labels.setdefault(filename, (owner, reference))
        for name, value in node.params.items():
            if name in node.file_features["input_hpc_staging_paths"] | node.file_features["output_hpc_staging_paths"] and contains_pending(value):
                raise ValueError(f"{owner}: HPC file parameter {name} must resolve during planning")
        for filename in node.paths("output_hpc_staging_paths", "output_hpc_download_paths"):
            remote = register(filename, owner)
            if filename in node.paths("output_hpc_download_paths") and (filename in node.paths("output_target_paths") if modern else not _within(filename, node.directories["temp_dir"])):
                downloads[filename] = remote
        if node.status == "loaded":
            needed.update(node.demanded_paths & set(node.paths("output_hpc_staging_paths")))
        else:
            for filename in [*node.paths("input_hpc_staging_paths"), *node.requirements]:
                register(filename, owner)
                if filename not in generated:
                    needed.add(filename)

    # Discovery is rerun remotely unless the recipe explicitly loads its context.
    # The importer records the actual files it read, including metadata/companions.
    needed.update(workflow.discovery_sources)
    needed.update(workflow.context_files.read_paths)
    discovery_labels = {name: load_plugin(staged[name]["plugin"]).discovery_input_parameter
                        for name in set(workflow.discovery_source_steps.values())}
    sources.update({filename: (step, "param:" + discovery_labels[step] if discovery_labels[step] else "discovery")
                    for filename, step in workflow.discovery_source_steps.items()})
    sources.update(products)
    sources.update({filename: (step, "core:load_context") for filename, step in workflow.context_files.read_steps.items()})

    # Explicit context files are normal declared data inputs/outputs, never hidden sidecars.
    context_inputs = set(workflow.context_files.read_paths)
    context_outputs = set(workflow.context_files.write_paths)
    step_contexts = {step["name"]: [] for step in workflow.steps}
    for node in workflow.nodes:
        if node.needed:
            step_contexts[node.step["name"]].append(node.context)
    for index in workflow.preflight_steps:
        step_contexts[workflow.steps[index]["name"]].extend(record["context"] for record in workflow.initial_records)
    for index, step in enumerate(workflow.steps):
        if not step["run"] or not step_contexts[step["name"]] and index not in workflow.preflight_steps:
            continue
        for context in step_contexts[step["name"]] or [workflow.initial_context]:
            for control in CONTEXT_CONTROLS:
                for template in step.get(control, {}):
                    filename = workflow.context_files.filename(template, context)
                    register(filename, f"{step['name']} core:{control}")
                    if control in SAVE_CONTROLS:
                        context_outputs.add(filename)
                    elif os.path.isfile(filename):
                        context_inputs.add(filename)
                        needed.add(filename)
    for filename in context_outputs:
        downloads[filename] = register(filename, "saved context")

    shared_names = [name for name, block in staged.items() if block.get("plugin") == "shared" and block.get("core:run", False)]
    if not shared_names:
        name = "hpc_settings"
        while name in staged:
            name = "_" + name
        staged[name] = {"plugin": "shared", "core:run": True}
        shared_names = [name]

    statistics_updates = {}
    for control in ("save_statistics_path", "load_statistics_path"):
        value = workflow.controls[control]
        if value is None:
            continue
        filename = path(workflow._resolve(value, workflow.initial_context), base_dir=config_dir)
        # Statistics are generated control artifacts with a default location, not scene data.
        remote = os.path.join(remote_output_dir if control.startswith("save") else remote_reference_dir, "statistics", Path(filename).name)
        if control.startswith("save"):
            downloads[filename] = remote
        elif os.path.isfile(filename):
            # Separate history input so an upload cannot overwrite remote appended results.
            remote = os.path.join(remote_reference_dir, "statistics", Path(filename).name)
            destinations[remote] = filename
        statistics_updates["core:" + control] = remote

    if workflow.controls["output_metadata_path"] is not None:
        for context in [record["context"] for record in workflow.records] or [workflow.initial_context]:
            try:
                filename = path(workflow._resolve(workflow.controls["output_metadata_path"], context), base_dir=config_dir)
            except Deferred:
                continue
            downloads[filename] = register(filename, "core:output_metadata_path")
            if not workflow.controls["delete_final_json_first"] and os.path.isfile(filename):
                needed.add(filename)

    uploads = {}
    for filename in sorted(needed):
        remote = register(filename, "input")
        if not os.path.exists(filename):
            raise FileNotFoundError(f"Required workflow input is missing: {filename}")
        uploads[filename] = remote
    history = workflow.controls["load_statistics_path"]
    if history is not None:
        local = path(workflow._resolve(history, workflow.initial_context), base_dir=config_dir)
        if os.path.isfile(local):
            uploads[local] = os.path.join(remote_reference_dir, "statistics", Path(local).name)
            sources[local] = ("statistics", "core:load_statistics_path")

    def rewrite_scalar(value):
        if not isinstance(value, str):
            return value
        for prefix in ("path:", "literal:", ""):
            if not value.startswith(prefix):
                continue
            text = value[len(prefix):]
            if text.startswith(("/", "~/", "./", "../")):
                normalized = path(text, base_dir=config_dir)
                transformed = remap_paths(normalized, roots)
                if transformed != normalized:
                    return (prefix or "path:") + transformed
                return value
        if value.startswith(("expr:", "literal:expr:")):
            # Only literal path prefixes inside expressions change; variable references stay intact.
            for local in sorted(roots, key=len, reverse=True):
                for quote in ("'", '"'):
                    value = value.replace(quote + local + "/", quote + roots[local].rstrip("/") + "/")
                    value = value.replace(quote + local + quote, quote + roots[local] + quote)
        return value

    def rewrite(value):
        if isinstance(value, dict):
            return {key: rewrite(item) for key, item in value.items()}
        if isinstance(value, list):
            return [rewrite(item) for item in value]
        return rewrite_scalar(value)

    staged = rewrite(staged)
    staged[shared_names[-1]].update(statistics_updates)
    metadata = workflow.controls["output_metadata_path"]
    if metadata is not None and not metadata.startswith(("expr:", "var:", "const:", "path:expr:", "path:var:", "path:const:")):
        filename = path(workflow._resolve(metadata, workflow.initial_context), base_dir=config_dir)
        staged[shared_names[-1]]["core:output_metadata_path"] = mapped(filename, "core:output_metadata_path")
    # Directory bindings replace their first definition. Files keep their names
    # and list shapes, and only need an assignment edit if normal evaluation
    # against the relocated import would produce a different value.
    for selector in mappings:
        if selector in file_bindings:
            continue
        pairs = [(context,
                  remap_paths(path(value, base_dir=config_dir), roots)) for context, value in bindings[selector]]
        assignment = context_path_value(pairs, workflow, roots)
        assigned = False
        for name, original in config.items():
            if selector in original and original.get("core:run", False):
                staged[name][selector] = assignment
                assigned = True
                break
        if not assigned:
            for step in workflow.steps:
                if not step["run"]:
                    continue
                plugin = load_plugin(step["plugin"])
                locations = (*plugin.temporary_directory_context_paths, *plugin.output_directory_context_paths)
                if selector.replace(":", ".", 1) in locations:
                    staged[step["name"]][selector] = assignment
                    parameter = plugin.directory_parameters.get(selector.replace(":", ".", 1))
                    if parameter:
                        staged[step["name"]]["param:" + parameter] = context_path_value(
                            pairs, workflow, roots, literal=True, plugin=plugin)
                    assigned = True
                    break
        if not assigned:
            # A source plugin may publish an ordinary root field without a YAML
            # assignment. Add the user's explicit mapping at that source step.
            for index in sorted(workflow.preflight_steps):
                step = workflow.steps[index]
                if load_plugin(step["plugin"]).scene_records_return and any(
                    workflow._known(record["context"], selector.replace(":", ".", 1))
                    for record in workflow.initial_records
                ):
                    staged[step["name"]][selector] = assignment
                    assigned = True
                    break
        if not assigned:
            raise ValueError(f"Mapped root {selector} needs an explicit var:/const: definition in the recipe")
    for name, run in workflow.discovery_runs.items():
        plugin = load_plugin(staged[name]["plugin"])
        overrides = plugin.stage_settings(
            settings=deepcopy(staged[name]), params=deepcopy(run["params"]), returned=deepcopy(run["returned"]),
            path_mappings=dict(roots), file_paths=set(files), discovery_paths=set(workflow.discovery_sources),
            config_dir=config_dir,
        )
        if not isinstance(overrides, dict) or any(
            not isinstance(key, str) or not key.startswith("param:") or not key[6:].isidentifier()
            for key in overrides
        ):
            raise ValueError(f"{name}: stage_settings must return a mapping of param:name overrides")
        staged[name].update(overrides)
    rewrite_file_assignments(staged, workflow, file_bindings, roots)
    if not modern:
        for step in workflow.steps:
            plugin = load_plugin(step["plugin"])
            if plugin.scene_records_return:
                for selector in plugin.temporary_directory_context_paths:
                    if selector.startswith("const."):
                        staged[step["name"]][selector.replace(".", ":", 1)] = "path:" + remote_temp_dir

    # Rewrite explicit context filenames (mapping keys), not their field selectors.
    for name, original in config.items():
        for control in CONTEXT_CONTROLS:
            key = "core:" + control
            if key in original:
                rewritten = {}
                for template, selection in original[key].items():
                    raw = template[5:] if template.startswith("path:") else template
                    if not raw.startswith(("expr:", "var:", "const:", "collect:")):
                        filename = workflow.context_files.filename(template, workflow.initial_context)
                        try:
                            target = "path:" + mapped(filename, key)
                        except ValueError:
                            target = template  # Unused/disabled context paths need no transfer coverage.
                        rewritten[target] = selection
                    else:
                        rewritten[rewrite_scalar(template)] = selection
                staged[name][key] = rewritten

    # Relocate only explicitly loaded snapshots. Their contents remain simple JSON;
    # there is no generated whole-workflow manifest or per-step frozen value table.
    for index, filename in enumerate(sorted(context_inputs)):
        if filename not in uploads:
            continue
        import json
        snapshot = validate_snapshot(json.loads(Path(filename).read_text()), filename)
        snapshot_roots = {**roots, **{remote: remote for remote in roots.values()}}
        relocated = remap_paths(snapshot, snapshot_roots)
        if "scenes" in relocated:
            relocated["scenes"] = {remap_paths(scene_id, snapshot_roots): value for scene_id, value in relocated["scenes"].items()}
        if relocated != snapshot:
            if context_staging_dir is None:
                context_staging_dir = tempfile.mkdtemp(prefix="vhr-staged-context-")
            copy = str(Path(context_staging_dir) / f"{index}-{Path(filename).name}")
            write_json(copy, relocated)
            uploads[copy] = uploads.pop(filename)
            sources[copy] = sources.get(filename, ("context", "core:load_context"))
            root_selectors[copy] = next((root_selectors[root] for root in sorted(root_selectors, key=len, reverse=True)
                                         if _within(filename, root)), "core:load_context")
    if upload_groups is not None:
        groups = {}
        for local, remote in sorted(uploads.items()):
            step, variable = sources.get(local, ("inputs", "paths"))
            variable = next((root_selectors[root] for root in sorted(root_selectors, key=len, reverse=True)
                             if _within(local, root)), variable)
            groups.setdefault((step, variable), []).append(local)
        upload_groups.extend({"step": step, "variable": variable, "files": files}
                             for (step, variable), files in groups.items())
    return staged, dict(sorted(uploads.items())), dict(sorted(downloads.items()))
