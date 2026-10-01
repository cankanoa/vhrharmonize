"""Edit changed YAML nodes without reserializing the surrounding document."""

from copy import deepcopy
import json
from pathlib import Path
import yaml
from yaml.nodes import MappingNode, SequenceNode, ScalarNode
from yaml.tokens import AliasToken, ScalarToken


def rewrite_yaml(source, original, updated):
    """Preserve comments, whitespace, key order and quoting outside changed nodes."""
    tree = yaml.compose(source)
    tokens = list(yaml.scan(source))
    changes = []
    newline = "\r\n" if "\r\n" in source else "\n"

    def scalar(value, style=None, indent=0):
        if isinstance(value, str):
            if style == "'":
                return "'" + value.replace("'", "''") + "'"
            if style == '"':
                return json.dumps(value, ensure_ascii=False)
        rendered = yaml.safe_dump(value, sort_keys=False, allow_unicode=True, default_flow_style=True, width=100000).rstrip()
        if rendered.endswith("\n..."):
            rendered = rendered[:-4]
        return rendered.replace("\n", newline + " " * indent)

    def fragment(value, indent):
        rendered = yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=100000).rstrip()
        return rendered.replace("\n", newline + " " * indent)

    def edit(node, old, new, key_node=None):
        if old == new:
            return
        # An alias node points at its anchor's definition. Override the alias
        # occurrence instead of accidentally modifying every use of that anchor.
        if key_node is not None and node.start_mark.index < key_node.end_mark.index:
            alias = next(token for token in tokens if isinstance(token, AliasToken) and token.start_mark.index >= key_node.end_mark.index)
            changes.append((alias.start_mark.index, alias.end_mark.index, scalar(new)))
            return
        if isinstance(node, MappingNode) and isinstance(old, dict) and isinstance(new, dict) and not node.flow_style:
            entries = {key.value: (key, value) for key, value in node.value}
            removed = [key for key in old if key not in new]
            added = [key for key in new if key not in old]
            # Context file maps often change only their filename keys.
            for before in list(removed):
                after = next((key for key in added if old[before] == new[key]), None)
                if after is not None and before in entries:
                    key, value = entries[before]
                    edit(key, before, after)
                    removed.remove(before)
                    added.remove(after)
            for name in old.keys() & new.keys():
                if old[name] == new[name]:
                    continue
                if name in entries:
                    key, value = entries[name]
                    edit(value, old[name], new[name], key)
                else:
                    added.append(name)  # inherited YAML merge key: add an override
            for name in removed:
                if name not in entries:
                    raise ValueError(f"Cannot delete inherited YAML key {name!r}")
                key, value = entries[name]
                start = source.rfind("\n", 0, key.start_mark.index) + 1
                end = value.end_mark.index
                if isinstance(value, ScalarNode):
                    line_end = source.find("\n", end)
                    line_end = len(source) if line_end < 0 else line_end + 1
                    tail = source[end:line_end]
                    comment = tail.find("#")
                    replacement = " " * key.start_mark.column + tail[comment:] if comment >= 0 else ""
                    changes.append((start, line_end, replacement))
                else:
                    changes.append((start, end, ""))
            if node is tree:
                # Step order is semantic: a new context loader must precede consumers.
                for name in list(added):
                    following = list(new)[list(new).index(name) + 1:]
                    successor = next((key for key in following if key in entries and key not in removed), None)
                    if successor is not None:
                        offset = entries[successor][0].start_mark.index
                        changes.append((offset, offset, fragment({name: new[name]}, 0) + newline))
                        added.remove(name)
            if added:
                indent = node.start_mark.column
                # Insert after the last original value, before following comments.
                last = node.value[-1][1] if node.value else node
                end = last.end_mark.index
                if isinstance(last, ScalarNode):
                    next_line = source.find("\n", end)
                    end = len(source) if next_line < 0 else next_line + 1
                prefix = "" if end == 0 or source[end - 1] == "\n" else newline
                content = prefix + " " * indent + fragment({name: new[name] for name in dict.fromkeys(added)}, indent) + newline
                changes.append((end, end, content))
            return
        if isinstance(node, SequenceNode) and isinstance(old, list) and isinstance(new, list) and len(old) == len(new):
            for child, before, after in zip(node.value, old, new):
                edit(child, before, after)
            return
        if isinstance(node, ScalarNode) and not isinstance(new, (dict, list)):
            token = next((t for t in tokens if isinstance(t, ScalarToken) and node.start_mark.index <= t.start_mark.index and t.end_mark.index <= node.end_mark.index), None)
            start, end = (token.start_mark.index, token.end_mark.index) if token else (node.start_mark.index, node.end_mark.index)
            replacement = scalar(new, node.style, node.start_mark.column)
            if node.style in {"|", ">"}:
                header = source[start:end].splitlines()[0]
                if "#" in header:
                    replacement += " " + header[header.index("#"):]
                replacement += newline
            changes.append((start, end, replacement))
        else:
            changes.append((node.start_mark.index, node.end_mark.index, scalar(new)))

    edit(tree, original, updated)
    last_start = len(source) + 1
    for start, end, value in sorted(changes, key=lambda item: (item[0], item[1]), reverse=True):
        if end > last_start:
            raise ValueError("Overlapping YAML edits")
        source = source[:start] + value + source[end:]
        last_start = start
    if yaml.safe_load(source) != updated:
        raise ValueError("Selective YAML edits did not reproduce the staged configuration")
    return source


def write_yaml_copy(source_path, destination, original, updated):
    with open(source_path, encoding="utf-8", newline="") as handle:
        source = handle.read()
    result = rewrite_yaml(source, original, updated)
    Path(destination).parent.mkdir(parents=True, exist_ok=True)
    with open(destination, "w", encoding="utf-8", newline="") as handle:
        handle.write(result)


def preparation_config(config, target):
    """Select one named cutoff without changing parameter expressions."""
    if target not in config or config[target].get("plugin") == "shared":
        raise ValueError(f"Unknown preparation step: {target}")
    result = deepcopy(config)
    after = False
    for name, settings in result.items():
        if settings.get("plugin") == "shared":
            continue
        if after:
            settings["core:run"] = False
        if name == target:
            settings["core:run"] = True
            after = True
    return result


def hpc_preparation_config(config, target):
    """Keep the full graph, but permit plugin calls only through the local cutoff."""
    if target is not None and (target not in config or config[target].get("plugin") == "shared"):
        raise ValueError(f"Unknown preparation step: {target}")
    result = deepcopy(config)
    skip = target is None
    for name, settings in result.items():
        if settings.get("plugin") == "shared":
            continue
        settings["core:skip_plugin_call"] = skip or settings.get("core:skip_plugin_call", False)
        if name == target:
            settings["core:run"] = True
            settings["core:require_outputs"] = True
            skip = True
    return result
