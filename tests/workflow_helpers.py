"""Small fixtures using the public prefixed recipe format."""

from pathlib import Path
from copy import deepcopy
import shutil
from vhrharmonize.plugins.base import FunctionPlugin, INPUT_PATH_FEATURES, OUTPUT_PATH_FEATURES
from vhrharmonize.workflow.staging import stage_workflow


def import_settings(source, tmp_path):
    return {
        "plugin": "import_files",
        "core:run": True,
        "param:search_glob": str(source),
        "var:mul": "returned:file_path",
        "var:basename": r"expr:$replace($split(var.file_path, '/')[-1], /\.[^.]*$/, '')",
        "var:filename": "expr:$split(var.file_path, '/')[-1]",
        "var:suffix": "",
        "const:output_dir": str(tmp_path / "output"),
        "const:temp_dir": str(tmp_path / "temp"),
    }


def copy_step(name, source, destination=None, *, suffix="", run=True, folder="temp_dir", require_outputs=False):
    return {
        "plugin": "file_source",
        "core:run": run,
        "core:require_outputs": "param:output_path" if require_outputs else False,
        "param:input_path": "var:" + source,
        "var:suffix": "expr:var.suffix & " + repr(suffix),
        "var:" + name: destination
        or f"expr:const.{folder} & '/' & var.basename & var.suffix & '.txt'",
        "param:output_path": "var:" + name,
    }


def install_function(
    monkeypatch, name, function, *, scope="var", input_paths=(), output_paths=(), **file_features
):
    from vhrharmonize.workflow import registry

    plugin = FunctionPlugin()
    plugin.function = lambda: function
    plugin.scope = scope
    for feature in INPUT_PATH_FEATURES:
        setattr(plugin, feature, frozenset(input_paths))
    for feature in OUTPUT_PATH_FEATURES:
        setattr(plugin, feature, frozenset(output_paths))
    for feature, value in file_features.items():
        setattr(plugin, feature, value)

    class Entry:
        def __init__(self):
            self.name = name

        def load(self):
            return lambda: plugin

    previous = registry._entries
    monkeypatch.setattr(registry, "_entries", lambda: [*previous(), Entry()])
    return plugin


def stage(recipe, tmp_path):
    recipe = deepcopy(recipe)
    shared = next((settings for settings in recipe.values()
                   if settings.get("plugin") == "shared" and settings.get("core:run")), None)
    if shared is None:
        shared = recipe["staging_roots"] = {"plugin": "shared", "core:run": True}
    shared["const:staging_root"] = str(tmp_path)
    mappings = {"const:staging_root": str(tmp_path / "remote/reference")}
    for reference, destination in (("const:output_dir", "output"), ("var:output_dir", "output"),
                                    ("var:relative_output_dir", "output"), ("const:temp_dir", "temp")):
        if any(reference in settings for settings in recipe.values()):
            mappings[reference] = str(tmp_path / "remote" / destination)
    return stage_workflow(
        recipe,
        config_dir=str(tmp_path),
        remote_output_dir=str(tmp_path / "remote/output"),
        remote_temp_dir=str(tmp_path / "remote/temp"),
        remote_reference_dir=str(tmp_path / "remote/reference"),
        path_mappings=mappings,
    )


def transfer(mapping):
    for source, target in mapping.items():
        Path(target).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


def context_controls(filename, selectors="defined"):
    """Opt a test recipe into explicit persistence of selected context fields."""
    return {"core:" + operation: {str(filename): selectors}
            for operation in ("load_context", "save_context")}
