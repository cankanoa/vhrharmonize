"""Discovery plugins can use their own fields, parameters and HPC rewrites."""

from copy import deepcopy
import glob
from pathlib import Path
import shutil

import pytest

from vhrharmonize.plugins.base import FunctionPlugin
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.staging import stage_workflow
from vhrharmonize.workflow.staging_files import per_scene_value
from vhrharmonize.workflow.values import remap_paths, resolve
from test_explicit_context import shared
from workflow_helpers import install_function, transfer


def discovery_recipe(tmp_path, monkeypatch):
    local = tmp_path / "local"
    for scene in ("a", "b"):
        source = local / scene / "image.txt"
        source.parent.mkdir(parents=True)
        source.write_text(scene)
    calls, hooks = [], []

    def discover(locations, products, workspace):
        calls.append(locations)
        patterns = [locations] if isinstance(locations, str) else locations
        items = []
        for filename in sorted({name for pattern in patterns for name in glob.glob(pattern, recursive=True)}):
            source = Path(filename)
            item = {"identity": {"key": source.parent.name}, "asset": {"uri": str(source)},
                    "originals": [str(source)], "metadata": {"label": source.read_text()}}
            item["folders"] = {"products": resolve(products, {"var": item}),
                               "workspace": resolve(workspace, {"var": item})}
            items.append(item)
        return {"result": {"items": items}, "globals": {"version": 1}}

    plugin = install_function(monkeypatch, "catalog_reader", discover,
        scene_records_return="result.items", scene_id_return="identity.key", scene_path_return="asset.uri",
        constant_values_return="globals", source_file_protection_paths_return="originals",
        output_directory_context_paths=("var.folders.products",),
        temporary_directory_context_paths=("var.folders.workspace",),
        directory_parameters={"var.folders.products": "products", "var.folders.workspace": "workspace"},
        discovery_input_parameter="locations")

    def stage_settings(**arguments):
        hooks.append(deepcopy(arguments))
        # Mutating the hook's input must not change the recipe or enable/disable steps.
        arguments["settings"]["core:run"] = False
        return {"param:locations": [remap_paths(item["asset"]["uri"], arguments["path_mappings"])
                for item in arguments["returned"]["result"]["items"]]}

    plugin.stage_settings = stage_settings
    recipe = {"settings": shared(**{"const:root": str(tmp_path)}),
              "discover_inputs": {"plugin": "catalog_reader", "core:run": True,
                  "param:locations": str(local / "**/image.txt"),
                  "param:products": r"literal:expr:$replace(var.asset.uri, /[^\/]+$/, '') & 'products'",
                  "param:workspace": r"literal:expr:$replace(var.asset.uri, /[^\/]+$/, '') & 'workspace'",
                  "var:image": "returned:asset.uri"},
              "copy": {"plugin": "file_source", "core:run": True,
                  "core:require_outputs": "param:output_path", "param:input_path": "var:image",
                  "param:output_path": "expr:var.folders.products & '/done.txt'"}}
    remote = tmp_path / "remote"
    mappings = {
        "var:image": "expr:'" + str(remote / "inputs") + "/' & var.identity.key",
        "var:folders.products": "expr:'" + str(remote / "products") + "/' & var.identity.key",
        "var:folders.workspace": "expr:'" + str(remote / "work") + "/' & var.identity.key",
    }
    return recipe, plugin, mappings, calls, hooks


@pytest.mark.parametrize("primary_path", ["asset.uri", None])
@pytest.mark.parametrize("context_setup", [False, True])
def test_custom_discovery_fields_stage_directories_and_execute_remotely(tmp_path, monkeypatch, primary_path, context_setup):
    recipe, plugin, mappings, calls, hooks = discovery_recipe(tmp_path, monkeypatch)
    if context_setup:
        recipe = {"setup": {"core:run": True, "const:tag": "ready"}, **recipe}
    plugin.scene_path_return = primary_path  # A stable declared ID can key per-scene mappings too.
    original = deepcopy(recipe)
    remote = tmp_path / "remote"
    staged, uploads, downloads = stage_workflow(recipe, config_dir=tmp_path,
        remote_work_dir=str(remote), path_mappings=mappings)
    assert recipe == original
    assert staged["discover_inputs"]["core:run"] is True
    assert len(hooks) == 1 and hooks[0]["returned"]["result"]["items"]
    assert staged["discover_inputs"]["param:products"].startswith("literal:expr:")
    assert staged["discover_inputs"]["param:workspace"].startswith("literal:expr:")
    assert set(uploads.values()) == {str(remote / f"inputs/{scene}/image.txt") for scene in ("a", "b")}
    assert set(downloads.values()) == {str(remote / f"products/{scene}/done.txt") for scene in ("a", "b")}
    transfer(uploads)
    shutil.rmtree(tmp_path / "local")
    workflow = Workflow(staged, config_dir=tmp_path)
    workflow.run()
    assert len(calls) == 2
    assert [Path(record["context"]["var"]["folders"]["workspace"]).parent for record in workflow.records] == [remote / "work"] * 2
    assert [(remote / f"products/{scene}/done.txt").read_text() for scene in ("a", "b")] == ["a", "b"]


@pytest.mark.parametrize("label", ["locations", None])
def test_discovery_upload_labels_use_the_plugin_declaration(tmp_path, monkeypatch, label):
    recipe, plugin, _, _, _ = discovery_recipe(tmp_path, monkeypatch)
    plugin.discovery_input_parameter = label
    recipe.pop("copy")
    groups = []
    stage_workflow(recipe, config_dir=tmp_path, remote_reference_dir="/remote/inputs",
                   remote_output_dir="/remote/products", remote_temp_dir="/remote/work", upload_groups=groups)
    assert {group["variable"] for group in groups} == {"param:locations" if label else "discovery"}
    assert {group["step"] for group in groups} == {"discover_inputs"}


@pytest.mark.parametrize("saved", [False, "all", "var.image"])
def test_custom_primary_path_satisfies_outputs_from_discovery_or_context(tmp_path, monkeypatch, saved):
    recipe, _, _, calls, hooks = discovery_recipe(tmp_path, monkeypatch)
    recipe["discover_inputs"]["core:satisfies"] = {"mask": "output_path"}
    finish = recipe.pop("copy")
    recipe["mask"] = {"plugin": "file_source", "core:run": True,
                      "param:input_path": "var:missing_upstream", "var:masked": "expr:const.root & '/' & var.identity.key & '.txt'",
                      "param:output_path": "var:masked"}
    finish["param:input_path"] = "var:masked"
    recipe["copy"] = finish
    if saved:
        snapshot = tmp_path / "context.json"
        # Keep only the alias for the primary path in the partial-context case.
        selectors = saved if saved == "all" else [saved, "var.folders", "var.identity"]
        recipe["discover_inputs"]["core:save_context"] = {str(snapshot): selectors}
        Workflow(recipe, config_dir=tmp_path)
        recipe["discover_inputs"]["core:load_context"] = {str(snapshot): selectors}
    workflow = Workflow(recipe, config_dir=tmp_path)
    workflow.run()
    assert workflow.counts()["mask"]["processing"] == 0
    assert len(calls) == 1 and hooks == []
    assert [(tmp_path / f"local/{scene}/products/done.txt").read_text() for scene in ("a", "b")] == ["a", "b"]


def test_custom_discovery_context_stages_without_rerunning_the_plugin(tmp_path, monkeypatch):
    recipe, _, mappings, calls, hooks = discovery_recipe(tmp_path, monkeypatch)
    snapshot = tmp_path / "context/saved.json"
    recipe["settings"]["const:context_dir"] = str(snapshot.parent)
    recipe["discover_inputs"].update({"core:save_context": {str(snapshot): "all"},
                                      "core:load_context": {str(snapshot): "all"}})
    Workflow(recipe, config_dir=tmp_path)
    original = snapshot.read_bytes()
    remote = tmp_path / "remote"
    mappings["const:context_dir"] = str(remote / "context")
    staged, uploads, _ = stage_workflow(recipe, config_dir=tmp_path, remote_work_dir=str(remote),
        path_mappings=mappings, context_staging_dir=tmp_path / "staged-context")
    assert len(calls) == 1 and hooks == []
    assert snapshot.read_bytes() == original
    transfer(uploads)
    shutil.rmtree(tmp_path / "local")
    Workflow(staged, config_dir=tmp_path).run()
    assert len(calls) == 1
    assert [(remote / f"products/{scene}/done.txt").read_text() for scene in ("a", "b")] == ["a", "b"]


def test_runtime_discovery_registers_satisfied_outputs_before_replanning(tmp_path, monkeypatch):
    recipe, _, _, calls, _ = discovery_recipe(tmp_path, monkeypatch)
    discovery = recipe.pop("discover_inputs")
    finish = recipe.pop("copy")
    pattern = discovery["param:locations"]
    install_function(monkeypatch, "resolve_catalog", lambda: pattern)
    recipe["catalog"] = {"plugin": "resolve_catalog", "core:run": True, "const:catalog": "returned:$"}
    discovery["param:locations"] = "const:catalog"
    discovery["core:satisfies"] = {"mask": "output_path"}
    recipe["discover_inputs"] = discovery
    recipe["mask"] = {"plugin": "file_source", "core:run": True,
                      "param:input_path": "var:missing_upstream", "var:masked": "expr:const.root & '/' & var.identity.key & '.txt'",
                      "param:output_path": "var:masked"}
    finish["param:input_path"] = "var:masked"
    recipe["copy"] = finish
    Workflow(recipe, config_dir=tmp_path).run()
    assert len(calls) == 1
    assert [(tmp_path / f"local/{scene}/products/done.txt").read_text() for scene in ("a", "b")] == ["a", "b"]


def test_satisfaction_requires_declared_discovery_capabilities(tmp_path, monkeypatch):
    recipe, plugin, _, _, _ = discovery_recipe(tmp_path, monkeypatch)
    plugin.scene_path_return = None
    recipe["discover_inputs"]["core:satisfies"] = {"copy": "output_path"}
    with pytest.raises(ValueError, match="scene_records_return and scene_path_return"):
        Workflow(recipe, config_dir=tmp_path)


@pytest.mark.parametrize("declaration, value", [
    ("scene_path_return", 1), ("scene_path_return", ""),
    ("directory_parameters", []), ("directory_parameters", {"var.missing": "products"}),
    ("discovery_input_parameter", "param:locations"), ("discovery_input_parameter", False),
])
def test_invalid_discovery_declarations_are_rejected(declaration, value):
    plugin = FunctionPlugin()
    setattr(plugin, declaration, value)
    with pytest.raises(ValueError, match=declaration):
        plugin.file_features()


@pytest.mark.parametrize("overrides", [None, {"core:run": False}])
def test_invalid_staging_hook_result_is_rejected(tmp_path, monkeypatch, overrides):
    recipe, plugin, mappings, _, _ = discovery_recipe(tmp_path, monkeypatch)
    plugin.stage_settings = lambda **kwargs: overrides
    with pytest.raises(ValueError, match="stage_settings must return"):
        stage_workflow(recipe, config_dir=tmp_path, remote_work_dir=str(tmp_path / "remote"), path_mappings=mappings)


def test_scene_lookup_supports_declared_numeric_identifiers():
    expression = per_scene_value([(1, ["/one"]), (2, [])], selector="var.identity.key", literal=True)[8:]
    assert resolve(expression, {"var": {"identity": {"key": 1}}}) == ["/one"]
    assert resolve(expression, {"var": {"identity": {"key": 2}}}) == []
