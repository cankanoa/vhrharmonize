"""The scene variable scope exists only after a scene-setting function returns."""

import pytest

from vhrharmonize.workflow.api import run_workflow
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.values import empty_context, evaluate_settings, resolve
from workflow_helpers import install_function, stage


@pytest.mark.parametrize(
    "value",
    [
        "var:name",
        "var:$",
        "var:",
        "expr:var",
        "expr:var.name ?? 'fallback'",
        "expr:$$.var.name",
        "expr:$exists(var)",
        "expr:$lookup($, 'var')",
        "expr:$lookup($, 'v' & 'ar') ?? {}",
        {"nested": ["var:name"]},
        "collect:$",
    ],
)
def test_resolver_rejects_scene_variables_before_initialization(value):
    with pytest.raises(ValueError, match="var is unavailable until a plugin initializes scenes"):
        resolve(value, empty_context(), records=[])


def test_initial_context_has_constants_and_no_implicit_var_object():
    context = empty_context()
    assert context == {"const": {}}
    params, _, context, _ = evaluate_settings(
        {
            "const:scale": 3,
            "const:nested": {"var": "ordinary JSON key"},
            "param:literal": "literal:var:name",
            "param:scale": "expr:const.scale * 2",
            "param:context": "expr:$",
            "param:nested": "expr:const.nested.var",
        },
        context,
    )
    assert params == {
        "literal": "var:name",
        "scale": 6,
        "context": context,
        "nested": "ordinary JSON key",
    }
    assert type(params["context"]) is dict
    with pytest.raises(ValueError, match="var is unavailable"):
        evaluate_settings({"var:name": "new"}, context)


@pytest.mark.parametrize(
    "settings",
    [
        {"var:name": "new"},
        {"const:name": "var:$"},
        {"const:name": "expr:var.name ?? 'fallback'"},
    ],
)
def test_shared_cannot_initialize_or_read_scene_variables(settings, tmp_path):
    with pytest.raises(ValueError, match="var is unavailable"):
        Workflow(
            {"shared": {"plugin": "shared", "core:run": True, **settings}}, config_dir=tmp_path
        )


@pytest.mark.parametrize("scope", ["scene", "aggregate"])
@pytest.mark.parametrize(
    "settings",
    [
        {"var:name": "new"},
        {"param:value": "var:$"},
        {"param:value": {"nested": "expr:var.missing ?? 3"}},
        {"param:value": "collect:$"},
        {"core:requires": "var:path"},
    ],
)
def test_enabled_steps_cannot_use_var_before_scene_setup(scope, settings, monkeypatch, tmp_path):
    install_function(monkeypatch, "example", lambda value=None: value, scope=scope)
    with pytest.raises(ValueError, match="var is unavailable"):
        Workflow(
            {"example": {"plugin": "example", "core:run": True, "core:require_outputs": True, **settings}}, config_dir=tmp_path
        )
    # Disabled steps do not resolve or assign variables.
    assert (
        Workflow(
            {"example": {"plugin": "example", "core:run": False, "core:require_outputs": True, **settings}}, config_dir=tmp_path
        ).run()
        == []
    )


def test_ordinary_steps_run_before_scenes_and_publish_constants(monkeypatch, tmp_path):
    seen = []
    install_function(monkeypatch, "setup", lambda scale: seen.append(scale) or {"scale": scale * 2})
    install_function(
        monkeypatch,
        "source",
        lambda scale: [{"n": scale}, {"n": scale + 1}],
        scene_records_return="$",
    )
    install_function(monkeypatch, "consume", lambda n: seen.append(n))
    workflow = Workflow(
        {
            "shared": {
                "plugin": "shared",
                "core:run": True,
                "core:log_to_console": False,
                "const:scale": 3,
            },
            "setup": {
                "plugin": "setup",
                "core:run": True, "core:require_outputs": True,
                "param:scale": "const:scale",
                "const:scale": "returned:scale",
            },
            "source": {
                "plugin": "source",
                "core:run": True, "core:require_outputs": True,
                "param:scale": "const:scale",
                "var:label": "expr:var.n + 10",
            },
            "consume": {"plugin": "consume", "core:run": True, "core:require_outputs": True, "param:n": "var:label"},
        },
        config_dir=tmp_path,
    )
    workflow.run()
    assert seen == [3, 16, 17]
    assert workflow.context["const"]["scale"] == 6


@pytest.mark.parametrize("during_planning", [False, True])
def test_source_cannot_read_var_until_its_return_is_mapped(during_planning, monkeypatch, tmp_path):
    install_function(
        monkeypatch,
        "source",
        lambda value: [{"n": value}],
        scene_records_return="$",
    )
    with pytest.raises(ValueError, match="var is unavailable"):
        Workflow(
            {"source": {"plugin": "source", "core:run": True, "core:require_outputs": True, "param:value": "var:$"}},
            config_dir=tmp_path,
        )


@pytest.mark.parametrize("items", [[], [{}]])
def test_empty_scene_list_and_empty_scene_dict_are_initialized(items, monkeypatch, tmp_path):
    seen = []
    install_function(monkeypatch, "source", lambda: items, scene_records_return="$")
    install_function(monkeypatch, "scene", lambda name: seen.append(name))
    install_function(
        monkeypatch, "aggregate", lambda scenes: seen.append(scenes), scope="aggregate"
    )
    Workflow(
        {
            "shared": {"plugin": "shared", "core:run": True, "core:log_to_console": False},
            "source": {"plugin": "source", "core:run": True, "core:require_outputs": True},
            "scene": {
                "plugin": "scene",
                "core:run": True, "core:require_outputs": True,
                "var:name": "ready",
                "param:name": "var:name",
            },
            "aggregate": {"plugin": "aggregate", "core:run": True, "core:require_outputs": True, "param:scenes": "collect:$"},
        },
        config_dir=tmp_path,
    ).run()
    assert seen == (["ready", [{"name": "ready"}]] if items else [[]])


def test_workflow_without_scenes_runs_normally_locally_and_on_hpc(monkeypatch, tmp_path):
    seen = []
    install_function(monkeypatch, "example", lambda value, context: seen.append((value, context)))
    recipe = {
        "shared": {
            "plugin": "shared",
            "core:run": True,
            "core:log_to_console": False,
            "const:value": 7,
        },
        "example": {
            "plugin": "example",
            "core:run": True, "core:require_outputs": True,
            "param:value": "const:value",
            "param:context": "expr:$",
        },
    }
    assert run_workflow(recipe, config_dir=tmp_path)["example"]["processing"] == 1
    recipe.setdefault("shared", {"plugin": "shared", "core:run": True}).update({"core:save_statistics_path": None, "core:load_statistics_path": None})
    staged, uploads, downloads = stage(recipe, tmp_path)
    assert not uploads and not downloads
    assert "restore_scenes" not in staged
    Workflow(staged, config_dir=tmp_path).run()
    assert [value for value, context in seen] == [7, 7]
    assert all("var" not in context for value, context in seen)


def test_cli_and_python_raise_the_same_var_initialization_error(monkeypatch, tmp_path):
    import yaml
    from vhrharmonize.cli.main import main

    install_function(monkeypatch, "example", lambda value: value)
    config = {"example": {"plugin": "example", "core:run": True, "core:require_outputs": True, "param:value": "var:$"}}
    filename = tmp_path / "recipe.yml"
    filename.write_text(yaml.safe_dump(config, sort_keys=False))
    with pytest.raises(ValueError, match="var is unavailable") as python_error:
        run_workflow(filename)
    with pytest.raises(ValueError) as cli_error:
        main(["workflow", "--config", str(filename), "--dry-run"])
    assert str(cli_error.value) == str(python_error.value)
