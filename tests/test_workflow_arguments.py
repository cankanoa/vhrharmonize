"""Only explicit and shared parameters supply function arguments."""

from datetime import date
from pathlib import Path

import pytest

from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import context_controls
from workflow_helpers import import_settings, install_function, stage, transfer


@pytest.mark.parametrize("scopes", [("const",), ("var",), ("const", "var")])
@pytest.mark.parametrize(
    "shared,params,expected",
    [
        ({}, {}, (7, None, None)),
        ({"gain": 2, "custom_nodata_value": -2}, {}, (2, -2, None)),
        ({"gain": 2, "custom_nodata_value": -2}, {"gain": 3, "custom_nodata_value": -3}, (3, -3, None)),
    ],
)
def test_context_names_cannot_fill_or_override_function_arguments(
    scopes, shared, params, expected, recipe, monkeypatch
):
    observed = []
    plugin = install_function(
        monkeypatch,
        "inspect",
        lambda gain=7, output_nodata=None, variables=None: observed.append(
            (gain, output_nodata, variables)
        ),
    )
    plugin.aliases = {"custom_nodata_value": "output_nodata"}
    for scope in scopes:
        for key, value in {
            "gain": 99,
            "custom_nodata_value": 99,
            "output_nodata": 99,
            "variables": {"gain": 99},
        }.items():
            owner = "import_files" if scope == "var" else "shared"
            recipe[owner][f"{scope}:{key}"] = value
    recipe["shared"].update({"param:" + k: v for k, v in shared.items()})
    recipe["inspect"] = {"plugin": 'inspect', "core:run": True, "core:require_outputs": True, **{"param:" + k: v for k, v in params.items()}}
    Workflow(recipe).run()
    assert observed == [expected]


def test_required_function_argument_is_not_filled_from_context(recipe, monkeypatch):
    install_function(monkeypatch, "inspect", lambda gain: gain)
    recipe["shared"]["const:gain"] = 2
    recipe["import_files"]["var:gain"] = 3
    recipe["inspect"] = {"plugin": 'inspect', "core:run": True, "core:require_outputs": True}
    with pytest.raises(TypeError, match="missing a required argument: 'gain'"):
        Workflow(recipe).run()


@pytest.fixture
def recipe(tmp_path):
    source = tmp_path / "source.txt"
    source.write_text("source")
    return {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False},
        "import_files": import_settings(source, tmp_path),
    }


def test_context_paths_are_not_planned_staged_or_passed(recipe, tmp_path):
    recipe["import_files"].update(
        {
            "var:input_path": "var:mul",
            "const:output_path": str(tmp_path / "unexpected.txt"),
        }
    )
    recipe["file_source"] = {"plugin": 'file_source', "core:run": True, "core:require_outputs": True}
    workflow = Workflow(recipe)
    assert workflow.nodes[0].params == {}
    recipe["shared"].update({"core:save_statistics_path": None, "core:load_statistics_path": None})
    _, uploads, downloads = stage(recipe, tmp_path)
    assert set(uploads) == {recipe["import_files"]["param:search_glob"]}  # Ordinary discovery runs remotely.
    assert downloads == {}
    with pytest.raises(TypeError, match="missing a required argument: 'input_path'"):
        workflow.run()
    assert not (tmp_path / "unexpected.txt").exists()


def test_shared_path_references_still_work_locally_and_on_hpc(recipe, tmp_path):
    recipe["shared"].update(
        {
            "param:input_path": "var:mul",
            "param:output_path": "expr:const.output_dir & '/copy.txt'",
        }
    )
    recipe["file_source"] = {"plugin": 'file_source', "core:run": True, "core:require_outputs": True}
    staged, uploads, _ = stage(recipe, tmp_path)
    transfer(uploads)
    Workflow(staged).run()
    Workflow(recipe).run()
    assert (tmp_path / "remote/output/copy.txt").read_text() == "source"
    assert (tmp_path / "output/copy.txt").read_text() == "source"


def test_unreferenced_matching_constant_does_not_require_its_producer(
    recipe, tmp_path, monkeypatch
):
    def producer(output_path):
        pytest.fail("This branch has already been satisfied by a cached consumer")

    install_function(
        monkeypatch, "producer", producer, scope="aggregate", output_paths={"output_path"}
    )
    install_function(
        monkeypatch, "cached", producer, scope="aggregate", output_paths={"output_path"}
    )
    observed = []
    install_function(monkeypatch, "inspect", lambda gain=7: observed.append(gain))
    recipe.update(
        {
            "producer": {"plugin": 'producer', 
                "core:run": True, "core:require_outputs": False,
                "param:output_path": "expr:const.temp_dir & '/data.txt'",
                "const:gain": "returned:$",
            },
            "cached": {"plugin": 'cached', 
                "core:run": True, "core:require_outputs": True,
                "param:output_path": "expr:const.output_dir & '/saved.txt'",
                "param:gain": "const:gain",
            },
            "inspect": {"plugin": 'inspect', "core:run": True, "core:require_outputs": True},
        }
    )
    saved = tmp_path / "output/saved.txt"
    saved.parent.mkdir()
    saved.write_text("saved")
    workflow = Workflow(recipe)
    workflow.run()
    assert workflow.nodes[0].status == "unused"
    assert observed == [7]


@pytest.mark.parametrize("defaulted", [False, True])
def test_context_argument_is_never_injected(recipe, monkeypatch, defaulted):
    observed = []
    function = (
        (lambda context=None: observed.append(context))
        if defaulted
        else (lambda context: observed.append(context))
    )
    install_function(monkeypatch, "inspect", function)
    recipe["inspect"] = {"plugin": 'inspect', "core:run": True, "core:require_outputs": True}
    if defaulted:
        Workflow(recipe).run()
        assert observed == [None]
    else:
        with pytest.raises(TypeError, match="missing a required argument: 'context'"):
            Workflow(recipe).run()
        assert observed == []


@pytest.mark.parametrize(
    "binding,scope",
    [
        ("const:$", "aggregate"),
        ("var:$", "scene"),
        ("expr:$", "aggregate"),
        ("const:", "aggregate"),
        ("var:", "scene"),
    ],
)
def test_explicit_scope_arguments_track_runtime_dependencies_locally_and_on_hpc(
    recipe, tmp_path, monkeypatch, binding, scope
):
    # This test exercises reuse/staging of retained intermediate context files.
    recipe["shared"]["core:delete_temp_steps_proactively"] = False

    def produce(output_path):
        Path(output_path).write_text("produced")
        return 4

    install_function(monkeypatch, "producer", produce, scope=scope, output_paths={"output_path"})
    observed = []
    install_function(monkeypatch, "consumer", lambda context: observed.append(context))
    namespace = "const" if scope == "aggregate" else "var"
    recipe.update(
        {
            "producer": {"plugin": 'producer', 
                "core:run": True, "core:require_outputs": True,
                "param:output_path": "expr:const.temp_dir & '/produced.txt'",
                namespace + ":gain": "returned:$",
            },
            "consumer": {"plugin": 'consumer', "core:run": True, "core:require_outputs": True, "param:context": binding},
        }
    )
    recipe["producer"].update(context_controls(tmp_path / "gain.json", namespace + ".gain"))
    workflow = Workflow(recipe)
    assert workflow.nodes[1].dependencies == {0}
    workflow.run()
    # Cached producers restore returned values for a whole-scope argument, too.
    resumed = Workflow(recipe)
    resumed.run()
    assert resumed.nodes[0].status == "loaded"
    staged, uploads, _ = stage(recipe, tmp_path)
    transfer(uploads)
    remote = Workflow(staged)
    remote.run()
    assert remote.nodes[0].status == "loaded"
    assert len(observed) == 3
    for value in observed:
        assert (value[namespace] if binding == "expr:$" else value)["gain"] == 4


def test_fetch_atmosphere_uses_only_explicit_or_shared_date(tmp_path, monkeypatch):
    from vhrharmonize.plugins import fetch_atmosphere

    calls = []
    monkeypatch.setattr(
        fetch_atmosphere,
        "fetch_power_atmosphere_for_bbox",
        lambda day_utc: calls.append(day_utc) or {},
    )
    plugin = fetch_atmosphere.FetchAtmosphere()
    arguments = {"output_path": str(tmp_path / "atmosphere.json")}
    with pytest.raises(ValueError):
        plugin.run(params=arguments, shared={})
    plugin.run(params=arguments, shared={"day_utc": "2020-01-01"})
    plugin.run(params={**arguments, "day_utc": "2021-01-01"}, shared={"day_utc": "2020-01-01"})
    assert calls == [date(2020, 1, 1), date(2021, 1, 1)]
