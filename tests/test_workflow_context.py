"""Workflow constants, independent scene variables, and their persisted context."""

from copy import deepcopy
import json
from pathlib import Path

import pytest

from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.values import (
    Deferred,
    Pending,
    evaluate_settings,
    expression,
    expression_names,
    resolve,
)
from workflow_helpers import context_controls
from workflow_helpers import import_settings, install_function, stage, transfer


@pytest.fixture
def recipe(tmp_path):
    for folder in ("a", "b"):
        source = tmp_path / folder / "same.txt"
        source.parent.mkdir()
        source.write_text(folder)
    return {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False, "const:scale": 2},
        "import_files": import_settings(tmp_path / "[ab]/same.txt", tmp_path),
    }


def test_namespaces_nested_values_and_ordered_updates():
    context = {
        "const": {"calibration": {"WV03": {"BAND_C": [0.905, -8.604]}}, "name": "global"},
        "var": {"sensor": "WV03", "name": "scene", "pending": Pending("var.pending")},
    }
    original = deepcopy(context)
    assert resolve("var:name", context) == "scene"
    assert resolve("const:name", context) == "global"
    assert resolve("const:calibration.WV03.BAND_C.0", context) == 0.905
    assert expression("$lookup(const.calibration, var.sensor).BAND_C[1]", context) == -8.604
    assert expression_names("const.calibration.WV03.BAND_C[0]") == {"const.calibration"}
    assert expression_names("var.name & const.name") == {"var.name", "const.name"}
    for source in ("var.pending", "$lookup(var, 'pending')", "$$"):
        with pytest.raises(Deferred):
            expression(source, context)
    params, _, updated, _ = evaluate_settings(
        {
            "const:calibration.WV03.BAND_C": [1.1, -2],
            "var:name": "expr:var.name & '_' & const.name",
            "param:pair": "const:calibration.WV03.BAND_C",
            "param:nested": {"name": "var:name", "gain": "expr:const.calibration.WV03.BAND_C[0]"},
        },
        context,
    )
    assert params == {"pair": [1.1, -2], "nested": {"name": "scene_global", "gain": 1.1}}
    assert updated["const"]["calibration"]["WV03"]["BAND_C"] == [1.1, -2]
    assert context == original


def test_aggregate_constants_flow_to_scenes_and_other_aggregates(recipe, monkeypatch):
    observed = []
    install_function(monkeypatch, "setup", lambda scale: {"scale": scale * 3}, scope="aggregate")
    install_function(monkeypatch, "next_setup", lambda scale: scale + 1, scope="aggregate")

    def measure(mul, scale, context, variables):
        observed.append((context, variables))
        return {"label": Path(mul).read_text(), "value": scale}

    install_function(monkeypatch, "measure", measure)
    install_function(monkeypatch, "summarize", lambda values: sum(values), scope="aggregate")
    recipe.update(
        {
            "setup": {"plugin": 'setup', 
                "core:run": True, "core:require_outputs": True,
                "param:scale": "const:scale",
                "const:scale": "returned:scale",
            },
            "next_setup": {"plugin": 'next_setup', 
                "core:run": True, "core:require_outputs": True,
                "param:scale": "const:scale",
                "const:scale": "returned:$",
            },
            "measure": {"plugin": 'measure', 
                "core:run": True, "core:require_outputs": True,
                "param:mul": "var:mul",
                "param:scale": "const:scale",
                "param:variables": "expr:$merge([const, var])",
                "param:context": "expr:$",
                "var:label": "returned:label",
                "var:value": "returned:value",
            },
            "summarize": {"plugin": 'summarize', 
                "core:run": True, "core:require_outputs": True,
                "param:values": "collect:value",
                "const:total": "returned:$",
            },
        }
    )
    workflow = Workflow(recipe)
    records = workflow.run()
    assert workflow.context["const"]["scale"] == 7
    assert workflow.context["const"]["total"] == 14
    assert len({r["id"] for r in records}) == 2
    assert [r["context"]["var"]["basename"] for r in records] == ["same", "same"]
    assert [r["context"]["var"]["label"] for r in records] == ["a", "b"]
    assert all(r["context"]["const"]["total"] == 14 for r in records)
    assert all(c["const"]["scale"] == v["scale"] == 7 for c, v in observed)
    assert all("label" not in c["var"] for c, _ in observed)


def test_shared_arguments_precedence_and_context_snapshots(recipe, monkeypatch):
    recipe["import_files"]["var:choice"] = "scene"
    recipe["shared"].update(
        {
            "const:choice": "constant",
            "const:options": {"gain": 1},
            "param:choice": "shared",
            "param:unsupported": 10,
        }
    )
    observed = []

    def inspect(choice, context, variables):
        observed.append((choice, context["const"]["options"]["gain"], variables["choice"]))
        context["const"]["options"]["gain"] = 99

    install_function(monkeypatch, "inspect", inspect)
    recipe.update({'inspect_1': {**({"core:run": True, "core:require_outputs": True, "param:variables": "var:$", "param:context": "expr:$"}), "plugin": 'inspect'}, 'inspect_2': {**({
            "core:run": True, "core:require_outputs": True,
            "param:choice": "explicit",
            "param:variables": "var:$",
            "param:context": "expr:$",
        }), "plugin": 'inspect'}})
    workflow = Workflow(recipe)
    workflow.run()
    assert observed == [("shared", 1, "scene"), ("explicit", 1, "scene")] * 2
    assert workflow.context["const"]["options"] == {"gain": 1}
    assert all(r["context"]["const"]["options"] == {"gain": 1} for r in workflow.records)


@pytest.mark.parametrize(
    "scope,assignment,error",
    [
        ("aggregate", "var:scale", "Invalid JSONata expression"),
    ],
)
def test_scope_write_validation_and_disabled_noop(recipe, monkeypatch, scope, assignment, error):
    install_function(monkeypatch, "example", lambda: None, scope=scope)
    recipe["example"] = {"plugin": 'example', "core:run": True, "core:require_outputs": True, assignment: "expr:invalid ! expression"}
    with pytest.raises(ValueError, match=error):
        Workflow(recipe)
    recipe["example"]["core:run"] = False
    Workflow(recipe).run()


@pytest.mark.parametrize(
    "template",
    [
        "var:basename",
        "expr:var.basename",
        "expr:$$.var.basename",
        "expr:var.basename ?? 'fallback'",
        "returned:$",
        "collect:basename",
        {"nested": ["var:basename"]},
        {"nested": ["returned:value"]},
    ],
)
def test_scene_constants_accept_scene_values_returns_and_collections(
    recipe, monkeypatch, template
):
    install_function(monkeypatch, "example", lambda: {"value": "shared"})
    recipe["example"] = {"plugin": 'example', "core:run": True, "core:require_outputs": True, "const:label": template}
    workflow = Workflow(recipe)
    workflow.run()
    expected = resolve(template, workflow.records[0]["context"],
                       records=[r["context"] for r in workflow.records], returned={"value": "shared"})
    assert workflow.context["const"]["label"] == expected
    assert all(r["context"]["const"]["label"] == expected for r in workflow.records)
    recipe["example"]["core:run"] = False
    workflow = Workflow(recipe)
    workflow.run()
    assert "label" not in workflow.context["const"]


@pytest.mark.parametrize("empty", [False, True])
def test_scene_constants_are_shared_once_in_order_and_survive_staging(
    recipe, tmp_path, monkeypatch, empty
):
    if empty:
        recipe["import_files"]["param:search_glob"] = str(tmp_path / "missing/*.txt")
    observed = []
    install_function(
        monkeypatch, "example", lambda before, after, token: observed.append((before, after, token))
    )
    totals = []
    install_function(
        monkeypatch,
        "summary",
        lambda scale, token: totals.append((scale, token)),
        scope="aggregate",
    )
    recipe.update(
        {
            "example": {"plugin": 'example', 
                "core:run": True, "core:require_outputs": True,
                "param:before": "const:scale",
                "const:scale": "expr:const.scale + 1",
                "param:after": "const:scale",
                "const:token": "expr:$random()",
                "param:token": "const:token",
                "const:options": {"value": "const:scale", "label": "constant text"},
                "const:options.label": "literal:var:literal_text",
                "const:list": ["alpha", "beta"],
                "const:filtered": "expr:const.list[$contains($, 'alpha')]",
                "const:from_root": "expr:$lookup($, 'const').scale",
            },
            "summary": {"plugin": 'summary', 
                "core:run": True, "core:require_outputs": True,
                "param:scale": "const:scale",
                "param:token": "const:token",
            },
        }
    )
    workflow = Workflow(recipe)
    workflow.run()
    context = workflow.context["const"]
    assert context["scale"] == context["from_root"] == 3
    assert context["options"] == {"value": 3, "label": "var:literal_text"}
    assert context["filtered"] == "alpha"
    assert observed == ([] if empty else [(2, 3, context["token"])] * 2)
    assert totals == [(3, context["token"])]
    observed.clear()
    totals.clear()
    staged, uploads, _ = stage(recipe, tmp_path)
    assert staged["example"]["const:options"]["label"] == "constant text"
    assert staged["example"]["const:options.label"] == "literal:var:literal_text"
    transfer(uploads)
    remote = Workflow(staged)
    remote.run()
    assert remote.context["const"]["scale"] == 3
    assert observed == ([] if empty else [(2, 3, remote.context["const"]["token"])] * 2)
    assert totals == [(3, remote.context["const"]["token"])]


@pytest.mark.parametrize("empty", [False, True])
def test_scene_constants_can_derive_from_runtime_constants(recipe, monkeypatch, empty):
    if empty:
        recipe["import_files"]["param:search_glob"] = "no-matching-files"
    install_function(monkeypatch, "setup", lambda: 5, scope="aggregate")
    observed = []
    install_function(monkeypatch, "example", lambda before, after: observed.append((before, after)))
    totals = []
    install_function(
        monkeypatch,
        "summary",
        lambda scale, doubled: totals.append((scale, doubled)),
        scope="aggregate",
    )
    recipe.update(
        {
            "setup": {"plugin": 'setup', "core:run": True, "core:require_outputs": True, "const:scale": "returned:$"},
            "example": {"plugin": 'example', 
                "core:run": True, "core:require_outputs": True,
                "param:before": "const:scale",
                "const:scale": "expr:const.scale + 1",
                "param:after": "const:scale",
                "const:doubled": "expr:const.scale * 2",
            },
            "summary": {"plugin": 'summary', 
                "core:run": True, "core:require_outputs": True,
                "param:scale": "const:scale",
                "param:doubled": "const:doubled",
            },
        }
    )
    workflow = Workflow(recipe)
    workflow.run()
    assert observed == ([] if empty else [(5, 6)] * 2)
    assert totals == [(6, 12)]
    assert workflow.context["const"]["doubled"] == 12


def test_scene_constants_and_returned_variables_resume_after_hpc_staging(
    recipe, tmp_path, monkeypatch
):
    calls = []

    def setup(output_path):
        calls.append("setup")
        Path(output_path).write_text("setup")
        return 5

    def process(output_path, scale):
        calls.append("process")
        Path(output_path).write_text(str(scale))
        return scale

    install_function(
        monkeypatch, "setup", setup, scope="aggregate", output_paths={"output_path"}
    )
    install_function(monkeypatch, "process", process, output_paths={"output_path"})
    totals = []
    install_function(
        monkeypatch,
        "summary",
        lambda doubled, values: totals.append((doubled, values)),
        scope="aggregate",
    )
    recipe.update(
        {
            "setup": {"plugin": 'setup', 
                "core:run": True, "core:require_outputs": True,
                "param:output_path": "expr:const.temp_dir & '/setup.txt'",
                "const:scale": "returned:$",
            },
            "process": {"plugin": 'process', 
                "core:run": True, "core:require_outputs": True,
                "const:scale": "expr:const.scale + 1",
                "const:doubled": "expr:const.scale * 2",
                "param:scale": "const:scale",
                "param:output_path": "expr:var.mul & '.processed'",
                "var:value": "returned:$",
            },
            "summary": {"plugin": 'summary', 
                "core:run": True, "core:require_outputs": True,
                "param:values": "collect:value",
                "param:doubled": "const:doubled",
            },
        }
    )
    recipe["import_files"]["param:scene_id"] = r"literal:expr:$replace(var.file_path, /^.*\/([^\/]+\/[^\/]+)$/, '$1')"
    recipe["setup"].update(context_controls(tmp_path / "scale.json", "const.scale"))
    recipe["process"].update(context_controls(tmp_path / "values.json", "var.value"))
    staged, uploads, downloads = stage(recipe, tmp_path)
    transfer(uploads)
    Workflow(staged).run()
    assert totals == [(12, [6, 6])]
    transfer({remote: local for local, remote in downloads.items()})
    totals.clear()
    local = Workflow(recipe)
    local.run()
    assert local.counts()["process"]["loaded"] == 2
    assert totals == [(12, [6, 6])]
    assert local.context["const"]["doubled"] == 12


def test_static_scene_constants_do_not_require_running_cached_functions(recipe, monkeypatch):
    def must_not_run(output_path):
        raise AssertionError("The cached scene function must not execute")

    install_function(monkeypatch, "example", must_not_run, output_paths={"output_path"})
    recipe["example"] = {"plugin": 'example', 
        "core:run": True, "core:require_outputs": True,
        "const:scale": "expr:const.scale + 1",
        "param:output_path": "expr:var.mul & '.cached'",
    }
    workflow = Workflow(recipe)
    for node in workflow.nodes:
        Path(node.params["output_path"]).write_text("cached")
    workflow.run()
    assert workflow.counts()["example"]["loaded"] == 2
    assert workflow.context["const"]["scale"] == 3


@pytest.mark.parametrize("value", ["var:basename", "expr:var.basename", "collect:basename"])
def test_import_constants_can_collect_newly_mapped_scenes(recipe, value):
    recipe["import_files"]["const:label"] = value
    workflow = Workflow(recipe)
    assert workflow.context["const"]["label"] == (["same", "same"] if value == "collect:basename" else "same")


def test_old_meta_assignment_is_rejected(recipe):
    recipe["import_files"]["meta:basename"] = "returned:basename"
    with pytest.raises(ValueError, match="Invalid setting"):
        Workflow(recipe)


def test_import_also_receives_supported_shared_parameters(recipe, tmp_path):
    recipe["shared"].update(
        {
            "param:search_glob": recipe["import_files"].pop("param:search_glob"),
            "param:create_metadata_json": {
                "source": {"to_json": "extra.json"},
                "other": {"path": "*.json"},
            },
            "param:custom_nodata_value": -9999,  # Unused by this plugin.
        }
    )
    for folder in ("a", "b"):
        (tmp_path / folder / "extra.json").write_text('{"gain":2}')
    recipe["import_files"]["var:companions"] = "returned:other"
    recipe["import_files"]["var:gain"] = "returned:source.gain"
    workflow = Workflow(recipe)
    assert workflow.records[0]["context"]["var"]["companions"] == [str(tmp_path / "a/extra.json")]
    assert workflow.records[0]["context"]["var"]["gain"] == 2
    recipe["import_files"]["param:create_metadata_json"] = {
        "source": {"to_json": "extra.json"}, "other": {"path": "missing.json"}
    }
    assert all(r["context"]["var"]["companions"] == [] for r in Workflow(recipe).records)


@pytest.mark.parametrize("empty", [False, True])
def test_explicit_constants_and_hpc_roundtrip(recipe, tmp_path, monkeypatch, empty):
    if empty:
        recipe["import_files"]["param:search_glob"] = str(tmp_path / "missing/*.txt")
    calls = []

    def setup(output_path, scale):
        calls.append(output_path)
        Path(output_path).write_text("setup")
        return {"scale": scale * 4, "path": output_path}

    install_function(
        monkeypatch, "setup", setup, scope="aggregate", output_paths={"output_path"}
    )
    observed = []
    install_function(
        monkeypatch, "consume", lambda settings: observed.append(settings), scope="aggregate"
    )
    recipe.update(
        {
            "setup": {"plugin": 'setup', 
                "core:run": True, "core:require_outputs": True,
                "param:output_path": "expr:const.output_dir & '/setup.txt'",
                "param:scale": "const:scale",
                "const:settings": "returned:$",
            },
            "consume": {"plugin": 'consume', "core:run": True, "core:require_outputs": True, "param:settings": "const:settings"},
        }
    )
    recipe["setup"].update(context_controls(tmp_path / "settings.json", "const.settings"))
    staged, uploads, downloads = stage(recipe, tmp_path)
    transfer(uploads)
    remote = Workflow(staged)
    remote.run()
    assert len(remote.records) == (0 if empty else 2)
    assert remote.context["const"]["settings"]["scale"] == 8
    transfer({remote: local for local, remote in downloads.items()})
    local = Workflow(recipe)
    assert local.counts()["setup"]["loaded"] == 1
    local.run()
    assert len(calls) == 1
    assert observed[-1] == {"scale": 8, "path": str(tmp_path / "remote/output/setup.txt")}  # Snapshots contain resolved values.
    assert local.context["const"]["settings"] == observed[-1]
    checkpoint = json.loads((tmp_path / "settings.json").read_text())
    assert "settings" in checkpoint["const"]
    # A subsequent upload must restore the same constants with remote paths again.
    staged, uploads, _ = stage(recipe, tmp_path)
    transfer(uploads)
    resumed = Workflow(staged)
    resumed.run()
    assert len(calls) == 1
    assert resumed.context["const"]["settings"]["path"] == str(tmp_path / "remote/output/setup.txt")


def test_hpc_rebases_constant_input_files(recipe, tmp_path, monkeypatch):
    reference = tmp_path / "reference.txt"
    reference.write_text("reference")
    recipe["shared"]["const:reference"] = str(reference)
    observed = []

    def inspect(input_path, context):
        observed.append((input_path, context["const"]["reference"]))

    install_function(monkeypatch, "inspect", inspect, input_paths={"input_path"})
    recipe["inspect"] = {"plugin": 'inspect', 
        "core:run": True, "core:require_outputs": True,
        "param:input_path": "const:reference",
        "param:context": "expr:$",
    }
    staged, uploads, _ = stage(recipe, tmp_path)
    transfer(uploads)
    Workflow(staged).run()
    assert observed == [(uploads[str(reference)], uploads[str(reference)])] * 2


def test_cleanup_is_generic_and_preserves_undeclared_cache(recipe, tmp_path, monkeypatch):
    recipe["import_files"]["param:search_glob"] = str(tmp_path / "a/same.txt")
    recipe["shared"]["core:delete_temp_steps_proactively"] = True

    def produce(input_path, output_path):
        Path(output_path).write_text(Path(input_path).read_text())
        Path(output_path).with_name("private-cache.txt").write_text("plugin-owned cache")

    install_function(
        monkeypatch,
        "produce",
        produce,
        input_paths={"input_path"},
        output_paths={"output_path"},
    )
    observed = []
    install_function(
        monkeypatch,
        "consume",
        lambda input_path: observed.append(Path(input_path).read_text()),
        input_paths={"input_path"},
    )
    recipe.update(
        {
            "produce": {"plugin": 'produce', 
                "core:run": True, "core:require_outputs": False,
                "param:input_path": "var:mul",
                "var:intermediate": "expr:const.temp_dir & '/intermediate.txt'",
                "param:output_path": "var:intermediate",
            },
            "consume": {"plugin": 'consume', "core:run": True, "core:require_outputs": True, "param:input_path": "var:intermediate"},
        }
    )
    Workflow(recipe).run()
    assert observed == ["a"]
    assert not (tmp_path / "temp/intermediate.txt").exists()
    assert (tmp_path / "temp/private-cache.txt").read_text() == "plugin-owned cache"
    assert (tmp_path / "a/same.txt").read_text() == "a"


def test_aliased_constant_still_requires_its_producer_when_another_consumer_is_cached(
    recipe, tmp_path, monkeypatch
):
    def setup(output_path):
        Path(output_path).write_text("setup")
        return -9999

    install_function(
        monkeypatch, "setup", setup, scope="aggregate", output_paths={"output_path"}
    )
    install_function(
        monkeypatch,
        "saved",
        lambda output_path, custom_nodata_value: None,
        scope="aggregate",
        output_paths={"output_path"},
    )
    observed = []
    plugin = install_function(
        monkeypatch, "inspect", lambda output_nodata: observed.append(output_nodata)
    )
    plugin.aliases = {"custom_nodata_value": "output_nodata"}
    recipe.update(
        {
            "setup": {"plugin": 'setup', 
                "core:run": True, "core:require_outputs": True,
                "param:output_path": "expr:const.temp_dir & '/setup.txt'",
                "const:custom_nodata_value": "returned:$",
            },
            "saved": {"plugin": 'saved', 
                "core:run": True, "core:require_outputs": True,
                "param:output_path": "expr:const.output_dir & '/saved.txt'",
                "param:custom_nodata_value": "const:custom_nodata_value",
            },
            "inspect": {"plugin": 'inspect', "core:run": True, "core:require_outputs": True, "param:custom_nodata_value": "const:custom_nodata_value"},
        }
    )
    (tmp_path / "output").mkdir()
    (tmp_path / "output/saved.txt").write_text("cached")
    workflow = Workflow(recipe).plan()
    assert workflow.nodes[0].status == "processing"
    assert workflow.nodes[1].status == "loaded"
    workflow.run()
    assert observed == [-9999, -9999]
