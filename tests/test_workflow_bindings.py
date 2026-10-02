"""Explicit typed values, JSONata, metadata ordering and local/remote execution."""

from pathlib import Path
import json
import shutil
import pytest
import yaml
from vhrharmonize.workflow.api import load_workflow
from vhrharmonize.workflow.config import load_config
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.values import resolve, expression, evaluate_settings, Pending
from workflow_helpers import context_controls
from workflow_helpers import import_settings, copy_step, install_function, stage, transfer


@pytest.fixture
def recipe(tmp_path):
    source = tmp_path / "inputs/scene.txt"
    source.parent.mkdir()
    source.write_text("original scene")
    return {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False,
                   "core:cleanup_dirs": ["const:temp_dir"]},
        "import_files": import_settings(source, tmp_path),
    }


def test_file_and_metadata_names_are_explicit(recipe):
    first = recipe["import_files"]
    Path(first["param:search_glob"]).with_suffix(".json").write_text('{"angle":25}')
    first.update(
        {
            "param:create_metadata_json": {"document": {"to_json": "literal:expr:$replace(var.file_path, '.txt', '.json')"}},
            "var:zenith": "var:document.angle",
        }
    )
    metadata = Workflow(recipe).records[0]["context"]["var"]
    assert metadata["mul"] == first["param:search_glob"] and metadata["zenith"] == 25
    assert not {"input", "input_stem", "raw", "source", "document_path"} & metadata.keys()
    first["var:input"] = first.pop("var:mul")
    assert Workflow(recipe).records[0]["context"]["var"]["input"] == first["param:search_glob"]


@pytest.mark.parametrize(
    "first_enabled,second_enabled,expected",
    [
        (True, True, "scene_one_two_final.txt"),
        (False, True, "scene_two_final.txt"),
        (True, False, "scene_one_final.txt"),
        (False, False, "scene_final.txt"),
    ],
)
def test_disabled_steps_require_explicit_links_and_do_not_change_names(
    recipe, first_enabled, second_enabled, expected
):
    recipe.update({'file_source_1': {**(copy_step("first", "mul", suffix="_one", run=first_enabled, require_outputs=True)), "plugin": 'file_source'}, 'file_source_2': {**(copy_step("second", "first" if first_enabled else "mul", suffix="_two", run=second_enabled, require_outputs=True)), "plugin": 'file_source'}, 'file_source_3': {**(copy_step(
            "final",
            "second" if second_enabled else "first" if first_enabled else "mul",
            suffix="_final",
            folder="output_dir",
         require_outputs=True)), "plugin": 'file_source'}})
    workflow = Workflow(recipe)
    workflow.run()
    final = Path(workflow.nodes[-1].params["output_path"])
    assert final.name == expected and final.read_text() == "original scene"
    assert workflow.records[0]["context"]["var"]["suffix"] == expected[5:-4]
    assert len(workflow.nodes) == 1 + first_enabled + second_enabled


@pytest.mark.parametrize("scope", ["var", "aggregate"])
def test_disabled_steps_do_not_load_resolve_export_or_change_metadata(recipe, monkeypatch, scope):
    from vhrharmonize.workflow import engine

    load = engine.load_plugin

    def guarded(name):
        assert name != "alignment"
        return load(name)

    monkeypatch.setattr(engine, "load_plugin", guarded)
    recipe["alignment"] = {"plugin": 'alignment', 
        "core:run": False, "core:require_outputs": True,
        "core:scope": scope,
        "var:disabled": "expr:var.invalid ! var.expression",
        "param:moving_image_path": "var:undefined",
        "var:suffix": "expr:var.suffix & '_disabled'",
    }
    workflow = Workflow(recipe)
    workflow.run()
    assert not workflow.nodes and workflow.records[0]["context"]["var"]["suffix"] == ""
    assert "disabled" not in workflow.records[0]["context"]["var"]


@pytest.mark.parametrize("cached", [False, True])
def test_disabled_step_does_not_export_a_file_even_if_it_exists(recipe, cached):
    if cached:
        saved = Path(recipe["import_files"]["const:temp_dir"]) / "scene_disabled.txt"
        saved.parent.mkdir()
        saved.write_text("cached")
    recipe.update({'file_source_1': {**(copy_step("disabled", "mul", suffix="_disabled", run=False, require_outputs=True)), "plugin": 'file_source'}, 'file_source_2': {**(copy_step("final", "disabled", suffix="_final", require_outputs=True)), "plugin": 'file_source'}})
    with pytest.raises(ValueError, match="Undefined variable field: disabled"):
        Workflow(recipe).run()


def test_cached_downstream_naming_stops_upstream_processing(recipe):
    recipe.update({'file_source_1': {**(copy_step("first", "mul", suffix="_one", require_outputs=True)), "plugin": 'file_source'}, 'file_source_2': {**(copy_step("last", "first", suffix="_two", folder="output_dir", require_outputs=True)), "plugin": 'file_source'}})
    root = Path(recipe["import_files"]["const:output_dir"])
    root.mkdir()
    (root / "scene_one_two.txt").write_text("saved")
    recipe["file_source_1"]["core:require_outputs"] = False
    workflow = Workflow(recipe)
    workflow.run()
    assert workflow.counts()["file_source_1"]["processing"] == 0
    assert not Path(recipe["import_files"]["const:temp_dir"]).exists()


def test_independent_naming_branch_and_explicit_filename(recipe):
    recipe["import_files"]["var:pan_suffix"] = ""
    branch = copy_step("branch", "mul", suffix="_pan", require_outputs=True)
    branch.pop("var:suffix")
    branch["var:pan_suffix"] = "expr:var.pan_suffix & '_pan'"
    branch["var:branch"] = "expr:const.temp_dir & '/override.txt'"
    recipe.update({'file_source_1': {**(copy_step("first", "mul", suffix="_one", require_outputs=True)), "plugin": 'file_source'}, 'file_source_2': {**(branch), "plugin": 'file_source'}, 'file_source_3': {**(copy_step("last", "first", suffix="_last", folder="output_dir", require_outputs=True)), "plugin": 'file_source'}})
    workflow = Workflow(recipe)
    workflow.run()
    assert Path(workflow.records[0]["context"]["var"]["branch"]).name == "override.txt"
    assert Path(workflow.records[0]["context"]["var"]["last"]).name == "scene_one_last.txt"


def test_plugin_selection_preserves_enabled_upstream_names(recipe, monkeypatch):

    def write(input_path, output_path):
        shutil.copy2(input_path, output_path)

    install_function(
        monkeypatch,
        "finish",
        write,
        input_paths={"input_path"},
        output_paths={"output_path"},
    )
    recipe["file_source"] = copy_step("first", "mul", suffix="_one", require_outputs=True)
    recipe["finish"] = copy_step("last", "first", suffix="_last", folder="output_dir", require_outputs=True)
    recipe["finish"]["plugin"] = "finish"
    with pytest.raises(ValueError, match="Unselected step"):
        load_workflow(recipe, plugin="finish").run()
    load_workflow(recipe, plugin="file_source").run()
    workflow = load_workflow(recipe, plugin="finish")
    workflow.run()
    assert Path(workflow.records[0]["context"]["var"]["last"]).name == "scene_one_last.txt"


def test_returned_objects_and_scalars_load_explicit_context(recipe, monkeypatch):
    calls = []

    def produce(output_path):
        calls.append(output_path)
        Path(output_path).write_text("result")
        return {"object": {"quality": 7}, "scalar": 3}

    install_function(monkeypatch, "produce", produce, output_paths={"output_path"})
    recipe["produce"] = {"plugin": 'produce', 
        "core:run": True, "core:require_outputs": True,
        "param:output_path": "expr:const.output_dir & '/data.txt'",
        "var:object": "returned:object",
        "var:scalar": "returned:scalar",
        "var:double": "expr:var.scalar * 2",
    }
    recipe["produce"].update(context_controls(Path(recipe["import_files"]["const:output_dir"]) / "metadata.json"))
    workflow = Workflow(recipe)
    workflow.run()
    assert workflow.records[0]["context"]["var"]["double"] == 6
    resumed = Workflow(recipe)
    assert resumed.counts()["produce"]["loaded"] == 1
    resumed.run()
    assert len(calls) == 1
    assert resumed.records[0]["context"]["var"]["object"] == {"quality": 7}
    assert resumed.records[0]["context"]["var"]["double"] == 6


def test_aggregate_reads_returned_variables_at_its_position(recipe, monkeypatch):
    Path(recipe["import_files"]["param:search_glob"]).with_name("second.txt").write_text("second")
    recipe["import_files"]["param:search_glob"] = str(
        Path(recipe["import_files"]["param:search_glob"]).parent / "*.txt"
    )
    install_function(monkeypatch, "measure", lambda basename: len(basename))
    observed = []

    def summary(values, records):
        observed.append(records)
        return sum(values)

    install_function(monkeypatch, "summary", summary, scope="aggregate")
    recipe["measure"] = {"plugin": 'measure', 
        "core:run": True, "core:require_outputs": True,
        "param:basename": "var:basename",
        "var:score": "returned:$",
    }
    recipe["summary"] = {"plugin": 'summary', 
        "core:run": True, "core:require_outputs": True,
        "param:values": "collect:score",
        "param:records": "collect:$",
        "const:total": "returned:$",
    }
    recipe["file_source"] = copy_step("final", "mul", folder="output_dir", require_outputs=True)
    workflow = Workflow(recipe)
    workflow.run()
    assert [r["context"]["const"]["total"] for r in workflow.records] == [11, 11]
    assert all(("final" not in r for r in observed[0]))


def test_compound_options_keep_nonfile_literals(recipe, monkeypatch):
    observed = []
    install_function(monkeypatch, "inspect", lambda mask: observed.append(mask))
    recipe["inspect"] = {"plugin": 'inspect', "core:run": True, "core:require_outputs": True, "param:mask": ["include", "var:mul", "image"]}
    Workflow(recipe).run()
    assert observed == [["include", recipe["import_files"]["param:search_glob"], "image"]]


@pytest.mark.parametrize("first_enabled", [False, True])
def test_hpc_preserves_naming_and_scalar_values(recipe, tmp_path, first_enabled):
    recipe["import_files"]["var:tag"] = "ordinary string"
    recipe.update({'file_source_1': {**(copy_step("first", "mul", suffix="_one", run=first_enabled, require_outputs=True)), "plugin": 'file_source'}, 'file_source_2': {**(copy_step(
            "final", "first" if first_enabled else "mul", suffix="_last", folder="output_dir"
        , require_outputs=True)), "plugin": 'file_source'}})
    staged, uploads, _ = stage(recipe, tmp_path)
    transfer(uploads)
    workflow = Workflow(staged)
    workflow.run()
    assert Path(workflow.records[0]["context"]["var"]["final"]).name == (
        "scene_one_last.txt" if first_enabled else "scene_last.txt"
    )
    assert workflow.records[0]["context"]["var"]["tag"] == "ordinary string"


@pytest.mark.parametrize(
    "heading", ["inputs", "outputs", "output", "metadata", "run", "control:run", "metadata:name"]
)
def test_old_nested_schema_is_rejected(recipe, heading):
    recipe["file_source"] = {heading: {}}
    with pytest.raises(ValueError, match="prefix|Invalid setting"):
        Workflow(recipe)


@pytest.mark.parametrize("cached", [False, True])
def test_hpc_transports_explicit_returned_context(recipe, tmp_path, monkeypatch, cached):
    calls = []

    def produce(input_path, output_path):
        calls.append(input_path)
        shutil.copy2(input_path, output_path)
        return {"factor": 5}

    install_function(
        monkeypatch,
        "produce",
        produce,
        input_paths={"input_path"},
        output_paths={"output_path"},
    )
    recipe["produce"] = {"plugin": 'produce', 
        "core:run": True, "core:require_outputs": True,
        "param:input_path": "var:mul",
        "var:product": "expr:const.output_dir & '/product.txt'",
        "param:output_path": "var:product",
        "var:factor": "returned:factor",
    }
    recipe["file_source"] = copy_step("final", "product", suffix="_final", folder="output_dir", require_outputs=True)
    recipe["file_source"]["var:seen"] = "expr:var.factor * 2"
    recipe["import_files"]["param:scene_id"] = "literal:expr:$split(var.file_path, '/')[-1]"
    recipe["produce"].update(context_controls(tmp_path / "factor.json", "var.factor"))
    if cached:
        load_workflow(recipe, plugin="produce").run()
        calls.clear()
    staged, uploads, downloads = stage(recipe, tmp_path)
    transfer(uploads)
    remote = Workflow(staged)
    remote.run()
    assert remote.records[0]["context"]["var"]["seen"] == 10
    assert len(calls) == (0 if cached else 1)
    transfer({r: l for l, r in downloads.items()})
    local = Workflow(recipe)
    local.run()
    assert local.records[0]["context"]["var"]["seen"] == 10


def test_hpc_rewrites_declared_files_in_compound_options(recipe, tmp_path, monkeypatch):
    observed = []
    install_function(monkeypatch, "inspect", lambda mask: observed.append(mask))
    recipe["inspect"] = {"plugin": 'inspect', 
        "core:run": True, "core:require_outputs": True,
        "core:requires": "var:mul",
        "param:mask": ["include", "var:mul", "image"],
    }
    staged, uploads, _ = stage(recipe, tmp_path)
    transfer(uploads)
    Workflow(staged).run()
    assert observed == [["include", uploads[recipe["import_files"]["param:search_glob"]], "image"]]


def test_cleanup_ignores_disabled_alias_declarations(recipe):
    recipe["shared"]["core:delete_temp_steps_proactively"] = True
    recipe.update({'file_source_1': {**(copy_step("first", "mul", suffix="_one", require_outputs=True)), "plugin": 'file_source'}, 'file_source_2': {**(copy_step("alias", "first", "var:first", run=False, require_outputs=True)), "plugin": 'file_source'}, 'file_source_3': {**(copy_step("final", "first", suffix="_final", folder="output_dir", require_outputs=True)), "plugin": 'file_source'}})
    recipe["file_source_1"]["core:require_outputs"] = False
    workflow = Workflow(recipe)
    workflow.run()
    assert not Path(workflow.records[0]["context"]["var"]["first"]).exists()
    assert Path(workflow.records[0]["context"]["var"]["final"]).exists()


def test_unquoted_yaml_and_jsonata_features(recipe, tmp_path):
    settings = yaml.safe_load(
        'var:suffix: expr:var.suffix & \'_aligned\'\nparam:angle: var:solar.zenith\nvar:average: expr:$average(var.samples)\nvar:selected: expr:var.items[enabled].name\nvar:object: >-\n  expr:{"angle": var.solar.zenith, "count": $count(var.samples)}\nvar:fallback: expr:var.missing ?? 12\nvar:literal: literal:expr:ordinary text\n'
    )
    params, _, meta, _ = evaluate_settings(
        settings,
        {
            "const": {},
            "var": {
                "suffix": "_ortho",
                "solar": {"zenith": 32},
                "samples": [2, 4],
                "items": [{"enabled": True, "name": "a"}, {"enabled": False, "name": "b"}],
            },
        },
    )
    meta = meta["var"]
    assert params == {"angle": 32}
    assert meta["suffix"] == "_ortho_aligned" and meta["average"] == 3
    assert meta["selected"] == "a" and meta["object"] == {"angle": 32, "count": 2}
    assert meta["fallback"] == 12 and meta["literal"] == "expr:ordinary text"


@pytest.mark.parametrize(
    "value,expected",
    [
        ("var:x", [1, 2]),
        ("expr:var.x", [1, 2]),
        ("expr:null", None),
        ("expr:true", True),
        ("literal:var:x", "var:x"),
        ("https://example.org", "https://example.org"),
    ],
)
def test_value_types_are_preserved(value, expected):
    assert resolve(value, {"const": {}, "var": {"x": [1, 2]}}) == expected


def test_missing_metadata_and_invalid_jsonata_have_python_errors():
    for value in ["var:missing", "expr:var.missing", "expr:2 +"]:
        with pytest.raises(ValueError):
            resolve(value, {"const": {}, "var": {}})
    with pytest.raises(ValueError, match="initialized scene records"):
        resolve("collect:x", {"const": {}, "var": {"x": 1}})


def test_assignments_are_ordered_and_can_update_themselves():
    params, _, metadata, _ = evaluate_settings(
        {"var:x": "expr:var.x + 2", "param:gain": "var:x", "var:y": "expr:var.x * 3"},
        {"const": {}, "var": {"x": 1}},
    )
    assert params == {"gain": 3} and metadata["var"] == {"x": 3, "y": 9}


def test_same_invocation_return_cannot_supply_its_inputs():
    for settings in [
        {"param:gain": "returned:gain"},
        {"var:x": "returned:x", "param:gain": "var:x"},
    ]:
        with pytest.raises(ValueError, match="cannot depend"):
            evaluate_settings(settings, {"const": {}, "var": {}})


def test_nested_updates_and_deferred_null_return(recipe, monkeypatch):
    install_function(monkeypatch, "produce", lambda: {"x": None})
    recipe["produce"] = {"plugin": 'produce', "core:run": True, "core:require_outputs": True, "var:stats.angle": 42, "var:stats.x": "returned:x"}
    workflow = Workflow(recipe)
    workflow.run()
    assert workflow.records[0]["context"]["var"]["stats"] == {"angle": 42, "x": None}


def test_unused_return_values_are_omitted_from_final_json(recipe, monkeypatch):

    def produce(input_path, output_path):
        pytest.fail("Cached downstream output should bypass this producer")

    install_function(
        monkeypatch,
        "produce",
        produce,
        input_paths={"input_path"},
        output_paths={"output_path"},
    )
    recipe["produce"] = {"plugin": 'produce', 
        "core:run": True, "core:require_outputs": False,
        "param:input_path": "var:mul",
        "var:intermediate": "expr:const.temp_dir & '/intermediate.txt'",
        "param:output_path": "var:intermediate",
        "var:unused": "returned:unused",
        "var:stats.known": 4,
        "var:stats.missing": "returned:missing",
    }
    recipe["file_source"] = copy_step("final", "intermediate", folder="output_dir", require_outputs=True)
    dest = Path(recipe["import_files"]["const:output_dir"]) / "scene.txt"
    dest.parent.mkdir()
    dest.write_text("saved")
    workflow = Workflow(recipe)
    workflow.run()
    metadata = workflow.records[0]["context"]["var"]
    assert "unused" not in metadata and metadata["stats"] == {"known": 4}
    json.dumps(metadata, allow_nan=False)


def test_downloaded_context_keeps_resolved_returned_paths(recipe, tmp_path):
    recipe["file_source"] = copy_step("product", "mul", folder="output_dir", require_outputs=True)
    recipe["file_source"]["var:saved"] = "returned:$"
    recipe["import_files"]["param:scene_id"] = "literal:expr:$split(var.file_path, '/')[-1]"
    recipe["file_source"].update(context_controls(tmp_path / "saved-path.json", "var.saved"))
    staged, uploads, downloads = stage(recipe, tmp_path)
    transfer(uploads)
    remote = Workflow(staged)
    remote.run()
    transfer({r: l for l, r in downloads.items()})
    local = Workflow(recipe)
    local.run()
    assert local.records[0]["context"]["var"]["saved"] == str(
        tmp_path / "remote/output/scene.txt"  # A copied snapshot keeps its resolved values.
    )
    assert local.counts()["file_source"]["loaded"] == 1


def test_shared_function_values_reference_initialized_scene_variables(recipe, monkeypatch):
    calls = []

    def compute(gain, variables):
        calls.append((gain, variables["category"]))
        return gain

    install_function(monkeypatch, "compute", compute)
    recipe["import_files"].update({"var:category": "sensor independent", "var:base_gain": 2})
    recipe["shared"].update(
        {
            "param:gain": "expr:var.base_gain * 3",
            "param:variables": "var:$",
        }
    )
    recipe.update({'compute_1': {**({"core:run": True, "core:require_outputs": True, "var:first": "returned:$"}), "plugin": 'compute'}, 'compute_2': {**({"core:run": True, "core:require_outputs": True, "param:gain": 11, "var:second": "returned:$"}), "plugin": 'compute'}})
    Workflow(recipe).run()
    assert calls == [(6, "sensor independent"), (11, "sensor independent")]


def test_import_directory_expressions_follow_prior_assignments(recipe, tmp_path):
    recipe["shared"]["const:project"] = "scene"
    recipe["import_files"]["const:output_dir"] = "path:expr:'./products/' & const.project"
    recipe["import_files"]["const:temp_dir"] = "path:expr:const.output_dir & '/work'"
    metadata = Workflow(recipe, config_dir=tmp_path).records[0]["context"]["const"]
    assert metadata["output_dir"] == str(tmp_path / "products/scene")
    assert metadata["temp_dir"] == str(tmp_path / "products/scene/work")
    assert not Path(metadata["temp_dir"]).exists()


def test_path_arguments_explicitly_linked_to_variables_are_planned_and_staged(recipe, tmp_path):
    recipe["import_files"]["var:input_path"] = "var:mul"
    recipe["import_files"]["var:output_path"] = "expr:const.output_dir & '/inherited.txt'"
    recipe["file_source"] = {"plugin": 'file_source', 
        "core:run": True, "core:require_outputs": True,
        "param:input_path": "var:input_path",
        "param:output_path": "var:output_path",
    }
    workflow = Workflow(recipe)
    assert workflow.nodes[0].params["input_path"] == recipe["import_files"]["param:search_glob"]
    staged, uploads, _ = stage(recipe, tmp_path)
    transfer(uploads)
    Workflow(staged).run()
    assert (tmp_path / "remote/output/inherited.txt").read_text() == "original scene"


def test_shared_param_namespace_cannot_change_runner_controls(recipe):
    recipe["shared"]["param:run_from_existing"] = False
    recipe["file_source"] = copy_step("copied", "mul", folder="output_dir", require_outputs=True)
    root = Path(recipe["import_files"]["const:output_dir"])
    root.mkdir()
    (root / "scene.txt").write_text("cached")
    workflow = Workflow(recipe)
    workflow.run()
    assert workflow.counts()["file_source"]["loaded"] == 1
    assert (root / "scene.txt").read_text() == "cached"


def test_imported_relative_companions_in_compound_options_work_locally_and_on_hpc(
    recipe, tmp_path, monkeypatch
):
    source = Path(recipe["import_files"]["param:search_glob"])
    auxiliary = source.with_name("mask.json")
    auxiliary.write_text("{}")
    recipe["import_files"]["param:create_metadata_json"] = {"mask": {"path": "mask.json"}}
    recipe["import_files"]["var:mask"] = "returned:mask.0"
    calls = []
    install_function(monkeypatch, "inspect", lambda mask: calls.append(mask))
    recipe["inspect"] = {"plugin": 'inspect', 
        "core:run": True, "core:require_outputs": True,
        "core:requires": "var:mask",
        "param:mask": ["include", "var:mask", "image"],
    }
    Workflow(recipe).run()
    assert calls == [["include", str(auxiliary), "image"]]
    calls.clear()
    staged, uploads, _ = stage(recipe, tmp_path)
    transfer(uploads)
    Workflow(staged).run()
    assert calls == [["include", uploads[str(auxiliary)], "image"]]
