from copy import deepcopy
import json
from pathlib import Path
import shutil

import pytest

from test_explicit_context import shared
from test_selective_hpc_staging import transfer
from workflow_helpers import install_function
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.staging import stage_workflow
from vhrharmonize.workflow.values import resolve


def file_recipe(tmp_path, monkeypatch):
    inputs = tmp_path / "local"
    inputs.mkdir()
    for number in (1, 2):
        folder = inputs / str(number)
        for kind in ("MUL", "PAN"):
            directory = folder / ("scene_" + kind)
            directory.mkdir(parents=True)
            image = directory / f"P{number}-{'M1' if kind == 'MUL' else 'P1'}BS-image.txt"
            image.write_text(kind)
            if kind == "MUL" or number == 1:
                image.with_suffix(".RPB").write_text("rpc")
            if kind == "MUL":
                image.with_suffix(".json").write_text(json.dumps({"number": number}))
    references = tmp_path / "local-references"
    references.mkdir()
    for name in ("dem", "reference"):
        (references / (name + ".txt")).write_text(name)

    seen = []

    def process(input_path, pan_path, companions, pan_companions, references, metadata, output_path):
        assert Path(input_path).read_text() == "MUL"
        assert Path(pan_path).read_text() == "PAN"
        assert len(companions) == 2 and all(Path(p).is_file() for p in companions)
        assert len(pan_companions) == (1 if metadata["number"] == 1 else 0)
        assert all(Path(p).is_file() for p in pan_companions)
        assert [Path(p).read_text() for p in references] == ["dem", "reference"]
        seen.append((input_path, pan_path, companions, pan_companions, references))
        Path(output_path).write_text(str(metadata["number"]))

    install_function(monkeypatch, "mapped_inputs", process,
        input_paths=("input_path", "pan_path", "companions", "pan_companions", "references"),
        output_paths=("output_path",))
    config = {
        "settings": shared(**{"const:output_dir": str(tmp_path / "output"),
                              "const:dem_path": str(references / "dem.txt"),
                              "const:reference_path": str(references / "reference.txt")}),
        "import": {"plugin": "import_files", "core:run": True,
            "param:search_glob": str(inputs / "**/scene_MUL/*.txt"),
            "param:scene_id": "literal:expr:$split(var.file_path, '/')[-1]",
            "param:create_metadata_json": {
                "mul_companions": {"path": "literal:expr:$replace(var.file_path, '.txt', '.{RPB,json}')"},
                "pan_companions": {"path": "literal:expr:$replace($replace($replace(var.file_path, '_MUL/', '_PAN/'), '-M1BS-', '-P1BS-'), '.txt', '.RPB')"},
                "metadata": {"to_json": "literal:expr:$replace(var.file_path, '.txt', '.json')"}},
            "var:mul": "returned:file_path", "var:current_image_paths": "returned:file_path",
            "var:pan": "expr:$replace($replace(var.mul, '_MUL/', '_PAN/'), '-M1BS-', '-P1BS-')"},
        "process": {"plugin": "mapped_inputs", "core:run": True, "core:require_outputs": "param:output_path",
            "param:input_path": "var:current_image_paths", "param:pan_path": "var:pan",
            "param:companions": "var:mul_companions", "param:pan_companions": "var:pan_companions",
            "param:references": ["const:dem_path", "const:reference_path"], "param:metadata": "var:metadata",
            "var:current_image_paths": "expr:const.output_dir & '/' & var.scene_id",
            "param:output_path": "var:current_image_paths"},
    }
    remote = tmp_path / "remote"
    mappings = {"var:current_image_paths": str(remote / "images"), "var:pan": str(remote / "images"),
                "var:mul_companions": str(remote / "images"), "var:pan_companions": str(remote / "images"),
                "const:dem_path": str(remote / "references"), "const:reference_path": str(remote / "references"),
                "const:output_dir": str(remote / "products")}
    return config, mappings, remote, seen


@pytest.mark.parametrize("separate_companions", [False, True])
def test_flattened_scene_files_companions_and_shared_reference_folder_run_remotely(tmp_path, monkeypatch, separate_companions):
    config, mappings, remote, seen = file_recipe(tmp_path, monkeypatch)
    if separate_companions:
        mappings["var:pan"] = str(remote / "pan")
        mappings["var:mul_companions"] = str(remote / "mul-sidecars")
        mappings["var:pan_companions"] = str(remote / "pan-sidecars")
    original = deepcopy(config)
    groups = []
    staged, uploads, downloads = stage_workflow(config, config_dir=tmp_path,
        remote_work_dir=str(remote), path_mappings=mappings, upload_groups=groups)
    assert config == original
    assert staged["process"] == config["process"]  # Later mutation keeps its expression.
    assert staged["import"]["var:current_image_paths"] == "returned:file_path"
    assert set(staged) == set(config)
    assert len(staged["import"]["param:search_glob"]) == 2
    assert len(uploads) == 11  # 4 images, 5 sidecars, 2 independent references.
    assert len(downloads) == 2
    assert all(Path(target).name == Path(source).name for source, target in uploads.items())
    labels = {filename: group["variable"] for group in groups for filename in group["files"]}
    for source in uploads:
        name = Path(source).name
        expected = ("const:dem_path" if name == "dem.txt" else "const:reference_path" if name == "reference.txt"
            else "var:pan_companions" if "P1BS" in name and name.endswith(".RPB")
            else "var:pan" if "P1BS" in name else "var:current_image_paths" if name.endswith(".txt")
            else "var:mul_companions")
        assert labels[source] == expected
    if not separate_companions:
        assert staged["import"]["param:create_metadata_json"] == config["import"]["param:create_metadata_json"]
        assert staged["import"]["var:pan"] == config["import"]["var:pan"]
    transfer(uploads)
    shutil.rmtree(tmp_path / "local")
    shutil.rmtree(tmp_path / "local-references")
    Workflow(staged, config_dir=tmp_path).run()
    assert len(seen) == 2
    assert sorted(path.read_text() for path in (remote / "products").iterdir()) == ["1", "2"]


def test_file_collision_is_rejected_but_alias_of_same_file_is_allowed(tmp_path, monkeypatch):
    config, mappings, remote, _ = file_recipe(tmp_path, monkeypatch)
    mappings["var:mul"] = mappings["var:current_image_paths"]
    _, uploads, _ = stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote), path_mappings=mappings)
    assert len(uploads) == 11
    other = tmp_path / "elsewhere/dem.txt"
    other.parent.mkdir()
    other.write_text("different")
    config["settings"]["const:reference_path"] = str(other)
    with pytest.raises(ValueError, match="HPC path collision"):
        stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote), path_mappings=mappings)


@pytest.mark.parametrize("value", [{"bad": "file"}, [["nested"]], [1], None])
def test_only_paths_and_flat_lists_are_accepted(tmp_path, monkeypatch, value):
    config, mappings, remote, _ = file_recipe(tmp_path, monkeypatch)
    config["settings"]["const:bad_paths"] = value
    mappings["const:bad_paths"] = str(remote / "other")
    with pytest.raises(ValueError, match="path or flat list of paths"):
        stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote), path_mappings=mappings)


def test_all_empty_companion_lists_are_valid(tmp_path, monkeypatch):
    config, mappings, remote, _ = file_recipe(tmp_path, monkeypatch)
    config["import"]["param:create_metadata_json"]["absent"] = {"path": "*.missing"}
    mappings["var:absent"] = str(remote / "sidecars")
    staged, _, _ = stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote), path_mappings=mappings)
    assert "var:absent" not in staged["import"]


def test_lookup_templates_preserve_lists_and_remote_home_independence():
    from vhrharmonize.workflow.staging_files import per_scene_value

    pairs = [("~/work/images/one.txt", ["~/work/sidecars/one.RPB"]), ("~/work/images/two.txt", [])]
    template = per_scene_value(pairs, selector="var.file_path", path_key=True, literal=True)[8:]
    assert resolve(template, {"var": {"file_path": "/home/clusteruser/work/images/one.txt"}}) == ["~/work/sidecars/one.RPB"]
    assert resolve(template, {"var": {"file_path": "/home/clusteruser/work/images/two.txt"}}) == []
    pairs = [("~/a/same.txt", ["~/a/one.RPB"]), ("~/b/same.txt", ["~/b/two.RPB"])]
    template = per_scene_value(pairs, selector="var.file_path", path_key=True, literal=True)[8:]
    assert resolve(template, {"var": {"file_path": "/home/clusteruser/b/same.txt"}}) == ["~/b/two.RPB"]


def test_constant_file_list_and_explicit_file_override_of_directory_root(tmp_path, monkeypatch):
    config, mappings, remote, seen = file_recipe(tmp_path, monkeypatch)
    config["settings"]["const:references"] = [config["settings"].pop("const:dem_path"),
                                               config["settings"].pop("const:reference_path")]
    config["process"]["param:references"] = "const:references"
    mappings.pop("const:dem_path")
    mappings.pop("const:reference_path")
    mappings["const:references"] = str(remote / "references")
    config["settings"]["const:local_root"] = str(tmp_path / "local")
    mappings["const:local_root"] = str(remote / "original-layout")
    staged, uploads, _ = stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote), path_mappings=mappings)
    assert isinstance(staged["settings"]["const:references"], list)
    assert all(str(remote / "original-layout") not in value for value in uploads.values())
    transfer(uploads)
    shutil.rmtree(tmp_path / "local")
    shutil.rmtree(tmp_path / "local-references")
    Workflow(staged, config_dir=tmp_path).run()
    assert len(seen) == 2


def test_explicit_saved_context_rebases_file_lists_without_rerunning_import(tmp_path, monkeypatch):
    config, mappings, remote, seen = file_recipe(tmp_path, monkeypatch)
    snapshot = tmp_path / "context/import.json"
    selectors = ["defined", "var.file_path", "var.scene_id", "var.mul_companions", "var.pan_companions", "var.metadata"]
    config["import"]["core:save_context"] = {str(snapshot): selectors}
    config["import"]["core:load_context"] = {str(snapshot): selectors}
    Workflow(config, config_dir=tmp_path)
    original = snapshot.read_bytes()
    config["settings"]["const:context_root"] = str(snapshot.parent)
    mappings["const:context_root"] = str(remote / "context")
    staged, uploads, _ = stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote),
        path_mappings=mappings, context_staging_dir=tmp_path / "staged-context")
    assert snapshot.read_bytes() == original
    assert staged["import"]["param:search_glob"] == config["import"]["param:search_glob"]
    transfer(uploads)
    shutil.rmtree(tmp_path / "local")
    shutil.rmtree(tmp_path / "local-references")
    Workflow(staged, config_dir=tmp_path).run()
    assert len(seen) == 2


def test_exact_file_import_escapes_glob_characters_and_preserves_yaml_comments(tmp_path):
    import yaml
    from vhrharmonize.workflow.yaml_document import rewrite_yaml

    source = tmp_path / "local/scene[1].txt"
    source.parent.mkdir()
    source.write_text("data")
    config = {"settings": shared(**{"const:output_dir": str(tmp_path / "output")}),
              "import": {"plugin": "import_files", "core:run": True,
                         "param:search_glob": str(source.parent / "*.txt"), "var:image": "returned:file_path"},
              "copy": {"plugin": "file_source", "core:run": True, "core:require_outputs": "param:output_path",
                       "param:input_path": "var:image", "param:output_path": "expr:const.output_dir & '/copy.txt'"}}
    remote = tmp_path / "remote"
    staged, uploads, _ = stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote),
        path_mappings={"var:image": str(remote / "images"), "const:output_dir": str(remote / "products")})
    text = "# recipe\n" + yaml.safe_dump(config, sort_keys=False)
    text = text.replace("  var:image: returned:file_path", "  var:image: returned:file_path # keep")
    rewritten = rewrite_yaml(text, config, staged)
    assert rewritten.startswith("# recipe\n") and "# keep" in rewritten
    assert yaml.safe_load(rewritten) == staged
    transfer(uploads)
    source.unlink()
    Workflow(staged, config_dir=tmp_path).run()
    assert (remote / "products/copy.txt").read_text() == "data"


def test_shared_metadata_rules_survive_step_specific_relocation(tmp_path, monkeypatch):
    config, mappings, remote, seen = file_recipe(tmp_path, monkeypatch)
    config["settings"]["param:create_metadata_json"] = config["import"].pop("param:create_metadata_json")
    mappings["var:mul_companions"] = str(remote / "sidecars")
    staged, uploads, _ = stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote), path_mappings=mappings)
    assert set(staged["import"]["param:create_metadata_json"]) == {"mul_companions", "pan_companions", "metadata"}
    transfer(uploads)
    shutil.rmtree(tmp_path / "local")
    Workflow(staged, config_dir=tmp_path).run()
    assert len(seen) == 2


def test_flattening_broad_companion_globs_keeps_each_scenes_own_files(tmp_path, monkeypatch):
    config, mappings, remote, seen = file_recipe(tmp_path, monkeypatch)
    config["import"]["param:create_metadata_json"]["mul_companions"] = {"path": "*.{RPB,json}"}
    staged, uploads, _ = stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote), path_mappings=mappings)
    assert staged["import"]["param:create_metadata_json"]["mul_companions"] != {"path": "*.{RPB,json}"}
    transfer(uploads)
    shutil.rmtree(tmp_path / "local")
    Workflow(staged, config_dir=tmp_path).run()
    assert len(seen) == 2
    for image, _, companions, _, _ in seen:
        assert all(Path(p).stem == Path(image).stem for p in companions)


def test_first_file_assignment_after_import_maps_all_scenes_and_missing_products(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    for name in ("one.txt", "two.txt"):
        (source / name).write_text(name)
    config = {"settings": shared(**{"const:products": str(tmp_path / "products")}),
              "import": {"plugin": "import_files", "core:run": True,
                         "param:search_glob": str(source / "*.txt"), "var:raw": "returned:file_path"},
              "copy": {"plugin": "file_source", "core:run": True, "core:require_outputs": "param:output_path",
                       "param:input_path": "var:raw",
                       "var:result": "expr:const.products & '/' & $split(var.file_path, '/')[-1]",
                       "param:output_path": "var:result"}}
    remote = tmp_path / "remote"
    staged, uploads, downloads = stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote),
        path_mappings={"var:raw": str(remote / "images"), "var:result": str(remote / "products")})
    assert len(uploads) == len(downloads) == 2
    assert all(Path(destination).parent == remote / "products" for destination in downloads.values())
    transfer(uploads)
    shutil.rmtree(source)
    Workflow(staged, config_dir=tmp_path).run()
    assert sorted(p.read_text() for p in (remote / "products").iterdir()) == ["one.txt", "two.txt"]


def test_destination_expressions_group_all_scene_inputs_and_keep_references_shared(tmp_path, monkeypatch):
    config, mappings, remote, seen = file_recipe(tmp_path, monkeypatch)
    template = "expr:" + json.dumps(str(remote / "inputs") + "/") + " & var.scene_id & '/'"
    for selector in mappings:
        if selector.startswith("var:"):
            mappings[selector] = template
    groups = []
    staged, uploads, downloads = stage_workflow(config, config_dir=tmp_path,
        remote_work_dir=str(remote / "run"), path_mappings=mappings, upload_groups=groups)
    assert len(uploads) == 11 and len(downloads) == 2
    for source, destination in uploads.items():
        if "local-references" in source:
            assert Path(destination).parent == remote / "references"
        else:
            scene_id = Path(source).name.split("-")[0] + "-M1BS-image.txt"
            assert Path(destination) == remote / "inputs" / scene_id / Path(source).name
        assert not destination.startswith("expr:")
    assert {group["variable"] for group in groups} == set(mappings) - {"const:output_dir"}
    assert staged["process"] == config["process"]
    transfer(uploads)
    shutil.rmtree(tmp_path / "local")
    shutil.rmtree(tmp_path / "local-references")
    Workflow(staged, config_dir=tmp_path).run()
    assert len(seen) == 2
    for image, pan, companions, pan_companions, references in seen:
        assert all(Path(p).parent == Path(image).parent for p in [pan, *companions, *pan_companions])
        assert all(Path(p).parent == remote / "references" for p in references)


@pytest.mark.parametrize("template", ["expr:[]", "expr:42", "expr:null", "expr:''", "expr:var.missing", "expr:~/inputs/var:scene_id/"])
def test_destination_expression_errors_name_the_mapping(tmp_path, monkeypatch, template):
    config, mappings, remote, _ = file_recipe(tmp_path, monkeypatch)
    mappings["var:current_image_paths"] = template
    with pytest.raises(ValueError, match="HPC path_mappings var:current_image_paths:.*destination"):
        stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote), path_mappings=mappings)


def test_remote_home_expression_is_not_expanded_on_the_preparing_computer():
    from vhrharmonize.workflow.staging_files import remote_directory

    assert remote_directory("expr:'~/inputs/' & var.scene_id & '/'", {"var": {"scene_id": "P004"}},
                            "var:images") == "~/inputs/P004"


def test_directory_destination_expressions_are_resolved_before_remote_import(tmp_path):
    source = tmp_path / "source"
    for name in ("one", "two"):
        (source / name).mkdir(parents=True)
        (source / name / "image.txt").write_text(name)
    config = {"settings": shared(),
        "import": {"plugin": "import_files", "core:run": True,
            "param:search_glob": str(source / "*/*.txt"),
            "param:scene_id": "literal:expr:$split(var.file_path, '/')[-2]",
            "var:raw": "returned:file_path"},
        "copy": {"plugin": "file_source", "core:run": True, "core:require_outputs": "param:output_path",
            "param:input_path": "var:raw", "param:output_path": "expr:var.output_dir & '/copy.txt'"}}
    remote = tmp_path / "remote"
    staged, uploads, downloads = stage_workflow(config, config_dir=tmp_path, remote_work_dir=str(remote),
        path_mappings={"var:raw": "expr:" + json.dumps(str(remote / "inputs") + "/") + " & var.scene_id",
                       "var:output_dir": "expr:" + json.dumps(str(remote / "products") + "/") + " & var.scene_id"})
    assert {str(remote / "products" / scene / "copy.txt") for scene in ("one", "two")} == set(downloads.values())
    transfer(uploads)
    shutil.rmtree(source)
    Workflow(staged, config_dir=tmp_path).run()
    assert (remote / "products/one/copy.txt").read_text() == "one"
    assert (remote / "products/two/copy.txt").read_text() == "two"


def test_prepare_resolves_run_id_before_scene_expressions_and_preserves_yaml(tmp_path, monkeypatch):
    import yaml
    from vhrharmonize.slurm import prepare_slurm_plan

    config, mappings, remote, _ = file_recipe(tmp_path, monkeypatch)
    for selector in mappings:
        if selector.startswith("var:"):
            mappings[selector] = "expr:'~/inputs/{run_id}/' & var.scene_id"
    workflow_file = tmp_path / "workflow.yml"
    source = "# original recipe\n" + yaml.safe_dump(config, sort_keys=False)
    workflow_file.write_text(source)
    job = tmp_path / "job.sbatch"
    job.write_text("#!/bin/sh\n")
    hpc = tmp_path / "hpc.yml"
    hpc.write_text("# HPC settings\n" + yaml.safe_dump(dict(workflow_config=str(workflow_file),
        slurm_start_file=str(job), run_id="R42", ssh_host="example.invalid", ssh_user="user",
        remote_work_dir="~/runs/{run_id}", remote_log_dir="~/runs/{run_id}/logs", path_mappings=mappings), sort_keys=False))
    original_hpc = hpc.read_text()
    plan = prepare_slurm_plan(str(hpc))
    assert workflow_file.read_text() == source and hpc.read_text() == original_hpc
    assert Path(plan["staged_workflow_file"]).read_text().startswith("# original recipe\n")
    assert Path(plan["staged_hpc_file"]).read_text().startswith("# HPC settings\n")
    for source, target in plan["uploaded_input_paths"].items():
        if "local-references" not in source:
            assert target.startswith("~/inputs/R42/P")
    assert plan["path_mappings"]["var:current_image_paths"] == "expr:'~/inputs/R42/' & var.scene_id"
