"""Plugins opt into each file operation independently."""

import json
from pathlib import Path

import pytest
import rasterio

from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import import_settings, install_function, stage, transfer, context_controls


@pytest.fixture
def recipe(tmp_path):
    source = tmp_path / "source" / "image.txt"
    source.parent.mkdir()
    source.write_text("source")
    return {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False, "core:show_progress": False, "core:save_statistics_path": None, "core:load_statistics_path": None},
        "import_files": import_settings(source, tmp_path),
    }


@pytest.mark.parametrize("direction", [None, "output"])
def test_path_resolution_does_not_enable_other_file_operations(
    recipe, tmp_path, monkeypatch, direction
):
    observed = []
    features = {direction + "_path_resolution_paths": {direction + "_path"}} if direction else {}
    install_function(
        monkeypatch, "example", lambda **kwargs: observed.append(kwargs), **features
    ).passthrough = True
    recipe["example"] = {"plugin": 'example', 
        "core:run": True, "core:require_outputs": True,
        "core:calculate_overviews": True,
        "param:input_path": "missing.txt",
        "param:output_path": "absent/output.tif",
    }
    Workflow(recipe, config_dir=tmp_path).run()
    expected = {"input_path": "missing.txt", "output_path": "absent/output.tif"}
    if direction:
        expected[direction + "_path"] = str(tmp_path / expected[direction + "_path"])
    assert observed == [expected]
    assert not (tmp_path / "output/absent").exists()


def test_input_checks_and_output_parent_creation_select_only_their_arguments(
    recipe, tmp_path, monkeypatch
):
    source = recipe["import_files"]["param:search_glob"]
    destination = tmp_path / "created/result.txt"
    ignored = tmp_path / "untouched/result.txt"
    observed = []
    install_function(
        monkeypatch,
        "example",
        lambda **kwargs: observed.append(destination.parent.is_dir()),
        input_existence_check_paths={"source"},
        output_parent_creation_paths={"destination"},
    ).passthrough = True
    recipe["example"] = {"plugin": 'example', 
        "core:run": True, "core:require_outputs": True,
        "param:source": source,
        "param:optional_input": str(tmp_path / "missing.txt"),
        "param:destination": str(destination),
        "param:other_output": str(ignored),
    }
    Workflow(recipe).run()
    assert observed == [True]
    assert not ignored.parent.exists()
    recipe["example"]["param:source"] = str(tmp_path / "missing.txt")
    with pytest.raises(FileNotFoundError, match="input does not exist"):
        Workflow(recipe).run()


@pytest.mark.parametrize(
    "reuse,validate,calls", [(False, True, 2), (True, False, 0), (True, True, 1)]
)
def test_validation_and_reuse_are_independent(
    recipe, tmp_path, monkeypatch, reuse, validate, calls
):
    output = tmp_path / "output.json"
    ignored = tmp_path / "ignored.json"
    output.write_text("broken")
    ignored.write_text("broken")
    observed = []

    def produce(output_path, other_output):
        observed.append(Path(output_path).read_text())
        Path(output_path).write_text("{}")

    install_function(
        monkeypatch,
        "example",
        produce,
        output_reuse_paths={"output_path"} if reuse else set(),
        output_validation_paths={"output_path"} if validate else set(),
    )
    recipe["example"] = {"plugin": 'example', 
        "core:run": True, "core:require_outputs": True,
        "param:output_path": str(output),
        "param:other_output": str(ignored),
    }
    Workflow(recipe).run()
    Workflow(recipe).run()
    assert len(observed) == calls
    if observed:
        assert observed[0] == "broken"  # Invalid removal was not enabled.
    assert ignored.read_text() == "broken"


def test_invalid_removal_and_post_validation_are_independent(recipe, tmp_path, monkeypatch):
    output = tmp_path / "invalid.json"
    ignored = tmp_path / "ignored.json"
    output.write_text("broken")
    ignored.write_text("broken")
    observed = []
    plugin = install_function(
        monkeypatch,
        "example",
        lambda output_path, ignored_path: observed.append(
            (Path(output_path).exists(), Path(ignored_path).exists())
        ),
        output_invalid_removal_paths={"output_path"},
    )
    recipe["example"] = {"plugin": 'example', 
        "core:run": True, "core:require_outputs": True,
        "param:output_path": str(output),
        "param:ignored_path": str(ignored),
    }
    Workflow(recipe).run()
    assert observed == [(False, True)]
    plugin.output_validation_paths = {"output_path"}
    with pytest.raises(RuntimeError, match="did not produce valid declared outputs"):
        Workflow(recipe).run()


@pytest.mark.parametrize("enabled", [False, True])
def test_overviews_only_touch_selected_outputs(
    recipe, tmp_path, monkeypatch, make_test_raster, enabled
):
    first, second = tmp_path / "first.tif", tmp_path / "second.tif"

    def produce(first_path, second_path):
        make_test_raster(Path(first_path), width=16, height=16)
        make_test_raster(Path(second_path), width=16, height=16)

    install_function(
        monkeypatch, "example", produce, output_overview_calculation_paths={"first_path"}
    )
    recipe["shared"]["param:window_scales"] = [2]
    recipe["example"] = {"plugin": 'example', 
        "core:run": True, "core:require_outputs": True,
        "core:calculate_overviews": enabled,
        "param:first_path": str(first),
        "param:second_path": str(second),
    }
    Workflow(recipe).run()
    with rasterio.open(first) as dataset:
        assert dataset.overviews(1) == ([2] if enabled else [])
    with rasterio.open(second) as dataset:
        assert dataset.overviews(1) == []


def test_cleanup_selects_files_after_their_consumers_finish(recipe, tmp_path, monkeypatch):
    first, second = tmp_path / "temp/first.txt", tmp_path / "temp/second.txt"

    def produce(first_path, second_path):
        first.parent.mkdir()
        Path(first_path).write_text("first")
        Path(second_path).write_text("second")

    observed = []
    install_function(
        monkeypatch,
        "produce",
        produce,
        output_dependency_paths={"first_path", "second_path"},
        output_temporary_cleanup_paths={"first_path"},
    )
    install_function(
        monkeypatch,
        "consume",
        lambda files: observed.append([Path(p).read_text() for p in files]),
        input_dependency_paths={"files"},
    )
    install_function(
        monkeypatch,
        "inspect",
        lambda reference: observed.append(Path(reference).read_text()),
        input_protection_paths={"reference"},
    )
    recipe["shared"]["core:delete_temp_steps_proactively"] = True
    recipe["produce"] = {"plugin": 'produce', 
        "core:run": True, "core:require_outputs": False,
        "param:first_path": str(first),
        "param:second_path": str(second),
    }
    recipe["consume"] = {"plugin": 'consume', "core:run": True, "core:require_outputs": True, "param:files": [str(first), str(second)]}
    recipe["inspect"] = {"plugin": 'inspect', "core:run": True, "core:require_outputs": True, "param:reference": str(first)}
    workflow = Workflow(recipe)
    assert workflow.nodes[-1].dependencies == set()  # Protection is not a dependency declaration.
    workflow.run()
    assert observed == [["first", "second"], "first"]
    assert not first.exists()
    assert second.exists()


@pytest.mark.parametrize("protect", [False, True])
def test_input_protection_is_independent_of_existence_checks(
    recipe, tmp_path, monkeypatch, protect
):
    reference = tmp_path / "reference.txt"
    reference.write_text("reference")
    alias = tmp_path / "alias.txt"
    alias.symlink_to(reference)
    observed = []
    install_function(
        monkeypatch,
        "example",
        lambda reference, output: observed.append(reference),
        input_existence_check_paths={"reference"},
        input_protection_paths={"reference"} if protect else set(),
        output_path_resolution_paths={"output"},
    )
    recipe["example"] = {"plugin": 'example', 
        "core:run": True, "core:require_outputs": True,
        "param:reference": str(reference),
        "param:output": str(alias),
    }
    if protect:
        with pytest.raises(ValueError, match="protected input"):
            Workflow(recipe).run()
        assert observed == []
    else:
        Workflow(recipe).run()
        assert observed == [str(reference)]
    assert reference.read_text() == "reference"


@pytest.mark.parametrize("target", [False, True])
def test_persistent_target_selection_controls_unused_branches(
    recipe, tmp_path, monkeypatch, target
):
    intermediate, final = tmp_path / "intermediate.txt", tmp_path / "final.txt"
    final.write_text("cached")
    observed = []
    install_function(
        monkeypatch,
        "produce",
        lambda output_path: observed.append("produce"),
        output_dependency_paths={"output_path"},

    )
    install_function(
        monkeypatch,
        "consume",
        lambda input_path, output_path: pytest.fail("Already cached"),
        input_dependency_paths={"input_path"},
        output_reuse_paths={"output_path"},
    )
    recipe["produce"] = {"plugin": 'produce', "core:run": True, "core:require_outputs": "param:output_path" if target else False, "param:output_path": str(intermediate)}
    recipe["consume"] = {"plugin": 'consume', 
        "core:run": True, "core:require_outputs": True,
        "param:input_path": str(intermediate),
        "param:output_path": str(final),
    }
    workflow = Workflow(recipe)
    assert workflow.nodes[1].dependencies == {0}
    workflow.run()
    assert observed == (["produce"] if target else [])


@pytest.mark.parametrize("check_collision", [False, True])
def test_collision_checks_are_separate_from_resolution(
    recipe, tmp_path, monkeypatch, check_collision
):
    observed = []
    install_function(
        monkeypatch,
        "example",
        lambda output_path: observed.append(output_path),
        output_path_resolution_paths={"output_path"},
        output_collision_check_paths={"output_path"} if check_collision else set(),
    )
    recipe.update({'example_1': {**({"core:run": True, "core:require_outputs": "param:output_path", "param:output_path": str(tmp_path / "result.txt")}), "plugin": 'example'}, 'example_2': {**({"core:run": True, "core:require_outputs": "param:output_path", "param:output_path": str(tmp_path / "result.txt")}), "plugin": 'example'}})
    if check_collision:
        with pytest.raises(ValueError, match="Output path collision"):
            Workflow(recipe)
    else:
        Workflow(recipe).run()
        assert len(observed) == 2


@pytest.mark.parametrize("checkpoint", [False, True])
def test_context_persistence_is_explicit(recipe, tmp_path, monkeypatch, checkpoint):
    first, second = tmp_path / "first.txt", tmp_path / "second.txt"
    observed = []

    def produce(first_path, second_path):
        observed.append(1)
        Path(first_path).write_text("first")
        Path(second_path).write_text("second")
        return 7

    install_function(
        monkeypatch,
        "example",
        produce,
        output_reuse_paths={"first_path"},

    )
    recipe["example"] = {"plugin": 'example', 
        "core:run": True, "core:require_outputs": "param:first_path",
        "param:first_path": str(first),
        "param:second_path": str(second),
        "var:gain": "returned:$",
    }
    if checkpoint:
        recipe["example"].update(context_controls(tmp_path / "gain.json", "var.gain"))
    install_function(monkeypatch, "use_gain", lambda gain: gain)
    recipe["consumer"] = {"plugin": "use_gain", "core:run": True, "core:require_outputs": True, "param:gain": "var:gain"}
    Workflow(recipe).run()
    resumed = Workflow(recipe)
    resumed.run()
    assert observed == ([1] if checkpoint else [1, 1])
    assert resumed.records[0]["context"]["var"]["gain"] == 7
    assert not Path(str(first) + ".context.json").exists()
    if checkpoint:
        assert (
            next(iter(json.loads((tmp_path / "gain.json").read_text())["scenes"].values()))["gain"] == 7
        )


def test_hpc_staging_and_download_selections_are_independent(recipe, tmp_path, monkeypatch):
    source = recipe["import_files"]["param:search_glob"]
    shared_file = tmp_path / "shared.txt"
    shared_file.write_text("shared")
    outputs = {
        name: str(tmp_path / "output" / (name + ".txt")) for name in ("cache", "result", "private")
    }
    install_function(
        monkeypatch,
        "example",
        lambda **kwargs: None,
        input_hpc_staging_paths={"source"},
        output_path_resolution_paths=set(outputs),
        output_hpc_staging_paths={"cache"},
        output_hpc_download_paths={"result"},
    ).passthrough = True
    recipe["example"] = {"plugin": 'example', 
        "core:run": True, "core:require_outputs": True,
        "param:source": source,
        "param:shared_file": str(shared_file),
        **{"param:" + k: v for k, v in outputs.items()},
    }
    staged, uploads, downloads = stage(recipe, tmp_path)
    assert set(uploads) == {source}
    assert downloads == {outputs["result"]: str(tmp_path / "remote/output/result.txt")}
    transfer(uploads)
    node = Workflow(staged).nodes[0]
    assert node.params == {
        "source": uploads[source],
        "shared_file": str(tmp_path / "remote/reference/shared.txt"),
        "cache": str(tmp_path / "remote/output/cache.txt"),
        "result": str(tmp_path / "remote/output/result.txt"),
        "private": str(tmp_path / "remote/output/private.txt"),
    }
