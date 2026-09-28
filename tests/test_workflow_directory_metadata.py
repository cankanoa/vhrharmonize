"""Plugin-owned roots, explicit final JSON destinations and first-write cleanup."""

import json
from pathlib import Path

import pytest

from vhrharmonize.plugins.import_files import import_files
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.metadata import FinalMetadataWriter
from workflow_helpers import import_settings, install_function, stage, transfer


def recipe_for(tmp_path):
    for name in ("a", "b"):
        (tmp_path / f"{name}.txt").write_text(name)
    return {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False},
        "import_files": import_settings(tmp_path / "*.txt", tmp_path),
    }


def test_metadata_appends_each_completion_and_clears_an_old_file_once(tmp_path, monkeypatch):
    destination = tmp_path / "all.json"
    destination.write_text("old invalid JSON should be discarded")
    observed = []

    def process(name):
        if name == "b":
            observed.extend(json.loads(destination.read_text()))
        return {"name": name, "done": True}

    install_function(monkeypatch, "process", process)
    recipe = recipe_for(tmp_path)
    recipe["shared"]["core:output_metadata_path"] = str(destination)
    recipe["process"] = {"plugin": 'process', "core:run": True, "param:name": "var:basename", "var:result": "returned:$"}
    Workflow(recipe).run()
    entries = json.loads(destination.read_text())
    assert len(observed) == 1  # The first scene was exported before the second finished.
    assert [entry["var"]["result"]["name"] for entry in entries] == ["a", "b"]
    assert all(entry["var"]["result"]["done"] for entry in entries)
    assert all(set(entry) == {"const", "var"} for entry in entries)
    Workflow(recipe).run()  # Each workflow invocation resets its first-write tracking.
    assert len(json.loads(destination.read_text())) == 2


def test_final_json_preserves_old_entries_when_delete_first_is_false(tmp_path, monkeypatch):
    destination = tmp_path / "all.json"
    destination.write_text('{"old": true}')
    recipe = recipe_for(tmp_path)
    recipe["shared"].update(
        {"core:output_metadata_path": str(destination), "core:delete_final_json_first": False}
    )
    Workflow(recipe).run()
    entries = json.loads(destination.read_text())
    assert entries[0] == {"old": True}
    assert [entry["var"]["basename"] for entry in entries[1:]] == ["a", "b"]


def test_scene_specific_metadata_destinations_use_arbitrary_yaml_variables(tmp_path, monkeypatch):
    install_function(monkeypatch, "process", lambda: 9)
    recipe = recipe_for(tmp_path)
    recipe["shared"]["core:output_metadata_path"] = "var:my_json"
    recipe["process"] = {"plugin": 'process', 
        "core:run": True,
        "var:my_json": "expr:const.output_dir & '/' & var.basename & '.json'",
        "var:score": "returned:$",
    }
    for name in ("a", "b"):
        target = tmp_path / "output" / f"{name}.json"
        target.parent.mkdir(exist_ok=True)
        target.write_text('[{"old": true}]')
    Workflow(recipe).run()
    for name in ("a", "b"):
        entries = json.loads((tmp_path / "output" / f"{name}.json").read_text())
        assert len(entries) == 1 and entries[0]["var"]["score"] == 9
        assert entries[0]["var"]["basename"] == name


def test_json_writer_tracks_normalized_and_symlinked_destinations(tmp_path):
    target = tmp_path / "results.json"
    target.write_text("[]")
    alias = tmp_path / "alias.json"
    alias.symlink_to(target)
    writer = FinalMetadataWriter()
    writer.append(str(target), {"scene": 1})
    writer.append(str(alias), {"scene": 2})
    assert json.loads(target.read_text()) == [{"scene": 1}, {"scene": 2}]
    assert len(writer.written_paths) == 1


def test_dry_run_does_not_clear_final_json(tmp_path, monkeypatch):
    install_function(monkeypatch, "process", lambda: pytest.fail("Not run during planning"))
    target = tmp_path / "results.json"
    target.write_text("old")
    recipe = recipe_for(tmp_path)
    recipe["shared"]["core:output_metadata_path"] = str(target)
    recipe["process"] = {"plugin": 'process', "core:run": True}
    Workflow(recipe).counts()
    assert target.read_text() == "old"


def test_roots_are_not_invented_by_core_and_cleanup_requires_registration(tmp_path, monkeypatch):
    install_function(monkeypatch, "plain", lambda: None)
    workflow = Workflow({"plain": {"plugin": 'plain', "core:run": True}}, config_dir=tmp_path)
    assert workflow.context["const"] == {}
    install_function(
        monkeypatch,
        "cleanup",
        lambda output_path: None,
        output_temporary_cleanup_paths={"output_path"},
    )
    config = {
        "shared": {"plugin": 'shared', "core:run": True, "const:temp_dir": str(tmp_path)},
        "cleanup": {"plugin": 'cleanup', "core:run": True, "param:output_path": str(tmp_path / "a.bin")},
    }
    with pytest.raises(ValueError, match="temporary_directory_context_paths"):
        Workflow(config)


def test_custom_scene_json_locations_drive_output_paths_and_cleanup(tmp_path, monkeypatch):
    install_function(
        monkeypatch,
        "source",
        lambda: [
            {
                "paths": {
                    "cache": str(tmp_path / name / "work"),
                    "products": str(tmp_path / name / "out"),
                }
            }
            for name in ("a", "b")
        ],
        scene_records_return="$",
        temporary_directory_context_paths=("var.paths.cache",),
        output_directory_context_paths=("var.paths.products",),
    )
    outputs = []
    install_function(
        monkeypatch,
        "produce",
        lambda output_path: outputs.append(output_path),
        output_path_resolution_paths={"output_path"},
        output_temporary_cleanup_paths={"output_path"},
    )
    Workflow(
        {
            "source": {"plugin": 'source', "core:run": True},
            "produce": {"plugin": 'produce', "core:run": True, "param:output_path": "image.tif"},
        }
    ).run()
    assert outputs == [str(tmp_path / name / "out/image.tif") for name in ("a", "b")]


def test_import_resolves_relative_paths_and_can_publish_scene_roots(tmp_path):
    for name in ("a", "b"):
        directory = tmp_path / name
        directory.mkdir()
        (directory / "image.tif").write_text("source")
        (directory / "document.json").write_text('{"gain": 2}')
    returned = import_files(
        str(tmp_path / "*/image.tif"),
        temp_dir_scope="var",
        output_dir_scope="var",
        temp_dir="work",
        output_dir="products",
        create_metadata_json={"metadata": {"to_json": "document.json"}},
    )
    assert returned["const"] == {}
    for name, scene in zip(("a", "b"), returned["scenes"]):
        assert scene["temp_dir"] == str(tmp_path / name / "work")
        assert scene["output_dir"] == str(tmp_path / name / "products")
        assert scene["metadata"]["gain"] == 2
        assert set((tmp_path / name).iterdir()) == {
            tmp_path / name / "image.tif", tmp_path / name / "document.json"
        }


@pytest.mark.parametrize("protect", [True, False])
def test_source_protection_has_an_explicit_shared_switch(tmp_path, monkeypatch, protect):
    raw = tmp_path / "original.txt"
    raw.write_text("original")
    install_function(
        monkeypatch,
        "source",
        lambda: [{"files": [str(raw)]}],
        scene_records_return="$",
        source_file_protection_paths_return="files",
    )
    install_function(
        monkeypatch,
        "replace",
        lambda output: Path(output).write_text("replacement"),
        output_path_resolution_paths={"output"},
    )
    config = {
        "shared": {"plugin": 'shared', "core:run": True, "core:protect_source_files": protect},
        "source": {"plugin": 'source', "core:run": True},
        "replace": {"plugin": 'replace', "core:run": True, "param:output": str(raw)},
    }
    if protect:
        with pytest.raises(ValueError, match="protected input"):
            Workflow(config).run()
        assert raw.read_text() == "original"
    else:
        Workflow(config).run()
        assert raw.read_text() == "replacement"


def test_hpc_rewrites_and_downloads_explicit_final_json_path(tmp_path, monkeypatch):
    install_function(monkeypatch, "process", lambda name: name)
    recipe = recipe_for(tmp_path)
    target = tmp_path / "combined.json"
    target.write_text('[{"previous": true}]')
    recipe["shared"].update(
        {"core:output_metadata_path": str(target), "core:delete_final_json_first": False}
    )
    recipe["process"] = {"plugin": 'process', 
        "core:run": True,
        "param:name": "var:basename",
        "var:processed": "returned:$",
    }
    staged, uploads, downloads = stage(recipe, tmp_path)
    assert str(target) in uploads and str(target) in downloads
    transfer(uploads)
    Workflow(staged).run()
    entries = json.loads(Path(downloads[str(target)]).read_text())
    assert entries[0] == {"previous": True} and len(entries) == 3


def test_final_metadata_restores_cached_results(tmp_path, monkeypatch):
    calls = []

    def produce(name, output_path):
        calls.append(name)
        Path(output_path).write_text(name)
        return {"score": 23}

    install_function(monkeypatch, "produce", produce, output_paths={"output_path"})
    recipe = recipe_for(tmp_path)
    destination = tmp_path / "final.json"
    recipe["shared"]["core:output_metadata_path"] = str(destination)
    recipe["produce"] = {"plugin": 'produce', 
        "core:run": True,
        "param:name": "var:basename",
        "param:output_path": "expr:const.output_dir & '/' & var.filename",
        "var:score": "returned:score",
    }
    Workflow(recipe).run()
    resumed = Workflow(recipe)
    resumed.run()
    assert calls == ["a", "b"]
    assert resumed.counts()["produce"]["loaded"] == 2
    assert [v["var"]["score"] for v in json.loads(destination.read_text())] == [23, 23]


def test_final_aggregate_metadata_includes_returned_constants(tmp_path, monkeypatch):
    install_function(monkeypatch, "aggregate", lambda names: {"names": names}, scope="aggregate")
    recipe = recipe_for(tmp_path)
    destination = tmp_path / "final.json"
    recipe["shared"]["core:output_metadata_path"] = str(destination)
    recipe["aggregate"] = {"plugin": 'aggregate', 
        "core:run": True,
        "param:names": "collect:basename",
        "const:summary": "returned:$",
    }
    Workflow(recipe).run()
    entries = json.loads(destination.read_text())
    assert [v["var"]["basename"] for v in entries] == ["a", "b"]
    assert all(v["const"]["summary"] == {"names": ["a", "b"]} for v in entries)


def test_discovery_only_hpc_exports_final_metadata(tmp_path):
    recipe = recipe_for(tmp_path)
    target = tmp_path / "combined.json"
    recipe["shared"]["core:output_metadata_path"] = str(target)
    staged, uploads, downloads = stage(recipe, tmp_path)
    assert str(target) not in uploads and str(target) in downloads
    Workflow(staged).run()
    assert len(json.loads(Path(downloads[str(target)]).read_text())) == 2


def test_runtime_directory_return_drives_temporary_cleanup(tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    install_function(
        monkeypatch,
        "directories",
        lambda: str(cache),
        temporary_directory_context_paths=("const.work",),
    )
    install_function(
        monkeypatch,
        "produce",
        lambda output_path: Path(output_path).write_text("made"),
        output_paths={"output_path"},
    )
    observed = []
    install_function(
        monkeypatch,
        "consume",
        lambda input_path: observed.append(Path(input_path).read_text()),
        input_paths={"input_path"},
    )
    recipe = {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False, "core:delete_temp_steps_proactively": True},
        "directories": {"plugin": 'directories', "core:run": True, "const:work": "returned:$"},
        "produce": {"plugin": 'produce', "core:run": True, "param:output_path": str(cache / "intermediate.txt")},
        "consume": {"plugin": 'consume', "core:run": True, "param:input_path": str(cache / "intermediate.txt")},
    }
    Workflow(recipe).run()
    assert observed == ["made"]
    assert not (cache / "intermediate.txt").exists()


def test_hpc_preserves_distinct_scene_directory_locations(tmp_path, monkeypatch):
    scenes = [
        {"paths": {"work": str(tmp_path / name / "temp"), "products": str(tmp_path / name / "out")}}
        for name in ("a", "b")
    ]
    install_function(
        monkeypatch,
        "source",
        lambda: scenes,
        scene_records_return="$",
        temporary_directory_context_paths=("var.paths.work",),
        output_directory_context_paths=("var.paths.products",),
    )
    install_function(
        monkeypatch,
        "produce",
        lambda output_path: Path(output_path).write_text("made"),
        output_paths={"output_path"},
    )
    recipe = {
        "shared": {"plugin": 'shared', "core:run": True, 
            "core:log_to_console": False,
            "core:output_metadata_path": "expr:var.paths.products & '/final.json'",
        },
        "source": {"plugin": 'source', "core:run": True},
        "produce": {"plugin": 'produce', "core:run": True, "param:output_path": "image.txt"},
    }
    staged, uploads, downloads = stage(recipe, tmp_path)
    assert not uploads and len(downloads) == 4
    assert len(set(downloads.values())) == 4
    Workflow(staged).run()
    for name in ("a", "b"):
        image = tmp_path / name / "out/image.txt"
        assert Path(downloads[str(image)]).read_text() == "made"
        result = tmp_path / name / "out/final.json"
        entry = json.loads(Path(downloads[str(result)]).read_text())[0]
        assert entry["var"]["paths"]["products"] == str(Path(downloads[str(image)]).parent)


def test_missing_declared_temp_directory_is_an_error(tmp_path, monkeypatch):
    install_function(
        monkeypatch,
        "source",
        lambda: [{}],
        scene_records_return="$",
        temporary_directory_context_paths=("var.work",),
    )
    install_function(
        monkeypatch,
        "produce",
        lambda output_path: None,
        output_temporary_cleanup_paths={"output_path"},
    )
    with pytest.raises(ValueError, match="temp_dir is required"):
        Workflow(
            {
                "source": {"plugin": 'source', "core:run": True},
                "produce": {"plugin": 'produce', "core:run": True, "param:output_path": str(tmp_path / "image")},
            }
        )
