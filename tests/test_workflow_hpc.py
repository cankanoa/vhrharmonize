from pathlib import Path
import shutil
import yaml
from vhrharmonize.slurm import prepare_slurm_plan
from vhrharmonize.workflow.config import load_config
from vhrharmonize.workflow.engine import Workflow


def test_hpc_prepares_ordered_workflow_and_executes_staged_recipe(tmp_path, make_test_raster):
    source = make_test_raster(tmp_path / "source" / "image.tif")
    recipe = tmp_path / "recipe.yml"
    recipe.write_text(
        yaml.safe_dump(
            {
                "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False},
                "import_files": {"plugin": 'import_files', 
                    "core:run": True,
                    "param:search_glob": str(source),
                    "var:raw": "returned:file_path",
                    "var:input": "returned:file_path",
                    "var:input_name": "expr:$split(var.file_path, '/')[-1]",
                    "var:input_stem": r"expr:$replace($split(var.file_path, '/')[-1], /\.[^.]*$/, '')",
                    "var:input_dir": r"expr:$replace(var.file_path, /[^\/]+$/, '')",
                    "const:output_dir": str(tmp_path / "out"),
                    "const:temp_dir": str(tmp_path / "temp"),
                },
                "file_source": {"plugin": 'file_source', 
                    "core:run": True,
                    "param:input_path": "var:raw",
                    "var:result": "expr:const.output_dir & '/' & var.input_name",
                    "param:output_path": "var:result",
                },
            },
            sort_keys=False,
        )
    )
    sbatch = tmp_path / "job.sbatch"
    sbatch.write_text(
        '#!/bin/bash\n#SBATCH --output=../logs/slurm-%j.out\nvhr workflow --config "$1"\n'
    )
    hpc = tmp_path / "hpc.yml"
    hpc.write_text(
        yaml.safe_dump(
            {
                "workflow_config": str(recipe),
                "slurm_start_file": str(sbatch),
                "run_id": "example",
                "ssh_host": "example.invalid",
                "ssh_user": "user",
                **{
                    f"remote_{key}_dir": str(tmp_path / "remote" / key)
                    for key in ["output", "log", "temp", "reference"]
                },
            },
            sort_keys=False,
        )
    )
    plan = prepare_slurm_plan(str(hpc))
    assert str(source) in plan["uploaded_input_paths"]
    assert plan["download_output_paths"][str(tmp_path / "out" / "image.tif")] == str(
        tmp_path / "remote" / "output" / "image.tif"
    )
    config = load_config(plan["staged_workflow_file"])
    assert list(config) == ["shared", "restore_scenes", "import_files", "file_source"]
    for local, remote in plan["uploaded_input_paths"].items():
        Path(remote).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(local, remote)
    Workflow(config).run()
    assert (tmp_path / "remote" / "output" / "image.tif").read_bytes() == source.read_bytes()


def test_hpc_stages_declared_rpc_companions_beside_the_raster(tmp_path, make_test_raster):
    from vhrharmonize.workflow.staging import stage_workflow

    source = make_test_raster(tmp_path / "source" / "image.tif")
    rpc = source.with_suffix(".RPB")
    rpc.write_text("example RPC coefficients")
    config = {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False},
        "import_files": {"plugin": 'import_files', 
            "core:run": True,
            "param:search_glob": str(source),
            "var:raw": "returned:file_path",
            "var:input": "returned:file_path",
            "var:input_name": "expr:$split(var.file_path, '/')[-1]",
            "var:input_stem": r"expr:$replace($split(var.file_path, '/')[-1], /\.[^.]*$/, '')",
            "var:input_dir": r"expr:$replace(var.file_path, /[^\/]+$/, '')",
            "param:create_metadata_json": {"rpc": {"path": "literal:expr:$replace(var.file_path, '.tif', '.{RPB,rpb}')"}},
            "const:output_dir": str(tmp_path / "out"),
            "const:temp_dir": str(tmp_path / "temp"),
        },
        "file_source": {"plugin": 'file_source', 
            "core:run": True,
            "core:requires": "var:rpc",
            "param:input_path": "var:raw",
            "var:result": "expr:const.output_dir & '/result.tif'",
            "param:output_path": "var:result",
        },
    }
    staged, uploads, _ = stage_workflow(
        config,
        config_dir=str(tmp_path),
        remote_output_dir=str(tmp_path / "remote/out"),
        remote_temp_dir=str(tmp_path / "remote/temp"),
        remote_reference_dir=str(tmp_path / "remote/ref"),
    )
    assert Path(uploads[str(source)]).parent == Path(uploads[str(rpc)]).parent
    for local, remote in uploads.items():
        Path(remote).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(local, remote)
    workflow = Workflow(staged)
    assert workflow.nodes[0].requirements == [uploads[str(rpc)]]
    workflow.run()
    assert (tmp_path / "remote/out/result.tif").exists()


def test_hpc_downloads_metadata_checkpoints_for_local_resume(tmp_path, make_test_raster):
    from vhrharmonize.workflow.staging import stage_workflow

    source = make_test_raster(tmp_path / "source.tif")
    config = {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False},
        "import_files": {"plugin": 'import_files', 
            "core:run": True,
            "param:search_glob": str(source),
            "var:raw": "returned:file_path",
            "var:input": "returned:file_path",
            "var:input_name": "expr:$split(var.file_path, '/')[-1]",
            "var:input_stem": r"expr:$replace($split(var.file_path, '/')[-1], /\.[^.]*$/, '')",
            "var:input_dir": r"expr:$replace(var.file_path, /[^\/]+$/, '')",
            "const:output_dir": str(tmp_path / "out"),
        },
        "file_source": {"plugin": 'file_source', 
            "core:run": True,
            "param:input_path": "var:raw",
            "var:result": "expr:const.output_dir & '/result.tif'",
            "param:output_path": "var:result",
            "var:saved_image": "returned:$",
        },
    }
    staged, uploads, downloads = stage_workflow(
        config,
        config_dir=str(tmp_path),
        remote_output_dir=str(tmp_path / "remote/out"),
        remote_temp_dir=str(tmp_path / "remote/temp"),
        remote_reference_dir=str(tmp_path / "remote/ref"),
    )
    checkpoint = str(tmp_path / "out/result.tif.context.json")
    assert checkpoint in downloads
    for local, remote in uploads.items():
        Path(remote).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(local, remote)
    Workflow(staged).run()
    for local, remote in downloads.items():
        Path(local).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(remote, local)
    assert Workflow(config).counts()["file_source"]["loaded"] == 1


def test_hpc_uploads_declared_directory_inputs(tmp_path, make_test_raster):
    from vhrharmonize.workflow.staging import stage_workflow
    from vhrharmonize.slurm import _stage_upload_tree

    source = make_test_raster(tmp_path / "source.tif")
    calibration = tmp_path / "calibration"
    (calibration / "empty").mkdir(parents=True)
    coefficients = calibration / "values.json"
    coefficients.write_text('{"gain": 1}')
    config = {
        "shared": {"plugin": 'shared', "core:run": True, "core:log_to_console": False},
        "import_files": {"plugin": 'import_files', 
            "core:run": True,
            "param:search_glob": str(source),
            "var:raw": "returned:file_path",
            "var:input": "returned:file_path",
            "var:input_name": "expr:$split(var.file_path, '/')[-1]",
            "var:input_stem": r"expr:$replace($split(var.file_path, '/')[-1], /\.[^.]*$/, '')",
            "var:input_dir": r"expr:$replace(var.file_path, /[^\/]+$/, '')",
            "const:output_dir": str(tmp_path / "out"),
        },
        "file_source": {"plugin": 'file_source', 
            "core:run": True,
            "param:input_path": "var:raw",
            "core:requires": str(calibration),
            "var:result": "expr:const.output_dir & '/result.tif'",
            "param:output_path": "var:result",
        },
    }
    remote_root = tmp_path / "remote"
    staged, uploads, _ = stage_workflow(
        config,
        config_dir=str(tmp_path),
        remote_output_dir=str(remote_root / "out"),
        remote_temp_dir=str(remote_root / "temp"),
        remote_reference_dir=str(remote_root / "ref"),
    )
    assert str(calibration) in uploads
    upload_tree = tmp_path / "upload_tree"
    _stage_upload_tree(
        [
            ("uploaded_input_paths", local, remote, str(Path(remote).relative_to(remote_root)))
            for local, remote in uploads.items()
        ],
        str(upload_tree),
    )
    shutil.copytree(upload_tree, remote_root)
    assert (Path(uploads[str(calibration)]) / "empty").is_dir()
    assert (
        Path(uploads[str(calibration)]) / "values.json"
    ).read_bytes() == coefficients.read_bytes()
    Workflow(staged).run()
    assert (remote_root / "out/result.tif").exists()
