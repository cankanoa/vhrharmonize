# HPC commands

The CLI is generated from the public Python functions in `vhrharmonize.slurm`. Flags match function arguments, so the configuration flag is `--config`. `hpc-prepare` and `hpc-upload` accept an `--overrides` JSON object or `@filename`. The former individual override flags are removed.

```python
from vhrharmonize import prepare_slurm_plan, upload_slurm_files, start_slurm_job

plan = prepare_slurm_plan("configs/example.hpc.yml", overrides={"run_id": "example"})
# Explicit separate operations when ready to use the cluster:
upload_slurm_files(plan["staged_hpc_file"])
start_slurm_job(plan["staged_hpc_file"])
```

Python and CLI callers use the same defaults and validation. Configuration, transfer and job-control failures propagate as Python exceptions; a CLI failure exits unsuccessfully with that same underlying error.

Use `configs/example.hpc.yml` and `configs/example.slurm.sbatch`. Set `workflow_config` to your ordered recipe and `staged_workflow_file` to the generated remote recipe location. There is no provider selector or upload-key list.

```bash
vhr hpc-prepare --config configs/example.hpc.yml
vhr hpc-upload --config configs/1.staged.hpc.yml
vhr hpc-start --config configs/1.staged.hpc.yml
vhr hpc-status --config configs/1.staged.hpc.yml
vhr hpc-progress --config configs/1.staged.hpc.yml
vhr hpc-download --config configs/1.staged.hpc.yml
```

`hpc-progress` fetches just the public progress snapshot and prints JSON when
available, without changing the staged YAML. Set `shared.core:report_progress: true`
and `shared.core:show_progress: false` to publish snapshots without a console dashboard. `hpc-status`
fetches the same data and renders it with Rich alongside scheduler status and
logs. See the [progress API](../api/progress.md) for Python callbacks and polling.

Preparation uses the same dependency plan as local execution. It stages required input files/directories and reusable outputs, including context checkpoints, and maps declared persistent outputs and their checkpoints for download. A cloud-masked image can therefore be uploaded without uploading or rebuilding its raw processing chain.

The generated `restore_scenes.param:records` contains per-image paths and their const/var contexts. `restore_scenes` is an ordinary registered scene-setting plugin. Original discovery steps are disabled, and arbitrary step names, explicit plugin selections, pluginless setup steps and processing order are preserved. This also supports custom scene-setting plugins; it does not depend on `import_files`. Scene updates, including later imports, whose required preceding processing has not finished are rejected before transfers. Staging materializes known variable assignments and file parameters so the remote job uses the same names without reparsing sensor documents. Adapters opt into `input_hpc_staging_paths`, `output_hpc_staging_paths` and `output_hpc_download_paths` independently; these select function argument names, while YAML still uses ordinary `param:` arguments. Input staging uploads required inputs; output staging rewrites destinations and uploads required cached products; download selection maps persistent products for retrieval without enabling upload. Path resolution is a separate adapter feature. Paths omitted from transfer selections remain the plugin or shared filesystem’s responsibility. Extra files, including paths embedded in compound options, use `core:requires`. Returned scalar/object values and their dependent assignments still resolve at runtime. File arguments selected for staging must be known during planning. Disabled steps contribute no transfers or assignments.

HPC staging rewrites the const/var directory locations registered by plugins, preserving initial workflow-wide constants in `shared` even when no scenes are imported. Distinct per-scene roots receive separate remote directories. The generated restore plugin carries these directory declarations forward. Constants declared by scene steps remain on those steps, with known values frozen once and references to aggregate returns deferred until execution. Both constant and scene return values survive checkpoint download/upload and path rebasing. Remote roots come from `remote_output_dir`, `remote_temp_dir`, and `remote_reference_dir`. The Slurm template runs `vhr workflow --config "$1"`. `prepare` creates local files only; `upload` and `start` remain separate operations.

An explicit `shared.core:output_metadata_path` is also rewritten and added to downloads, including destinations outside the registered output root. Shared destinations remain shared and scene-specific destinations remain separate. With `core:delete_final_json_first: false`, an existing local JSON is included in uploads so the remote run can append its completions; with the default `true`, each remote destination is cleared only on its first write. Metadata destinations must resolve during preparation.

`override_download_conflict` accepts `no`, `yes`, or `validate`. Directory outputs use recursive rsync and preserve per-file conflict handling. `stop` cancels the recorded job; `close` closes its SSH multiplex connection.

SpectralMatch stages require SpectralMatch 1.6+ on the cluster. Explicit raster output lists are staged and downloaded one file at a time; tiled merge directories are transferred recursively. Configure native Dask settings under `shared`: `param:concurrent_processing_backend: dask`, `param:image_threads: null`, and `param:dask_scheduler: [file, /remote/scheduler.json]`. Functions accepting these names inherit them automatically. The final single-file merge explicitly overrides its tile-scheduling settings with null. VHR’s `seamline_metadata` function also needs its `param:dask_scheduler_file` or `param:dask_scheduler_address` set when using the shared Dask backend. Core scene scheduling uses its own `core:` settings and does not configure SpectralMatch automatically.

Multiple `import_files` steps accumulate scenes and metadata. HPC snapshots preserve the combined collection and disable the imports already evaluated locally. New files are not sent through earlier processing stages.

Custom scene IDs remain unchanged during staging and restoration. An ID is rebased only when its value is itself an explicitly staged file path; IDs containing slashes or `~` are otherwise treated as ordinary identifiers.

`core:processing_direction` is preserved in the staged workflow, including per-step overrides. Vertical scene processing uses the same file dependencies and collection barriers on the cluster as locally.
