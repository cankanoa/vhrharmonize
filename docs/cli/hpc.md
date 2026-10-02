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

Use `configs/example.hpc.yml` and `configs/example.slurm.sbatch`. Set `workflow_config` to your ordered recipe and `staged_workflow_file` to the generated remote recipe location. Use `path_mappings` to associate workflow variables with remote directories.

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
fetches the same data and renders it with `prompt_toolkit` alongside scheduler status and
logs. See the [progress API](../api/progress.md) for Python callbacks and polling.

Preparation uses the same requested-output dependency plan as local execution. It uploads required inputs, valid reusable products, and explicitly loaded context files. Generated outputs are downloaded when requested through `core:require_outputs` and enabled by the plugin's download declaration.

The source YAML is preserved. Preparation and remote copies change only selected values or entries, retaining comments, quotes, ordering, and expressions. There is no `restore_scenes` section, embedded scene collection, frozen `staged_*` variable table, or automatically discovered `.context.json` sidecar.

## Upload progress

`hpc-upload` uses the same terminal display as workflow progress, with a total row
and one row per step and mapped variable (for example, `alignment` / `var:relative_output_dir`).
Each row shows completed/total files, processed/total bytes, speed, status, ETA, and
a byte-based bar at the right. Green means transferred, purple means already
current on HPC, and gray means remaining. File names are not printed. Directories
are counted by their contained files; duplicate destinations count once.

Set `show_progress: true` in the **HPC YAML** (the default). Redirected output gets
one final table; interactive terminals update the same box using native scrollback.
Disable it for a short completion message:

```bash
vhr hpc-upload --config configs/1.staged.hpc.yml --overrides '{"show_progress": false}'
```

Transfers still use batched rsync. Speed is the average logical file bytes processed
per second, including rsync's reconstruction of changed files, rather than measured
SSH network traffic. ETA uses that rate and remaining bytes; it shows `TBD` until
there is a rate. Unchanged files are credited when their batch succeeds. See
[rsync's progress semantics](https://download.samba.org/pub/rsync/rsync.1#opt--progress).
Ctrl-C cancels the transfer and restores the terminal; SSH can still read terminal
input for authentication.

New preparation records `upload_groups` in the staged HPC YAML for the step/root
labels. Older staged files still upload with `inputs` / `controls` groups; run
`hpc-prepare` again to get step attribution. Python callers can consume the same
data independently of the display using
`upload_slurm_files(config, progress_callback=callback)`; see the
[upload callback API](../api/progress.md#upload-progress).

## Local preparation through a named step

```yaml
run_to_step_before_prepare: discover_inputs # null prevents local plugin calls
```

Core writes a `.prepare.yml` copy beside the staged workflow. Steps after the cutoff receive `core:skip_plugin_call: true`; existing enable/disable choices remain intact. Enabled steps through the cutoff run locally, including explicit context saves. Their resulting state is handed directly to staging. The remote copy restores that state and disables the completed prefix. For a separate local run, use `vhr workflow --config workflow.yml --run-to-step discover_inputs`.

With `debug_logs: true` in the HPC YAML, preparation logs its start, the local
execution cutoff, and the start of staging. Discovery-only local runs keep their
logs and progress snapshots but omit the empty workflow dashboard.

The cutoff names a workflow step, independently of its `plugin:` selection.
Both sensor examples use `discover_inputs` for that step; changing its plugin
does not require changing the HPC cutoff. Discovery-specific path rewrites use
the optional [plugin staging interface](../getting-started/adding-plugins.md#custom-discovery-staging).

## Directory, file, and file-list mappings

Map path-valued workflow references to remote folders:

```yaml
remote_work_dir: ~/koa_scratch/vhrharmonize/run_{run_id}
path_mappings:
  var:current_image_paths: "expr:'~/koa_scratch/vhrharmonize/inputs/' & var.scene_id"
  var:pan: "expr:'~/koa_scratch/vhrharmonize/inputs/' & var.scene_id"
  var:metadata_path: "expr:'~/koa_scratch/vhrharmonize/inputs/' & var.scene_id"
  var:mul_companions: "expr:'~/koa_scratch/vhrharmonize/inputs/' & var.scene_id"
  var:pan_companions: "expr:'~/koa_scratch/vhrharmonize/inputs/' & var.scene_id"
  const:output_dir: ~/koa_scratch/vhrharmonize/run_{run_id}/products
  var:relative_output_dir: ~/koa_scratch/vhrharmonize/run_{run_id}/products
  const:temp_dir: ~/koa_scratch/vhrharmonize/run_{run_id}/temp
  const:dem_path: /shared/reference
  const:reference_path: /shared/reference
```

| Selected value | Remote placement |
| --- | --- |
| Directory | Replace the root and preserve the relative layout of its contents. |
| File | Place the file directly in the destination folder, retaining its filename. |
| Flat list of paths | Apply the corresponding rule to each element; empty companion lists are valid. |

For example, a directory mapping preserves `scenes/P004/image.tif`, while an
individual file mapping to the scene's folder puts that image at `inputs/P004/image.tif`. There are
no generated hashes or extra root basenames. Several file mappings can share a
folder, so the DEM and alignment reference can both use `/shared/reference`.
Two different files resolving to the same complete remote filename are an error;
repeated references to the same file transfer it only once. Explicit file mappings
take precedence over enclosing directory mappings; nested directories use the
most specific root.

Mappings use the first resolved values across scenes. `var:current_image_paths`
therefore selects the initial imported images, and later processing assignments
keep their normal expressions. The initial assignment retains its string/list
shape. For flattened imports, staging selectively replaces the source glob with
the selected remote image paths and adjusts companion/metadata rules when needed.
Portable expressions remain unchanged; irregular per-scene file associations use
small path lookups. Prepared metadata is stored once in a separate JSON control
file rather than embedded in the workflow YAML.

Destination expressions use the existing JSONata syntax: quote literal path text,
join it with `&`, and reference fields as `var.scene_id` or `const.name` inside
the expression. `expr:~/inputs/var:scene_id/` is not valid expression syntax.
An expression is evaluated in the original context of each selected path value,
before paths are relocated. All entries in a scene's file list use that scene's
destination. It must resolve to a nonempty directory string. `~` remains unexpanded
until use on HPC, and `{run_id}` is substituted before expression evaluation.

The WorldView example maps MUL images, PAN images, both RPC companion lists, and
the source IMD metadata paths into `inputs/<scene_id>/`. Set the local discovery
pattern directly in `param:search_glob`; no `const:input_root` is needed. The two reference
files have separate constants and share a remote reference folder rather than
being copied for each scene. Uploads under the common input directory remain batched.

`run_to_step_before_prepare` is the local plugin-call boundary. Preparation runs
one workflow: steps through the named cutoff may call plugins, and later steps
receive `core:skip_plugin_call: true` in the generated `.prepare.yml`. They remain
in the dependency plan, but do not call plugins, compute overviews, or save
processing results locally. `null` permits no local plugin calls. The alternate
spelling `run_to_stop_before_prepare` is accepted; conflicting settings are rejected.

The prepared constants and initialized `var` records are handed directly to
staging and written to one relocated `prepared.json` control file. The remote
recipe restores that state and disables the locally completed prefix. Discovery
is not repeated, even when the original recipe does not load a saved context.
Explicit context files remain supported. Preparation retains intermediates needed
for upload instead of deleting them during cleanup.

Paths derivable from configuration and prepared values are checked normally.
Unresolved paths remain needed, with their expressions intact for remote execution.
Dynamically named outputs require a known mapped output directory for downloading.
An unresolved external input must have explicit `core:requires` dependencies or
be resolved by advancing the local cutoff. If record discovery is still pending,
advance the cutoff through that producer so staging can determine per-record inputs.
Constants-only workflows need no record producer and can use a null cutoff.

Directory values are ordinary YAML assignments; the importer has no temp/output
directory parameters. `core:cleanup_dirs` selects the corresponding remapped roots
for cleanup during remote execution. Optional explicit context files under
`const:temp_dir` use its mapping too. If a load path depends on `const:temp_dir`,
define that constant in `shared` or an earlier step so it exists before loading.
The examples do not need `project_root`. Additional context files outside mapped
directories need their own mapping; statistics and generated workflow/sbatch
control files are staged separately.

Required scene/context paths must be covered; missing coverage is an error.
Mapping a directory to itself, or a file to its existing parent folder, keeps its
path but does not imply shared storage or automatically disable transfers.

Mappings choose placement, not which files to transfer. A mapped directory is uploaded recursively only when the directory itself is a required input. Declared plugin inputs/outputs and `core:requires` supply file dependencies; files read by imports supply discovery dependencies. Unused processing steps contribute no transfers. The older three-directory Python/HPC interface remains available for existing callers; new configurations should use explicit root mappings.

## Explicit metadata and discovered products

Use [context operations](../configuration/workflow-config.md#explicit-context-files) to save selected fields and load them before planning. A missing load file warns only when console logs are enabled and otherwise continues; it never triggers a producer or metadata preprocessing automatically. Required missing values can still prevent a later function from running.

An import with an explicit context load covering its declared scene assignments can reuse those records without scanning the original files. Otherwise discovery runs normally and its source images and metadata are included in transfers. Loaded snapshots are copied for remote path rebasing; the local JSON is untouched. Constants are stored once and scene values are keyed by scene ID. Downloads preserve the resolved values in JSON; they do not implicitly translate returned paths back to local paths. Use stable scene IDs and select portable metadata fields for snapshots that will be loaded on both hosts.

```yaml
import_cloudmasked:
  plugin: import_files
  core:run: true
  core:satisfies: {cloud_mask: output_raster_path}
  param:search_glob: /data/project/products/*_cloudmasked.tif
  # Supply a scene-ID rule matching the other imports.
  var:cloudmasked_file: returned:file_path
```

`core:satisfies` associates each imported image with the named step's declared output parameter. It does not supply unrelated metadata or mark all of that step's outputs complete. Existing output validation/reuse rules still apply.

The Slurm template runs `vhr workflow --config "$1"`. Preparation creates local artifacts and may execute the explicitly requested preparation cutoff; uploading and submitting remain separate commands.

`core:save_statistics_path` is a separate remote append destination
registered for download, so staging cannot replace remote history with an older local
copy. See [statistics on HPC](../api/statistics.md#hpc-files) and
[saved file formats](../saved-file-formats.md).
