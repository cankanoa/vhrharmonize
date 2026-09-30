# CLI commands

All commands use one executable. Workflow, plugin and HPC commands take a configuration file:

```bash
vhr workflow --config configs/example.worldview.yml --dry-run
vhr hpc-prepare --config configs/example.hpc.yml
vhr hpc-start --config configs/1.staged.hpc.yml
vhr cloud_mask --config configs/example.worldview.yml
vhr statistics --input-path history.jsonl --output-path summary.json
```

`workflow` runs the enabled steps in YAML order. A plugin command runs only that
plugin's enabled named steps; other processing steps provide existing paths and
explicit context loads. Missing required upstream outputs are errors. Set
`core:run: true` in the selected steps: choosing a command does not override their
enable/disable settings. Enabled upstream steps retain their configured filename suffixes; YAML-disabled steps are ignored entirely. Multiple named steps selecting the same plugin run in order. Enabled pluginless setup steps also run; the top-level step name does not select a CLI command.

Plugin commands use the same recipe and Python API as the workflow, including
path resolution, const/var context, output validation and reuse. They accept `--dry-run`
and `--config-dir`. Processing parameters belong in the YAML.

```python
from vhrharmonize import run_plugin

run_plugin("cloud_mask", "configs/example.worldview.yml", dry_run=True)
```

Run `vhr --help` to list installed plugin commands. Registration automatically adds
the corresponding subcommand; no CLI wrapper is needed. The built-in names are
`import_files`, `file_source`, `fetch_atmosphere`, `fetch_dem`,
`atmospheric_correction`, `orthorectification`, `pansharpen`, `cloud_mask`,
`alignment`, and `seamline_metadata`, plus the [19 individual SpectralMatch functions](../api/plugins-spectralmatch.md). For example, `vhr global_regression --config recipe.yml` runs enabled steps selecting `plugin: global_regression`. The old `spectralmatch` pipeline command is removed. `plugin: shared` configures core and does not create a command.

HPC commands are `hpc-prepare`, `hpc-upload`, `hpc-start`, `hpc-status`, `hpc-stop`,
`hpc-close`, and `hpc-download`. Prepare uses an HPC recipe; subsequent commands
normally use the generated staged HPC YAML. See [HPC](hpc.md).

Options and help are generated from Python signatures and docstrings. Validation
and execution errors come from those Python APIs. The old `vhr-*` executables and
`--config-path` flag are removed. `python -m vhrharmonize.cli` is equivalent to `vhr`.

`statistics` summarizes recorded core timings using pandas Table Schema JSON.
Use `--per-run` to separate executions or `--run-id ID` to select one. Enable raw
recording with `shared.core:save_statistics_path` in the workflow YAML. It defaults
to `statistics.jsonl` beside the YAML; `shared.core:load_statistics_path` defaults to
the same file and seeds runtime estimates. See
[Statistics](../api/statistics.md) for append behavior and measurement boundaries.

Run through a named step locally with `vhr workflow --config workflow.yml --run-to-step discover_inputs`. Later steps are disabled in a copy; the original YAML is unchanged.
