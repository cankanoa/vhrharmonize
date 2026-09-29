# HPC status and workflow progress

`shared.core:show_progress` defaults to `true`, enabling snapshots and terminal output for HPC runs. To publish snapshots without the job's console display, set `core:show_progress: false` and `core:report_progress: true` before preparing and uploading the workflow. Alongside Slurm status and logs, `vhr hpc-status --config <staged.hpc.yml>` prints the latest saved workflow dashboard once, including separate Unused, Loaded, Done, Run and All columns with percentages of All, per-step ETA, active operations and recent messages.

The job updates `<remote_workflow_config>.progress.json` atomically while it runs, even without an interactive terminal. The status command shows the snapshot timestamp and ignores snapshots from a different Slurm job ID. Before the workflow starts, or when progress is disabled, ordinary Slurm status and logs remain available. The fetched snapshot is saved under `workflow_progress` in the staged HPC YAML.

Apps can call `get_slurm_progress(staged_hpc_file)` to fetch only the data, without
printing or modifying the staged file. `vhr hpc-progress --config <staged.hpc.yml>`
exposes that snapshot as JSON. The `prompt_toolkit` renderer consumes the same public [progress API](progress.md)
for local updates and remote status; the current snapshot schema is version 2.

::: vhrharmonize.slurm
