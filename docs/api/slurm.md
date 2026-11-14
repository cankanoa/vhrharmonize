# HPC status and workflow progress

`shared.core:show_progress` defaults to `true`, enabling the workflow terminal
summary in local and HPC runs. Progress is kept in memory; no progress JSON files
are written or fetched. `vhr hpc-status --config <staged.hpc.yml>` reports scheduler
status, upload results, and job logs. Applications running a workflow can receive
updates through the [progress callback API](progress.md).

::: vhrharmonize.slurm
