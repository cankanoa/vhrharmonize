"""Unified commands generated from the Python workflow and HPC APIs."""
from collections.abc import Mapping
from pathlib import Path

from vhrharmonize import slurm
from vhrharmonize.workflow.api import run_plugin, run_workflow
from vhrharmonize.workflow.registry import plugin_names
from .functions import commands_cli


def _plugin_command(name):
    def command(config: str | Path | Mapping, *, config_dir: str | None = None,
                dry_run: bool = False):
        """Run enabled PLUGIN steps from a workflow recipe.

        Args:
            config: Workflow YAML filename.
            config_dir: Optional base directory for relative paths.
            dry_run: Inspect the plan without processing outputs.
        """
        return run_plugin(name, config, config_dir=config_dir, dry_run=dry_run)
    command.__doc__ = command.__doc__.replace("PLUGIN", name)
    return command


def main(argv=None):
    commands = {
        "workflow": run_workflow,
        "hpc-prepare": slurm.prepare_slurm_plan,
        "hpc-upload": slurm.upload_slurm_files,
        "hpc-start": slurm.start_slurm_job,
        "hpc-status": slurm.update_status_slurm_file,
        "hpc-stop": slurm.stop_slurm_job,
        "hpc-close": slurm.close_hpc_connection,
        "hpc-download": slurm.download_slurm_outputs,
    }
    for name in plugin_names():
        if name in commands:
            raise ValueError(f"Plugin name {name!r} conflicts with a built-in vhr command")
        commands[name] = _plugin_command(name)
    return commands_cli(commands, argv, prog="vhr", description="Run workflows, plugins and HPC jobs.")
