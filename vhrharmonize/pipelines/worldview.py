"""WorldView pipeline entrypoints built on packaged CLI implementation."""

from typing import Dict, List, Optional

from vhrharmonize.cli import worldview as worldview_cli


def run_worldview_pipeline(config_yaml: Optional[str] = None, extra_args: Optional[List[str]] = None) -> int:
    """Run the WorldView pipeline through the canonical CLI entrypoint."""
    argv: List[str] = []
    if config_yaml:
        argv.extend(["--config-yaml", config_yaml])
    if extra_args:
        argv.extend(extra_args)
    return worldview_cli.main(argv)


def run_worldview_from_config(config: Dict) -> int:
    """Run with an in-memory config mapping (translated to CLI args)."""
    argv: List[str] = []
    for key, value in config.items():
        arg = f"--{key.replace('_', '-')}"
        if isinstance(value, list):
            for item in value:
                argv.extend([arg, str(item)])
        elif isinstance(value, bool):
            if value:
                argv.append(arg)
        elif value is not None:
            argv.extend([arg, str(value)])
    return worldview_cli.main(argv)
