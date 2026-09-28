# vhr workflow

This command is generated from the public `run_workflow` function's signature and docstring. `--config`, `--config-dir` and `--dry-run` correspond directly to its Python arguments. Defaults and validation live in the Python API.

```python
from vhrharmonize import load_workflow, run_workflow

counts = run_workflow("configs/example.worldview.yml", dry_run=True)
run_workflow("configs/example.worldview.yml")

# An in-memory config works too; choose its relative-path base explicitly.
workflow = load_workflow(config_dict, config_dir="/data/project")
workflow.plan()
workflow.run()
context = workflow.records[0]["context"]  # {"const": {...}, "var": {...}}
scene_variables = context["var"]
constants = workflow.context["const"]
```

Run an ordered, sensor-neutral plugin recipe:

```bash
vhr workflow --config configs/example.worldview.yml --dry-run
vhr workflow --config configs/example.worldview.yml
vhr workflow --config configs/example.planet.yml
```

`--dry-run` prints existing output counts and actual pending work without running ordinary processing steps or cleanup. Scene-setting functions run to discover their records; their own import metadata saves may occur. Final processing JSON is only written by a real run. All step options live in the YAML. See [workflow configuration](../configuration/workflow-config.md) for paths, metadata mapping, plugin registration and resume behavior.

The former `vhr-worldview` command and flat configurations are not supported in version 3.
