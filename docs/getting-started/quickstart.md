# Quickstart

Install the package with the preprocessing extras you need, then edit `configs/example.worldview.yml` or `configs/example.planet.yml`. Set the import glob, metadata paths/mappings, output roots, reference images and enabled plugins.

```bash
vhr workflow --config configs/example.worldview.yml --dry-run
vhr workflow --config configs/example.worldview.yml
```

The ordered YAML controls all step parameters and path connections. Existing downstream outputs automatically satisfy dependencies. See [workflow configuration](../configuration/workflow-config.md) and [HPC execution](../cli/hpc.md).
