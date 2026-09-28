# vhrharmonize: VHR Satellite Imagery Preprocessing Library
[![PyPI version](https://img.shields.io/pypi/v/vhrharmonize.svg)](https://pypi.org/project/vhrharmonize/)
[![PyPI Downloads](https://static.pepy.tech/personalized-badge/vhrharmonize?period=total&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=GREEN&left_text=PyPI+Downloads)](https://pepy.tech/projects/vhrharmonize)
[![Your-License-Badge](https://img.shields.io/badge/License-MIT-green)](#)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15311571.svg)](https://doi.org/10.5281/zenodo.15311571)

---

## Overview

vhrharmonize is an open-source Python library, CLI, and QGIS plugin suite for preprocessing very high resolution (VHR) satellite imagery into analysis-ready products. Ordered YAML plugins define sensor-specific workflows. WorldView and Planet examples use the same runner; metadata mapping and sensor calibration live in the recipe.

---

## Features

- Atmospheric correction workflows ([Py6S](https://github.com/robintw/Py6S) default, [FLAASH](https://github.com/envi-idl/envipyengine) optional backend)
- RPC orthorectification ([Orthority](https://github.com/leftfield-geospatial/orthority))
- Pansharpening ([Orthority](https://github.com/leftfield-geospatial/orthority))
- Optional cloud masking ([OmniCloudMask](https://github.com/DPIRD-DMA/OmniCloudMask))
- Pairwise alignment ([coregix](https://github.com/iosefa/coregix))
- Individual [SpectralMatch](https://github.com/spectralmatch/spectralmatch) functions for coregistration, radiometric matching, seamlines, masking and mosaics
- Generic file discovery, IMD/JSON/XML/YAML parsing, and configurable metadata mapping
- CLI and library-first interfaces
- Automated SLURM processing for distributed High Performance Computing processing

---

## Installation

See the [installation docs](https://vhrharmonize.sefa.ai/getting-started/installation/) for detailed installation instructions or simply install like this: 

```bash
conda create -n vhrharmonize -c conda-forge gdal sixs python=3.11
conda activate vhrharmonize
pip install "vhrharmonize[defaults]"
```

## Getting Started
For an overview of using the library see the [quickstart docs](https://vhrharmonize.sefa.ai/getting-started/quickstart/). Edit a recipe such as [configs/example.worldview.yml](configs/example.worldview.yml), then run:

```bash
vhr workflow --config configs/example.worldview.yml
```
Preview actual pending work before running:

```bash
vhr workflow --config configs/example.worldview.yml --dry-run
```

Version 3 uses unique step names as top-level YAML keys, in execution order. Each step selects its implementation with `plugin: name`. A step without a plugin only assigns context; `plugin: shared` supplies workflow defaults. `param:` keys pass function arguments, `var:` keys assign per-scene variables, `const:` keys assign workflow-wide values, and `core:` keys control execution. `expr:` values use JSONata. The former workflow list, output mappings and WorldView-specific CLI are removed.
See [workflow configuration](https://vhrharmonize.sefa.ai/configuration/workflow-config/) and
[HPC staging](https://vhrharmonize.sefa.ai/cli/hpc/).

The Python API owns execution, defaults and validation. CLIs are generated from
function signatures and docstrings; processing functions live in their plugin modules.

```python
from vhrharmonize import run_workflow, prepare_slurm_plan

counts = run_workflow("configs/example.worldview.yml", dry_run=True)
plan = prepare_slurm_plan("configs/example.hpc.yml", overrides={"run_id": "example"})
```

SpectralMatch algorithms are imported from the installed package. The WorldView recipe documents 19 non-statistics plugins and selects matching, Markov seamlines, masking and a final merge. Alternatives remain commented; function settings are inherited from ordinary shared parameters, and VHR calculates overviews.

Every command uses `vhr <command> --config <yaml>`. Plugin options live in the
workflow YAML. Enable each desired step with `core:run: true`; omitted steps default
to disabled. Set directories through `import_files.param:output_dir` and
`param:temp_dir`; the importer publishes them into `const` by default, or `var` with
`param:directory_scope: var`. The default temp value `sys` creates a system temporary directory.
Plugins declare the context locations core uses for directory features.
`shared.core:output_metadata_path` optionally appends final const/var snapshots to a
YAML-selected JSON file; `core:delete_final_json_first` defaults to true and clears
each destination only on its first write during the run.
Values such as `var:mul` and `const:band_wavelengths_um` read the context, `returned:bytes` selects a function result,
and `expr:var.suffix & '_aligned'` computes a value. JSONata sees one object with `const` and `var` namespaces. Context fields reach function arguments only through explicit `param:` references. Argument precedence is per-plugin `param:`, then shared `param:`, then the function default. Disabled steps do nothing.

```bash
vhr --help
vhr workflow --help
vhr alignment --config configs/example.worldview.yml
vhr hpc-prepare --help
```

A plugin command runs enabled occurrences of that plugin using existing upstream
outputs. See [adding plugins](https://vhrharmonize.sefa.ai/getting-started/adding-plugins/) to register a
Python implementation and automatically expose its `plugin:` value and CLI command.

To use on a super computer (slurm):
```
vhr hpc-prepare --config configs/example.hpc.yml # Create staged HPC/workflow/slurm files
vhr hpc-upload --config configs/1.staged.hpc.yml # Upload required files
vhr hpc-start --config configs/1.staged.hpc.yml # Submit the job
vhr hpc-status --config configs/1.staged.hpc.yml # Print logs and job status
vhr hpc-download --config configs/1.staged.hpc.yml # Download declared outputs

# Other commands:
vhr hpc-stop --config configs/1.staged.hpc.yml # Cancel the submitted job
vhr hpc-close --config configs/1.staged.hpc.yml # Close the SSH multiplex connection

# Or all together:
vhr hpc-prepare --config configs/example.hpc.yml && vhr hpc-upload --config configs/1.staged.hpc.yml && vhr hpc-start --config configs/1.staged.hpc.yml && vhr hpc-status --config configs/1.staged.hpc.yml

# Helpful commands:
# Create preview image and download
conda install gdal
gdal raster resize --size 1%,1% -r average --co TILED=YES --co COMPRESS=DEFLATE "input.tif" "preview.tif"
rsync -avP user@ip:preview.tif .preview.tif
```

In the HPC YAML, `override_download_conflict: validate` (the default) skips existing local files that pass the workflow's output validation and downloads missing or invalid files. Use `no` to skip all existing files, or `yes` to always overwrite. Raster validation uses the same GDAL readability checks as processing steps; non-raster files, including logs, are checked only for existence. The setting is copied into the staged HPC YAML used by `download`; edit that staged file to change the policy for an already prepared run.

## Contributing

We welcome all contributions! We appreciate any feedback, suggestions, or pull requests to improve this project. See the [contributing docs](https://vhrharmonize.sefa.ai/getting-started/contributing/).

---

## License

This project is licensed under the MIT License. See the [LICENSE](https://github.com/cankanoa/vhrharmonize/blob/main/LICENSE) for details.
