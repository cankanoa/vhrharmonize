# vhr-flaash

## Overview

`vhr-flaash` runs the shared FLAASH wrapper from a supplied parameter object and writes the executed parameter dump.

## Installation

```bash
pip install "vhrharmonize[flaash]"
```

## Usage

```bash
vhr-flaash \
  --params-json-file outputs/flaash_params.json \
  --output-params-path outputs/flaash_params_used.json \
  --envi-engine-path /path/to/taskengine.exe
```

- Required: `--output-params-path`, `--envi-engine-path`
- Provide exactly one parameter source: `--params-json` or `--params-json-file`
- Optional path behavior: `--convert-paths-for-windows`, `--delete-output-before-run`
- The parameter object must include `INPUT_RASTER.url` and `OUTPUT_RASTER_URI`
