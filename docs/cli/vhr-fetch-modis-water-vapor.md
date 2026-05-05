# vhr-fetch-modis-water-vapor

## Overview

`vhr-fetch-modis-water-vapor` scans WorldView scenes, looks up MODIS atmosphere data in Google Earth Engine, and writes JSON and CSV reports.

## Installation

```bash
pip install "vhrharmonize"
```

## Usage

```bash
vhr-fetch-modis-water-vapor \
  --input-dir /data/worldview_batch \
  --ee-project your-ee-project \
  --output-json outputs/modis_water_vapor_results.json \
  --output-csv outputs/modis_water_vapor_results.csv
```

- Required: `--input-dir`
- Scene filtering: `--filter-basename`
- Earth Engine setup: `--env-file`, `--ee-project`, `--authenticate`
- Search controls: `--hours-window`, `--terra-collection`, `--aqua-collection`, `--band-name`, `--aod-band-candidates`, `--reduce-scale-m`
- Scaling controls: `--aod-scale-factor`, `--modis-scale-factor`, `--modtran-baseline-water-vapor`
- Outputs: `--output-json`, `--output-csv`
