# vhr-orthorectification

## Overview

`vhr-orthorectification` exposes the orthorectification helper utilities: direct RPC orthorectification, QGIS GCP text to CSV, and QGIS GCP text to Orthority-compatible GeoJSON.

## Installation

```bash
pip install "vhrharmonize[orthorectification]"
```

## Usage

```bash
vhr-orthorectification orthorectify \
  --input-image-path /data/input.tif \
  --output-image-path /data/output_ortho.tif \
  --dem-image-path /data/dem.tif \
  --output-epsg 6635
```

```bash
vhr-orthorectification qgis-gcps-to-csv \
  --input-gcp-path /data/qgis_gcps.txt \
  --output-csv-path /data/gcps.csv \
  --output-epsg 6635
```

```bash
vhr-orthorectification qgis-gcps-to-geojson \
  --input-image-path /data/input.tif \
  --qgis-gcp-file-path /data/qgis_gcps.txt \
  --file-name scene_name \
  --dem-file-path /data/dem.tif \
  --output-geojson-path /data/gcps.geojson
```

- `orthorectify`: `--input-image-path`, `--output-image-path`, `--dem-image-path`, `--output-epsg`, optional `--gcp-geojson-file-path`, `--output-nodata-value`, `--dtype`, `--output-resolution`, `--output-resolution-x`, `--output-resolution-y`
- `qgis-gcps-to-csv`: `--input-gcp-path`, `--output-csv-path`, optional `--output-epsg`
- `qgis-gcps-to-geojson`: `--input-image-path`, `--qgis-gcp-file-path`, `--file-name`, `--dem-file-path`, `--output-geojson-path`, optional `--force-positive-pixel-values`
