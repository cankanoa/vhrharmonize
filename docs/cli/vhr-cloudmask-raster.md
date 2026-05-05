# vhr-cloudmask-raster

## Overview

`vhr-cloudmask-raster` runs OmniCloudMask on one raster, writes a binary mask, and writes a masked raster.

## Installation

```bash
pip install "vhrharmonize[cloud]"
```

## Usage

```bash
vhr-cloudmask-raster \
  --input-raster /data/pansharpened.tif \
  --output-raster /data/pansharpened_cloudmasked.tif \
  --output-mask /data/pansharpened_cloudmask.tif \
  --buffer-pixels 3
```

- Required: `--input-raster`
- Optional outputs: `--output-raster`, `--output-mask`
- Mask values: `--nodata-value`, `--cloud-classes`, `--buffer-pixels`
- Band mapping: `--red-band-index`, `--green-band-index`, `--nir-band-index`
- Model kwargs passthrough: `--omnicloud-kwargs-json`
