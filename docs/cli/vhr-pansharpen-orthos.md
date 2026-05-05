# vhr-pansharpen-orthos

## Overview

`vhr-pansharpen-orthos` pansharpens an orthorectified multispectral raster with an orthorectified panchromatic raster.

## Installation

```bash
pip install "vhrharmonize[pansharpen]"
```

## Usage

```bash
vhr-pansharpen-orthos \
  --mul-ortho /data/mul_ortho.tif \
  --pan-ortho /data/pan_ortho.tif \
  --output /data/pansharpened.tif \
  --nodata-value -9999
```

- Required: `--mul-ortho`, `--pan-ortho`, `--output`
- Optional output nodata: `--nodata-value`
