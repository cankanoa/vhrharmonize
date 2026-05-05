# vhr-align-image-pair

## Overview

`vhr-align-image-pair` aligns one raster to another with `coregix` and writes a single aligned output raster.

## Installation

```bash
pip install "vhrharmonize[align]"
```

## Usage

```bash
vhr-align-image-pair \
  --moving-image /data/moving.tif \
  --fixed-image /data/fixed.tif \
  --output-image /data/aligned.tif \
  --split-factor 2 \
  --trim-edge-invalid \
  --edge-trim-depth 8 \
  --edge-trim-invalid-below -3000
```

- Required: `--moving-image`, `--fixed-image`, `--output-image`
- Band selection: `--band-index`, `--moving-band-index`, `--fixed-band-index`
- Nodata and overlap: `--moving-nodata`, `--fixed-nodata`, `--output-nodata`, `--min-valid-fraction`
- Temp and logging: `--temp-dir`, `--keep-temp-dir`, `--log-to-console`
- Alignment behavior: `--use-edge-proxies/--no-use-edge-proxies`, `--split-factor`, `--clip-fixed-to-moving/--no-clip-fixed-to-moving`, `--output-on-moving-grid/--no-output-on-moving-grid`, `--solve-resolution`
- Edge trimming: `--trim-edge-invalid/--no-trim-edge-invalid`, `--edge-trim-depth`, `--edge-trim-detection-band-index`, `--edge-trim-invalid-below`, `--edge-trim-invalid-above`
- Masking: `--enforce-mutual-valid-mask/--no-enforce-mutual-valid-mask`
