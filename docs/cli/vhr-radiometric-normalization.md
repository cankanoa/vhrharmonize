# vhr-radiometric-normalization

## Overview

`vhr-radiometric-normalization` runs the shared SpectralMatch wrapper directly. Any `--match-*` argument is forwarded into the upstream SpectralMatch pipeline after the `match_` prefix is stripped.

## Installation

```bash
pip install "vhrharmonize[radiometric-normalization]"
```

## Usage

```bash
vhr-radiometric-normalization \
  --input-image /data/image_a.tif \
  --input-image /data/image_b.tif \
  --output-image /data/normalized.tif
```

- Required: `--input-image` (repeatable), `--output-image`
- Shared runtime controls: `--temp-dir`, `--delete-temp-dir/--no-delete-temp-dir`, `--debug-logs/--no-debug-logs`, `--cache`
- Shared data controls: `--nodata-value`, `--window-size`, `--image-threads`, `--io-threads`, `--tile-threads`, `--calculation-dtype`, `--output-dtype`, `--save-as-cog/--no-save-as-cog`
- Extra SpectralMatch kwargs: `--extra-kwargs-json` and any `--match-*` argument
