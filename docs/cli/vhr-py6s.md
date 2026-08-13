# vhr-py6s

## Overview

`vhr-py6s` runs Py6S-only atmospheric correction on discovered WorldView scenes. It can optionally orthorectify the corrected raster afterward.

## Installation

```bash
pip install "vhrharmonize[py6s]"
```

This installs Py6S and a packaged 6S executable (from [6s-bin](https://github.com/brianschubert/6s-bin)) as Python ependencies. If the packaged executable is unavailable or fails on your platform, run `pip install Py6S`, then follow the [Py6S installation guide](https://py6s.readthedocs.io/en/latest/installation.html) to install 6S manually. Pass its executable path using `--py6s-executable` if its not included in sys.path.

## Usage

```bash
vhr-py6s \
  --input-dir /data/worldview_batch \
  --output-dir /data/py6s_only \
  --output-suffix _py6s \
  --py6s-atmosphere-profile user \
  --py6s-aerosol-profile maritime \
  --py6s-auto-atmos-source nasa_power
```

```bash
vhr-py6s \
  --input-dir /data/worldview_batch \
  --output-dir /data/py6s_only \
  --output-suffix _py6s \
  --epsg 6635 \
  --dem-file-path /data/dem.tif
```

- Discovery and outputs: `--config-yaml`, `--input-dir`, `--filter-basename`, `--output-dir`, `--output-suffix`
- Optional orthorectification after Py6S: `--epsg`, `--dem-file-path`, `--keep-intermediate-py6s`
- Py6S atmosphere controls: `--py6s-atmosphere-profile`, `--py6s-aerosol-profile`, `--py6s-aot550`, `--py6s-visibility`, `--py6s-water-vapor`, `--py6s-ozone`
- Output scaling: `--py6s-output-scale-factor`, `--py6s-output-dtype`, `--py6s-executable`
- Calibration toggles: `--py6s-use-imd-radiance-calibration/--no-py6s-use-imd-radiance-calibration`, `--py6s-use-worldview-gain-offset-adjustment/--no-py6s-use-worldview-gain-offset-adjustment`
- Auto atmosphere controls: `--py6s-auto-atmos-source`, `--py6s-auto-atmos-grid-size`, `--py6s-auto-atmos-search-days`, `--py6s-auto-atmos-timeout-s`, `--py6s-auto-atmos-power-endpoint`
