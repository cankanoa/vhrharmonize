# vhrharmonize: VHR Satellite Imagery Preprocessing Library

[![License: MIT](https://img.shields.io/badge/License-MIT-green)](#)

> [!IMPORTANT]
> This library is under active development and may change.


---

## Overview

`vhrharmonize` is an open-source Python library and CLI suite for preprocessing very high resolution (VHR) satellite imagery into analysis-ready products.

Current tested/supported end-to-end workflow: **WorldView**.  
Additional providers and sensors (for example, Planet) are being added.

---

## Features

- Atmospheric correction workflows (Py6S default, FLAASH optional backend)
- RPC orthorectification
- Pansharpening
- Optional cloud masking with OmniCloudMask
- Tile utilities and tile-aware pairwise alignment (elastix)
- WorldView scene discovery, IMD parsing, and standardized metadata mapping
- CLI and library-first interfaces

---

## Installation

### 1. Requirements

- Python `>=3.9`
- GDAL-compatible environment
- `6S` executable for Py6S workflows (default)
- ENVI Task Engine only for FLAASH workflows

### 2. Create environment (example)

```bash
conda create -n vhrharmonize python=3.11 -c conda-forge
conda activate vhrharmonize
```

### 3. Install package

```bash
git clone https://github.com/cankanoa/vhrharmonize.git
cd vhrharmonize
pip install -e .
```

Optional extras by step:

```bash
pip install -e ".[fetch-atmosphere]"            # requests
pip install -e ".[cloud]"                       # omnicloudmask + scipy + gdal
pip install -e ".[py6s]"                        # Py6S + tqdm
pip install -e ".[flaash]"                      # envipyengine + gdal + tqdm
pip install -e ".[orthorectification]"          # orthority + gdal + pyproj
pip install -e ".[pansharpen]"                  # orthority
pip install -e ".[align]"                       # itk-elastix
pip install -e ".[radiometric-normalization]"   # spectralmatch
pip install -e ".[docs]"
pip install -e ".[all]"
```

For conda users, install 6S separately when using Py6S:

```bash
mamba install -n vhrharmonize -c conda-forge py6s sixs
```

## Getting Started

For usage docs:

- `docs/getting-started/quickstart.md`
- `docs/cli/overview.md`

Primary CLI commands:

```bash
vhr-worldview --help
vhr-fetch-modis-water-vapor --help
vhr-flaash --help
vhr-cloudmask-raster --help
vhr-pansharpen-orthos --help
vhr-align-image-pair --help
vhr-orthorectification --help
vhr-radiometric-normalization --help
vhr-py6s --help
```

Docs site:

```bash
python -m mkdocs serve -f docs/mkdocs.yml
```

## Contributing

We welcome all contributions! To get started:
1. Fork the repository and create a new feature branch.
2. Make your changes.
3. Open a Pull Request against the main repository.

We appreciate any feedback, suggestions, or pull requests to improve this project.

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
