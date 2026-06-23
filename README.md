# vhrharmonize: VHR Satellite Imagery Preprocessing Library

[![License: MIT](https://img.shields.io/badge/License-MIT-green)](#)

> [!IMPORTANT]
> This library is under active development and may change.


---

## Overview

`vhrharmonize` is an open-source Python library and CLI suite for preprocessing very high resolution (VHR) satellite imagery into analysis-ready products. Current supported end-to-end workflows: **WorldView-3 B1 imagery**. Additional providers and sensors (for example, Planet) will be added.

---

## Features

- Atmospheric correction workflows ([Py6S](https://github.com/robintw/Py6S) default, [FLAASH](https://github.com/envi-idl/envipyengine) optional backend)
- RPC orthorectification ([Orthority](https://github.com/leftfield-geospatial/orthority))
- Pansharpening ([Orthority](https://github.com/leftfield-geospatial/orthority))
- Optional cloud masking ([OmniCloudMask](https://github.com/DPIRD-DMA/OmniCloudMask))
- Pairwise alignment ([coregix](https://github.com/iosefa/coregix))
- Relative Radiometric Normalization ([spectralmatch](https://github.com/spectralmatch/spectralmatch))
- WorldView scene discovery, IMD parsing, and standardized metadata mapping
- CLI and library-first interfaces

---

## Installation

See [docs/getting-started/installation.md](docs/getting-started/installation.md) for detailed installation instructions or simply install like this: 

```bash
conda create -n vhrharmonize -c conda-forge py6s sixs gdal python=3.11
conda activate vhrharmonize
pip install vhrharmonize[defaults]
```
## Getting Started
For an overview of using the library see [docs/getting-started/quickstart.md](docs/getting-started/quickstart.md). The CLI can be usd by passing in arguments from a yaml file like this one [docs/configs/example.worldview.yml](docs/configs/example.worldview.yml) and running:

```bash
vhr-worldview --config-yaml example.worldview.yml
```
Or pass in arguments directly from the command line:

```bash
vhr-worldview \
  --input-file-glob "/data/worldview/**/*.TIF" \
  --output-dir ../../processed \
  --run-alignment \
  --alignment-fixed-image /data/reference.tif
```

For detailed arguments use:

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

To use on a super computer (slurm):
```
vhr-hpc prepare --config configs/example.hpc.yml
vhr-hpc upload --config configs/1.staged.hpc.yml
vhr-hpc start --config configs/1.staged.hpc.yml
vhr-hpc status --config configs/1.staged.hpc.yml
vhr-hpc download --config configs/1.staged.hpc.yml
```

## Contributing

We welcome all contributions! We appreciate any feedback, suggestions, or pull requests to improve this project. See [docs/getting-started/contributing.md](docs/getting-started/contributing.md).

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
