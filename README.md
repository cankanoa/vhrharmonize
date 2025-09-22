# vhrharmonize: automated open source library to preprocess and organize satellite imagery

[![License: MIT](https://img.shields.io/badge/License-MIT-green)](#)

> [!IMPORTANT]
> This library is under active development and may change.


---

## Overview

vhrharmonize is an open-source Python library designed to assist in preprocessing satellite imagery into analysis-ready data products. It provides ready-to-use functions for performing atmospheric correction, orthorectification, pansharpening, and Relative Radiometric Normalization (spectralharmonize). Alongside these core tools, it offers a standardized structure and automated pipeline to help process large volumes of imagery consistently and efficiently.

The goal of vhrharmonize is to provide an automated, modular workflow for preprocessing satellite imagery—regardless of sensor or processing level—into spectrally and geometrically consistent outputs. With the goal of supporting a wide range of inputs for now it just supportss WorldView imagery.

---

## Features

- **End-to-End Processing:** Automated python scripts help to automate processing.
- **Atmospheric Correction:** Uses ENVI's FLAASH to convert raw DNs into reflectance.
- **GCP-Based Orthorectification:** Supports refined RPC correction using custom or QGIS GCPs.
- **Pansharpening:** Fuses multispectral and panchromatic images for sharper detail.
- **Batch-Friendly:** Automatically processes all images in a directory tree.
- **Customizable:** Allows image filtering, metadata overrides, and DEM selection.

---

## Installation

### 1. Manual System Requirements

- **ENVI ≥ 5.7** (for FLAASH and task engine support; with atmospheric correction module)
- Optional: **QGIS** (for GCP creation)

### 2. Library System Requirements
- **GDAL ≥ 3.4**
- **Python ≥ 3.8**

```bash
conda create -n vhrharmonize python=3.10 "gdal=3.10.2" -c conda-forge
conda activate spectralmatch
```

### 3. Install the package

```bash
git clone https://github.com/yourusername/vhrharmonize.git
cd vhrharmonize
pip install -e .
```

## Getting Started

To see how to run the full processing pipeline, start with the example script at 'docs/examples/worldview_b1.py'.

## Contributing

We welcome all contributions! To get started:
1. Fork the repository and create a new feature branch.
2. Make your changes.
3. Open a Pull Request against the main repository.

We appreciate any feedback, suggestions, or pull requests to improve this project.

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE.md) for details.