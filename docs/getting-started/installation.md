# Installation

## Install with PyPI

1. Create and activate a conda environment.

```bash
conda create -n vhrharmonize -c conda-forge py6s sixs gdal python=3.11
conda activate vhrharmonize
```

2. Install the package with default dependencies.

```bash
pip install vhrharmonize[default]
```

3. Install specific dependencies as needed.

```bash
pip install "vhrharmonize[cloud]"
pip install "vhrharmonize[py6s]"
pip install "vhrharmonize[flaash]"
pip install "vhrharmonize[orthorectification]"
pip install "vhrharmonize[pansharpen]"
pip install "vhrharmonize[align]"
pip install "vhrharmonize[radiometric-normalization]"
pip install "vhrharmonize[docs]"
pip install "vhrharmonize[all]"
```

## Install from source

1. Create and activate a conda environment.

```bash
conda create -n vhrharmonize -c conda-forge py6s sixs gdal python=3.11
conda activate vhrharmonize
```

2. Clone the repository.

```bash
git clone https://github.com/cankanoa/vhrharmonize.git
cd vhrharmonize
```

3. Install the default packages.

```bash
pip install -e '.[default]'
```

4. Install specific dependencies as needed.

```bash
pip install -e ".[cloud]"
pip install -e ".[py6s]"
pip install -e ".[flaash]"
pip install -e ".[orthorectification]"
pip install -e ".[pansharpen]"
pip install -e ".[align]"
pip install -e ".[radiometric-normalization]"
pip install -e ".[docs]"
pip install -e ".[all]"
```

5. Verify the entry points if desired.

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
