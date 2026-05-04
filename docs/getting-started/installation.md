# Installation

## System Requirements

- Python `>=3.9`
- A GDAL-compatible environment for raster-heavy steps
- `6S` executable only when using Py6S workflows
- ENVI Task Engine only when using FLAASH workflows

The base install now keeps only widely shared Python dependencies. Step-specific
packages are installed through extras.

Recommended conda setup for Py6S + 6S:

```bash
mamba install -n vhrharmonize -c conda-forge py6s sixs
```

## Install in Editable Mode

```bash
git clone https://github.com/cankanoa/vhrharmonize.git
cd vhrharmonize
pip install -e .
```

With cloud masking extras:

```bash
pip install -e ".[cloud]"
```

With external atmosphere fetch extras:

```bash
pip install -e ".[fetch-atmosphere]"
```

With Py6S atmospheric correction extras:

```bash
pip install -e ".[py6s]"
```

With FLAASH extras:

```bash
pip install -e ".[flaash]"
```

With orthorectification extras:

```bash
pip install -e ".[orthorectification]"
```

With pansharpening extras:

```bash
pip install -e ".[pansharpen]"
```

With radiometric normalization extras:

```bash
pip install -e ".[radiometric-normalization]"
```

With `coregix` alignment extras:

```bash
pip install -e ".[align]"
```

With everything:

```bash
pip install -e ".[all]"
```

With docs extras:

```bash
pip install -e ".[docs]"
```

## Verify CLI Entry Points

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

## Build Docs Locally

From repository root:

```bash
mkdocs serve -f docs/mkdocs.yml
```

Or build static output:

```bash
mkdocs build -f docs/mkdocs.yml
```
