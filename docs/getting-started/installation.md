# Installation

## System Requirements

- Python `>=3.9`
- GDAL-compatible environment
- `6S` executable for Py6S workflows (default atmospheric method)
- ENVI Task Engine only if using FLAASH workflows

Optional:

- OmniCloudMask support: install with `.[cloud]`
- Documentation tooling: install with `.[docs]`
- Elastix registration support: install with `.[elastix]`

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

With Py6S atmospheric correction extras (not required; included in base install):

```bash
pip install -e ".[py6s]"
```

With spectralmatch normalization extras:

```bash
pip install -e ".[spectralmatch]"
```

With elastix registration extras:

```bash
pip install -e ".[elastix]"
```

If you created your env from `environment.yml`, install elastix support separately:

```bash
pip install itk-elastix>=0.19.2
```

With docs extras:

```bash
pip install -e ".[docs]"
```

## Verify CLI Entry Points

```bash
vhr-worldview --help
vhr-fetch-modis-water-vapor --help
vhr-cloudmask-raster --help
vhr-pansharpen-orthos --help
vhr-align-image-pair --help
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
