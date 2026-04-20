# CLI Overview

The project publishes these CLI entry points via `pyproject.toml`:

- `vhr-worldview`: full-scene workflow orchestration
- `vhr-fetch-modis-water-vapor`: per-scene MODIS water-vapor lookup via Earth Engine
- `vhr-flaash`: direct shared FLAASH runner
- `vhr-cloudmask-raster`: standalone OmniCloudMask runner on a raster
- `vhr-pansharpen-orthos`: standalone pansharpening from orthorectified inputs
- `vhr-align-image-pair`: pairwise elastix alignment (tile-aware by default)
- `vhr-orthorectification`: direct orthorectification tools
- `vhr-radiometric-normalization`: direct spectralmatch wrapper
- `vhr-py6s`: Py6S-only atmospheric correction on discovered scenes
- `vhr-planet`: Planet pipeline scaffold (placeholder for upcoming implementation)

## General Pattern

- Use `--config-yaml` where supported for reproducibility.
- Store generated reports/artifacts under `outputs/`.
- Keep long-lived configuration under `configs/`.

## Help

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
