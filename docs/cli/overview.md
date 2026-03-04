# CLI Overview

The project publishes four primary CLI entry points via `pyproject.toml`:

- `vhr-worldview`: full-scene workflow orchestration
- `vhr-fetch-modis-water-vapor`: per-scene MODIS water-vapor lookup via Earth Engine
- `vhr-cloudmask-raster`: standalone OmniCloudMask runner on a raster
- `vhr-pansharpen-orthos`: standalone pansharpening from orthorectified inputs
- `vhr-align-image-pair`: pairwise elastix alignment (tile-aware by default)
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
vhr-cloudmask-raster --help
vhr-pansharpen-orthos --help
vhr-align-image-pair --help
vhr-py6s --help
```
