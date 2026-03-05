# `vhr-cloudmask-raster`

Standalone cloud-mask runner for an existing raster using OmniCloudMask.

## Usage

```bash
vhr-cloudmask-raster \
  --input-raster /data/out/scene_final.tif \
  --output-raster /data/out/scene_final_cloudmasked.tif \
  --output-mask /data/out/scene_final_cloudmask.tif \
  --buffer-pixels 3
```

## Key Options

- `--cloud-classes`: comma-separated cloud/cloud-shadow classes
- `--red-band-index`, `--green-band-index`, `--nir-band-index`
- `--omnicloud-kwargs-json`: pass model kwargs directly
