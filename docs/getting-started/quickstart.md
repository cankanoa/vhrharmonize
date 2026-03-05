# Quickstart

## 1. Start from the default config

Create a local config from the tracked template, then adjust paths for your machine.

```bash
cp configs/worldview.example.yml configs/worldview.yml
```

## 2. Run full-scene preprocessing

`vhr-worldview` runs the main end-to-end workflow for each discovered scene:

- scene/file discovery and metadata loading
- atmospheric correction (`py6s`, `flaash`, or `none`; default `py6s`)
- RPC orthorectification (MS and PAN)
- pansharpening
- optional cloud masking (if enabled in args/config)
- final output write to your configured output directory

```bash
vhr-worldview --config-yaml configs/worldview.yml
```

## Optional standalone tools

### Fetch MODIS water vapor (optional)

Fetches per-scene atmospheric water vapor from Google Earth Engine MODIS collections, using scene footprint + acquisition time, and writes JSON/CSV reports.

```bash
vhr-fetch-modis-water-vapor \
  --input-dir /data/worldview_batch \
  --output-json outputs/modis_water_vapor_results.json \
  --output-csv outputs/modis_water_vapor_results.csv
```

### Cloud mask an existing raster (optional)

```bash
vhr-cloudmask-raster \
  --input-raster /data/out/scene_final.tif \
  --buffer-pixels 3
```

### Pansharpen existing orthorectified rasters (optional)

```bash
vhr-pansharpen-orthos \
  --mul-ortho /data/mul_ortho.tif \
  --pan-ortho /data/pan_ortho.tif \
  --output /data/out/pansharp.tif
```

### Align one image to another with elastix (optional)

```bash
vhr-align-image-pair \
  --moving-image /data/image_a_cloudmasked.tif \
  --fixed-image /data/image_b_lidar.tif \
  --output-image /data/out/image_a_aligned_to_b.tif
```

### Run Py6S atmospheric correction only (optional)

```bash
vhr-py6s \
  --input-dir /data/worldview_batch \
  --output-dir /data/out/py6s_only \
  --output-suffix _py6s \
  --py6s-atmosphere-profile user \
  --py6s-auto-atmos-source nasa_power
```
