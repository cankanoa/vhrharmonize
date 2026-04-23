# Quickstart

## 1. Start from the default config

Create a local config from the tracked template, then adjust paths for your machine.

```bash
cp configs/worldview.example.yml configs/worldview.yml
```

If you keep the default `dem_file_path: online`, set an OpenTopography API key first:

```bash
export OPENTOPOGRAPHY_API_KEY=your_key_here
```

## 2. Run full-scene preprocessing

`vhr-worldview` runs the main end-to-end workflow for each discovered scene:

- recursive TIFF discovery, WorldView scene grouping, and metadata loading
- optional atmosphere fetch
- atmospheric correction (`py6s`, `flaash`, or `none`; default `py6s`)
- RPC orthorectification (MS and PAN)
- pansharpening
- optional cloud masking (if enabled in args/config)
- optional alignment
- optional radiometric normalization
- step outputs written according to each `save_*` setting
- final raster taken from the last enabled raster step

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

### Run FLAASH directly (optional)

```bash
vhr-flaash \
  --params-json-file outputs/flaash_params.json \
  --output-params-path outputs/flaash_params_used.json \
  --envi-engine-path /path/to/taskengine.exe
```

### Run radiometric normalization directly (optional)

```bash
vhr-radiometric-normalization \
  rrn \
  --input-image /data/image_a.tif \
  --input-image /data/image_b.tif \
  --output-image /data/out/normalized.tif
```
