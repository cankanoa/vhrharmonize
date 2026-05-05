# Quickstart

## WorldView-3
WorldView-3 is the only currently streamlined end-to-end sensor workflow. For easy setup,use the example config as your starting point at [configs/worldview.example.yml](configs/worldview.example.yml)

### Full workflow
Assuming you already have the library installed.
1. Obtain the example config:
   2. If you have installed from Pypi as a Python library, download this file: [configs/worldview.example.yml](configs/worldview.example.yml). 
   3. If you have cloned the repository, navigate to the configs/worldview.example.yml file.
2. Run the config from the command line with:

```bash
vhr-worldview --config-yaml configs/worldview.yml
```
Or specify the parameters directly in the command line:
```bash
# Example: override a few settings at runtime
vhr-worldview \
  --config-yaml configs/worldview.yml \
  --input-file-glob "/data/worldview/**/*.TIF" \
  --output-dir ../../processed \
  --run-alignment \
  --alignment-fixed-image /data/reference.tif
```

```bash
# Example: run with a local DEM and limit cloudier scenes
vhr-worldview \
  --dem-file-path /data/dem.tif \
  --max-cloud-cover-to-process 50 \
  --concurrent-processing 4
```

## Individual steps

```bash
# Fetch MODIS water vapor reports
vhr-fetch-modis-water-vapor \
  --input-dir /data/worldview_batch \
  --ee-project your-ee-project \
  --output-json outputs/modis.json \
  --output-csv outputs/modis.csv
```

```bash
# Run FLAASH from a prepared parameter file
vhr-flaash \
  --params-json-file outputs/flaash_params.json \
  --output-params-path outputs/flaash_params_used.json \
  --envi-engine-path /path/to/taskengine.exe
```

```bash
# Mask an existing raster
vhr-cloudmask-raster \
  --input-raster /data/pansharpened.tif \
  --output-raster /data/pansharpened_cloudmasked.tif \
  --output-mask /data/pansharpened_cloudmask.tif
```

```bash
# Pansharpen orthorectified rasters
vhr-pansharpen-orthos \
  --mul-ortho /data/mul_ortho.tif \
  --pan-ortho /data/pan_ortho.tif \
  --output /data/pansharpened.tif
```

```bash
# Align one raster to another
vhr-align-image-pair \
  --moving-image /data/moving.tif \
  --fixed-image /data/fixed.tif \
  --output-image /data/aligned.tif
```

```bash
# Orthorectify directly
vhr-orthorectification orthorectify \
  --input-image-path /data/input.tif \
  --output-image-path /data/output_ortho.tif \
  --dem-image-path /data/dem.tif \
  --output-epsg 6635
```

```bash
# Run Py6S-only processing
vhr-py6s \
  --input-dir /data/worldview_batch \
  --output-dir /data/py6s_only \
  --output-suffix _py6s
```

```bash
# Run radiometric normalization directly
vhr-radiometric-normalization \
  --input-image /data/image_a.tif \
  --input-image /data/image_b.tif \
  --output-image /data/normalized.tif
```
