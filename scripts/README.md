# WorldView CLI Workflow

This folder provides a command-line wrapper for the workflow in `docs/examples/worldview_b1.py`.

## Script

- `scripts/worldview_cli.py`
- `scripts/fetch_modis_water_vapor_gee.py`
- `scripts/run_flaash_matrix.py`
- `scripts/run_tile_elastix_registration.py`
- `scripts/run_tiled_image_lidar_elastix.py`

## What It Does

For each scene discovered from `Root*.txt` files, it performs:

1. Footprint conversion (`.shp` -> `.gpkg`)
2. Optional LiDAR tile-stat raster generation from EPT (mean of configured dimension by tile bounds)
3. FLAASH atmospheric correction (unless skipped)
4. Default-RPC orthorectification for multispectral and panchromatic images (no GCP refinement)
5. Pansharpening
6. Tile clipping from a tile `gpkg`
7. Optional cloud masking per tile (so cloud mask runs on small tile rasters, not full scene)
8. Writing final tile outputs only (intermediates stay in scratch)

Orthorectification is default RPC only in this wrapper.

## LiDAR Mean Raster From EPT

If `--lidar-ept-path` is set, the CLI writes one raster per tile using `pyforestscan.handlers.read_lidar` with EPT bounds set to each buffered tile extent, so it does not load the full EPT.

- Statistic: mean of `--lidar-dimension` (default `Intensity`)
- Grid size: `--lidar-resolution` (required when LiDAR step enabled)
- Extent: tile polygon bounds after `--tile-buffer-m`

Example:

```bash
python scripts/worldview_cli.py \
  --config-yaml scripts/worldview_config.yml \
  --lidar-ept-path /data/lidar/ept.json \
  --lidar-dimension Intensity \
  --lidar-resolution 1.0 \
  --lidar-output-dir /data/out/lidar_mean
```

## Standalone Elastix Tile Registration

Use this step after tiles already exist to register each moving imagery tile to its matching fixed intensity tile.

- Fixed: intensity tiles (`<tile_id>_<dimension><intensity_suffix>.tif`)
- Moving: imagery tiles (`<tile_id>_<scene><moving_suffix>.tif`)
- Final aligned output: written next to each moving tile as `<moving_basename>_aligned.tif` (configurable suffix)
- Work output: per-pair elastix folders (transforms/intermediates) under `--output-dir`

Required environment variables (used by elastix wrapper):

```bash
export ELASTIX_EXE="/mnt/c/.../elastix.exe"
export TRANSFORMIX_EXE="/mnt/c/.../transformix.exe"
```

Run with YAML:

```bash
python scripts/run_tile_elastix_registration.py \
  --config-yaml scripts/elastix_tile_config.yml
```

Or run directly:

```bash
python scripts/run_tile_elastix_registration.py \
  --moving-dir /mnt/d/test_data/out/tiles \
  --intensity-dir /mnt/d/test_data/out/lidar_mean \
  --output-dir /mnt/d/test_data/out/elastix_reg_work \
  --moving-suffix _final \
  --intensity-dimension Intensity \
  --intensity-suffix _mean \
  --moving-band-index 1 \
  --aligned-suffix _aligned
```

Notes:

- Moving band defaults to `1` (WV3 coastal blue).
- Existing registration outputs are skipped unless `--overwrite` is set.
- The script derives `--epsg`, `--ullr`, and `--res` from each fixed intensity tile and passes them to the elastix wrapper.

## Image + EPT + Elastix (No FLAASH)

Use this for a direct workflow on a large input image and EPT source when you do not need FLAASH/orthorectify/pansharpen:

1. Clip image to tiles
2. Build tile-wise LiDAR mean rasters from EPT bounds
3. Align each image tile to its intensity tile using elastix

Script:

- `scripts/run_tiled_image_lidar_elastix.py`

Config template for your current dataset:

- `scripts/puuanahulu_tiling_elastix.yml`

Run:

```bash
export ELASTIX_EXE="/mnt/c/.../elastix.exe"
export TRANSFORMIX_EXE="/mnt/c/.../transformix.exe"

python scripts/run_tiled_image_lidar_elastix.py \
  --config-yaml scripts/puuanahulu_tiling_elastix.yml
```

Outputs:

- image tiles: `image_tile_output_dir/<tile_id>_<image_basename>_final.tif`
- intensity tiles: `lidar_tile_output_dir/<tile_id>_Intensity_mean.tif`
- aligned tiles (next to image tiles): `..._final_aligned.tif` (all bands)
- elastix intermediates: `elastix_work_dir/...`

Robustness notes:

- The script skips tiles with too few valid pixels in moving/fixed rasters (`min_valid_pixels`, default `100`).
- Elastix failures are logged per tile and processing continues for remaining tiles.
- Registration is estimated from `moving_band_index` (default band 1), then the rigid transform is applied to all bands.

## Requirements

- Python environment with package dependencies installed (see project `pyproject.toml`)
- ENVI Task Engine path for FLAASH runs
- DEM in WGS84 ellipsoidal height that covers all scenes
- For built-in cloud masking: `pip install -e ".[cloud]"`

## Basic Usage

```bash
python scripts/worldview_cli.py \
  --input-dir /data/worldview_batch \
  --dem-file-path /data/dem/dem_wgs84_ellipsoidal.tif \
  --envi-engine-path "C:/Program Files/NV5/ENVI56/taskengine.exe" \
  --epsg 4326
```

## MODIS Water Vapor Helper (GEE)

Use this script to fetch per-scene MODIS water vapor from Google Earth Engine and compute a recommended `WATER_VAPOR_PRESET` value:

Set your GEE project id in `scripts/.env`:

```dotenv
GEE_PROJECT=your-ee-project-id
```

```bash
python scripts/fetch_modis_water_vapor_gee.py \
  --input-dir /data/worldview_batch \
  --output-json scripts/modis_water_vapor_results.json \
  --output-csv scripts/modis_water_vapor_results.csv
```

Defaults use daily MODIS atmosphere collections:
- Terra: `MODIS/061/MOD08_D3`
- Aqua: `MODIS/061/MYD08_D3`
- Band: `Atmospheric_Water_Vapor_Mean`
- Time window: `24` hours

Write computed values directly into each `Root_WV.txt` (under `ParamsOverridesPerPhoto`):

```bash
python scripts/fetch_modis_water_vapor_gee.py \
  --input-dir /data/worldview_batch \
  --write-root-overrides
```

Write both water-vapor and atmosphere/aerosol overrides:

```bash
python scripts/fetch_modis_water_vapor_gee.py \
  --input-dir /data/worldview_batch \
  --write-root-overrides \
  --write-atmosphere-overrides
```

CLI still overrides `.env`:

```bash
python scripts/fetch_modis_water_vapor_gee.py \
  --input-dir /data/worldview_batch \
  --ee-project my-other-project
```

## Minimal FLAASH Matrix Test

Initial defaults used by the original workflow were:
- `MODTRAN_ATM = Mid-Latitude Summer`
- `MODTRAN_AER = Maritime`
- `USE_AEROSOL = Disabled`
- `DEFAULT_VISIBILITY = None`

Run a minimal aerosol/atmosphere test matrix on one scene:

```bash
python scripts/run_flaash_matrix.py \
  --input-dir /mnt/d/20171019_35cm_WV03_BAB_050311750010 \
  --photo-basename 17OCT19211717-M1BS-050311750010_01_P001 \
  --dem-file-path /mnt/d/dem/Hawaii_SRTM_GL1Ellip.tif \
  --envi-engine-path "/mnt/c/Program Files/NV5/ENVI60/IDL90/bin/bin.x86_64/taskengine.exe" \
  --scratch-dir /mnt/d/test_data/flaash_matrix \
  --dem-ground-percentile 100 \
  --output-csv scripts/flaash_matrix_results.csv
```

The CSV reports sampled negative percentages per band for each scenario.

If needed, run interactive GEE auth once:

```bash
python scripts/fetch_modis_water_vapor_gee.py \
  --input-dir /data/worldview_batch \
  --authenticate
```

Use a YAML config file:

```yaml
# scripts/worldview_config.yml
input_dir:
  - /data/worldview_batch_a
  - /data/worldview_batch_b
dem_file_path: /data/dem/dem_wgs84_ellipsoidal.tif
envi_engine_path: /opt/envi/taskengine.exe
epsg: 4326
nodata_value: -9999
cloud_mask_method: omnicloudmask
cloud_buffer_pixels: 3
tile_gpkg: /data/tiles/tiles.gpkg
tile_layer: tiles
tile_id_field: tile_id
tile_buffer_m: 100
tile_output_dir: /data/out/tiles
tile_output_suffix: _final
scratch_dir: /mnt/d/vhr_scratch
lidar_ept_path:
lidar_dimension: Intensity
lidar_resolution: 1.0
lidar_output_dir: /data/out/lidar_mean
lidar_output_suffix: _mean
lidar_overwrite: false
```

```bash
python scripts/worldview_cli.py --config-yaml scripts/worldview_config.yml
```

CLI args still override YAML values:

```bash
python scripts/worldview_cli.py \
  --config-yaml scripts/worldview_config.yml \
  --epsg 6635
```

Tile mode (write only final per-tile outputs):

```bash
python scripts/worldview_cli.py \
  --input-dir /data/worldview_batch \
  --dem-file-path /data/dem/dem_wgs84_ellipsoidal.tif \
  --envi-engine-path "/opt/envi/taskengine.exe" \
  --tile-gpkg /data/tiles/tiles.gpkg \
  --tile-layer tiles \
  --tile-id-field tile_id \
  --tile-buffer-m 100 \
  --tile-output-dir /data/out/tiles \
  --tile-output-suffix _final
```

Process multiple roots:

```bash
python scripts/worldview_cli.py \
  --input-dir /data/worldview_batch_a \
  --input-dir /data/worldview_batch_b \
  --dem-file-path /data/dem/dem_wgs84_ellipsoidal.tif \
  --envi-engine-path "/opt/envi/taskengine.exe"
```

Filter to specific basenames:

```bash
python scripts/worldview_cli.py \
  --input-dir /data/worldview_batch \
  --dem-file-path /data/dem/dem_wgs84_ellipsoidal.tif \
  --envi-engine-path "/opt/envi/taskengine.exe" \
  --filter-basename 17DEC08211758-M1BS-016445319010_01_P003 \
  --filter-basename 17DEC08211758-P1BS-016445319010_01_P003
```

Skip FLAASH and use an existing multispectral `.dat`:

```bash
python scripts/worldview_cli.py \
  --input-dir /data/worldview_batch \
  --dem-file-path /data/dem/dem_wgs84_ellipsoidal.tif \
  --skip-flaash \
  --existing-flaash-input /data/existing/Mul_FLAASH/example.dat
```

Run optional cloud masking after pansharpening:

```bash
python scripts/worldview_cli.py \
  --input-dir /data/worldview_batch \
  --dem-file-path /data/dem/dem_wgs84_ellipsoidal.tif \
  --envi-engine-path "/opt/envi/taskengine.exe" \
  --cloud-mask-method omnicloudmask \
  --cloud-buffer-pixels 3
```

Optional custom cloud classes and model kwargs:

```bash
python scripts/worldview_cli.py \
  --input-dir /data/worldview_batch \
  --dem-file-path /data/dem/dem_wgs84_ellipsoidal.tif \
  --envi-engine-path "/opt/envi/taskengine.exe" \
  --cloud-mask-method omnicloudmask \
  --cloud-mask-classes 1,2,3 \
  --cloud-mask-omnicloud-kwargs-json '{"batch_size": 8}'
```

## CLI Arguments

- `--input-dir` (required, repeatable): root folder(s) to scan for `Root*.txt`
- `--dem-file-path` (required): DEM path
- `--config-yaml` (optional): YAML config file with argument keys (supports `snake_case` or `kebab-case`)
- `--envi-engine-path` (required unless `--skip-flaash`): ENVI task engine executable path
- `--epsg` (default: `4326`): output EPSG for orthorectified imagery
- `--nodata-value` (default: `-9999`): output NoData value
- `--dtype` (default: `int16`): GDAL data type for orthorectification output
- `--flaash-dem-ground-percentile` (default: `50`): DEM percentile used for `GROUND_ELEVATION` (median by default)
- `--flaash-modtran-atm` (default: `Mid-Latitude Summer`): base atmospheric model (overridable in `Root_WV.txt`)
- `--flaash-modtran-aer` (default: `Maritime`): base aerosol model (overridable in `Root_WV.txt`)
- `--flaash-use-aerosol` (default: `Disabled`): base aerosol retrieval setting
- `--flaash-default-visibility` (optional): base visibility in km
- `--footprint-epsg` (default: `4326`): CRS assigned to Maxar footprints before writing GeoPackage
- `--filter-basename` (optional, repeatable or comma-separated): process only matching basenames
- `--tile-gpkg` (required): tile polygons GeoPackage for clipping outputs
- `--tile-layer` (optional): layer name in the tile GeoPackage
- `--tile-id-field` (default: `tile_id`): attribute used in output filenames
- `--tile-buffer-m` (default: `100`): buffer distance applied to tile polygons before overlap/clipping (map units; meters for EPSG:6635)
- `--tile-output-dir` (default: `Tile_Output`): directory for final tile outputs
- `--tile-output-suffix` (default: `_final`): suffix for final output filenames
- `--lidar-ept-path` (optional): EPT metadata path (`ept.json`) for LiDAR mean tile rasters
- `--lidar-srs` (optional): LiDAR SRS override (auto-read from EPT when omitted)
- `--lidar-dimension` (default: `Intensity`): LiDAR dimension aggregated as mean
- `--lidar-resolution` (required with `--lidar-ept-path`): output raster resolution in LiDAR CRS units
- `--lidar-output-dir` (default: `<tile-output-dir>/lidar_mean`): output directory for LiDAR mean rasters
- `--lidar-output-suffix` (default: `_mean`): suffix for LiDAR mean output filenames
- `--lidar-nodata` (default: `-9999`): NoData value for LiDAR mean rasters
- `--lidar-overwrite` (flag): overwrite existing LiDAR mean tiles (default behavior is skip existing)
- `--scratch-dir` (default: `/tmp`): temporary workspace for intermediates
- `--skip-flaash` (flag): skip FLAASH stage
- `--existing-flaash-input` (required when `--skip-flaash`): multispectral `.dat` used as FLAASH output substitute
- `--cloud-mask-method` (optional): built-in cloud masking mode (`omnicloudmask`)
- `--cloud-mask-red-band-index` (default: `5`): 1-based red band index for omnicloudmask
- `--cloud-mask-green-band-index` (default: `3`): 1-based green band index for omnicloudmask
- `--cloud-mask-nir-band-index` (default: `7`): 1-based NIR band index for omnicloudmask
- `--cloud-mask-classes` (default: `1,2,3`): classes treated as cloud/cloud-shadow
- `--cloud-buffer-pixels` (default: `0`): dilation buffer around cloud regions
- `--cloud-mask-omnicloud-kwargs-json` (optional): JSON kwargs forwarded to `predict_from_array`
- `--cloud-mask-output-suffix` (default: `_cloudmasked`): output filename suffix for masked image
- `--cloud-mask-mask-suffix` (default: `_cloudmask`): output filename suffix for binary cloud mask
- `--cloud-mask-command` (optional): custom shell command hook for per-tile masking (kept for non-omnicloud tools)

## Notes

- The script now runs through a temporary scratch workspace and writes final outputs only.
- It expects Maxar-style filenames and sidecar files (`.IMD`, `.TIF`, `.shp`) under each scene root.
- Scenes are skipped early if they do not intersect any tile polygon.
- LiDAR EPT reads are bounded per buffered tile (`read_lidar(..., bounds=...)`) to avoid loading the full EPT.
- Existing LiDAR tile outputs are skipped unless `--lidar-overwrite` is set.
- When running FLAASH through Windows ENVI from WSL, set `--scratch-dir` on a Windows-mounted path (for example `/mnt/d/...`), not `/tmp`.
- Built-in omnicloudmask writes two outputs when enabled:
  - binary mask (`*_cloudmask.tif` by default)
  - masked image (`*_cloudmasked.tif` by default)
- Cloud-mask command placeholders:
  - `{input}`: per-tile clipped raster path
  - `{output}`: cloud-masked output path
  - `{scene_root}`: scene folder containing `Root*.txt`
  - `{image_basename}`: multispectral photo basename
