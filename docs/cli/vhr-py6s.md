# `vhr-py6s`

Run **Py6S atmospheric correction only** on WorldView scenes, with optional post-Py6S orthorectification.

This command discovers scenes the same way as `vhr-worldview`, then writes a corrected MS raster per scene.

## Example

```bash
vhr-py6s \
  --input-dir /data/worldview_batch \
  --output-dir /data/out/py6s_only \
  --output-suffix _py6s \
  --py6s-atmosphere-profile user \
  --py6s-aerosol-profile maritime \
  --py6s-auto-atmos-source nasa_power

# Optional: produce projected output (for example EPSG:6635)
vhr-py6s \
  --input-dir /data/worldview_batch \
  --output-dir /data/out/py6s_only \
  --output-suffix _py6s \
  --epsg 6635 \
  --dem-file-path /data/dem.tif
```

## Key Inputs

- `--input-dir` (repeatable): input root(s) containing WorldView scene structure
- `--filter-basename`: optional specific scene basename(s)
- `--output-dir`: output folder for Py6S rasters and metadata reports
- `--output-suffix`: output filename suffix
- `--epsg`: optional projected output EPSG; enables post-Py6S orthorectification
- `--dem-file-path`: DEM required when `--epsg` is set
- `--keep-intermediate-py6s`: keep raw Py6S output and write an additional `*_ortho.tif`

## Py6S Parameters

- `--py6s-atmosphere-profile`
- `--py6s-aerosol-profile`
- `--py6s-aot550`
- `--py6s-visibility` (overrides AOT when set)
- `--py6s-water-vapor`
- `--py6s-ozone`
- `--py6s-use-imd-radiance-calibration`
- `--py6s-use-worldview-gain-offset-adjustment`
- `--py6s-auto-atmos-source {none,nasa_power}`

## Outputs

For each scene:

- `<basename><suffix>.tif`: Py6S-corrected raster
  - without `--epsg`: output is RPC/GCP-referenced (not a projected orthorectified grid)
  - with `--epsg`: output is orthorectified to the requested EPSG
- `<basename><suffix>_metadata.json`: parameters and effective values used
