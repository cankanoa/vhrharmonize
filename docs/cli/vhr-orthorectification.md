# `vhr-orthorectification`

Run shared orthorectification utilities directly from the CLI.

Install extras with:

```bash
pip install -e ".[orthorectification]"
```

## Subcommands

- `orthorectify`
- `qgis-gcps-to-csv`
- `qgis-gcps-to-geojson`

## Example

```bash
vhr-orthorectification orthorectify \
  --input-image /data/input.tif \
  --output-image /data/out/output_ortho.tif \
  --dem-file-path /data/dem.tif \
  --output-epsg 6635
```
