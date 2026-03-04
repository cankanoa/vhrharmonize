# `vhr-fetch-modis-water-vapor`

Fetches MODIS water-vapor estimates per scene and computes recommended `WATER_VAPOR_PRESET` values.

## Typical Usage

```bash
vhr-fetch-modis-water-vapor \
  --input-dir /data/worldview_batch \
  --output-json outputs/modis_water_vapor_results.json \
  --output-csv outputs/modis_water_vapor_results.csv
```

## Earth Engine Setup

- Set project id via `--ee-project` or an env file (`--env-file`)
- Optional first-time auth with `--authenticate`

## Output

- JSON report (`--output-json`)
- CSV report (`--output-csv`)
- Optional write-back into `Root_WV.txt`:
  - `--write-root-overrides`
  - `--write-atmosphere-overrides`

## Notes

- Defaults target daily MODIS Terra/Aqua atmosphere collections.
- AOD and water-vapor scale factors are configurable.
