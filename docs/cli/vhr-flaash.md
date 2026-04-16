# `vhr-flaash`

Run the shared FLAASH implementation directly from the CLI.

Install extras with:

```bash
pip install -e ".[flaash]"
```

## Usage

```bash
vhr-flaash \
  --params-json-file outputs/flaash_params.json \
  --output-params-path outputs/flaash_params_used.json \
  --envi-engine-path /path/to/taskengine.exe
```

## Notes

- pass exactly one of `--params-json` or `--params-json-file`
- the parameter object must include `INPUT_RASTER.url` and `OUTPUT_RASTER_URI`
- this command is a thin wrapper around the shared FLAASH runner in `preprocess.atmospheric_correction`
