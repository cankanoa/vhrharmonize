# `vhr-worldview`

Runs full-scene WorldView preprocessing:

1. Recursive scene discovery from `--input-file-glob`
2. WorldView MUL/PAN pairing plus IMD parsing
3. Optional atmosphere fetch
4. Atmospheric correction (`flaash`, `py6s`, or `none`)
5. Orthorectification (default RPC)
6. Pansharpening
7. Optional cloud masking
8. Optional alignment to a fixed/reference raster
9. Optional radiometric normalization
10. Write scene metadata beside the last enabled step output

## Typical Usage

Create local config first:

```bash
cp configs/worldview.example.yml configs/worldview.yml
```

```bash
vhr-worldview --config-yaml configs/worldview.yml
```

Run with Py6S atmospheric correction:

```bash
vhr-worldview \
  --config-yaml configs/worldview.yml \
  --atmospheric-method py6s \
  --py6s-aot550 0.2 \
  --py6s-water-vapor 2.5 \
  --py6s-ozone 0.3
```

If 6S is not on `PATH`, pass it explicitly:

```bash
vhr-worldview \
  --config-yaml configs/worldview.yml \
  --atmospheric-method py6s \
  --py6s-executable /path/to/sixsV1.1
```

## Required Inputs

- `--input-file-glob` (repeatable) unless set in config
- `--dem-file-path` unless set in config
- `--envi-engine-path` only when `--atmospheric-method=flaash` and not using `--skip-flaash`

## Common Options

- `--output-dir`: base output folder used when a step has `save_<step>=output`
  - relative paths resolve from the MUL image folder
  - if unset, the default output base is `../Processed` from the MUL image folder
- `--skip-existing`: skip scenes when the final output TIFF already exists
- `--last-run-step`: resume the raster chain from an existing step output
- `--atmospheric-method`: choose `flaash`, `py6s`, or `none`
- `--run-fetch-atmosphere` / `--no-run-fetch-atmosphere`
- `--save-fetch-atmosphere`: `temp`, `output`, or a custom folder
- `--save-atmospheric-correction`: `temp`, `output`, or a custom folder
- `--save-orthorectification`: `temp`, `output`, or a custom folder
- `--save-pansharpen`: `temp`, `output`, or a custom folder
- `--save-cloud-mask`: `temp`, `output`, or a custom folder
- `--save-alignment`: `temp`, `output`, or a custom folder
- `--save-radiometric-normalization`: `temp`, `output`, or a custom folder
- `--skip-flaash` + `--existing-flaash-input`: resume from existing FLAASH output
- `--py6s-*`: Py6S atmosphere/aerosol/output controls
- `--py6s-visibility`: optional visibility mode (km), used instead of `--py6s-aot550`
- `--py6s-executable`: explicit path to 6S binary if not in `PATH`
- `--py6s-use-imd-radiance-calibration` / `--no-py6s-use-imd-radiance-calibration`: toggle WorldView IMD-based DN->radiance before Py6S (default: enabled)
- `--py6s-use-worldview-gain-offset-adjustment` / `--no-py6s-use-worldview-gain-offset-adjustment`: toggle built-in WV02/WV03 gain/offset table (default: enabled)
- `--py6s-auto-atmos-source`: `none` or `nasa_power` auto-fetch for `aot550/water_vapor/ozone` when profile is `user`
- `--py6s-auto-atmos-grid-size`: sample-grid size over scene bbox for NASA auto mode
- `--py6s-auto-atmos-search-days`: +/- day search window for NASA auto mode
- `--existing-mul-ortho-input` + `--existing-pan-ortho-input`: resume at pansharpen step
- `--run-atmospheric-correction` / `--no-run-atmospheric-correction`: enable/skip atmospheric correction stage
- `--run-pansharpen` / `--no-run-pansharpen`: enable/skip pansharpen stage
- `--run-cloud-mask` / `--no-run-cloud-mask`: enable or skip cloud masking without editing YAML
- `--cloud-mask-method omnicloudmask`: run built-in cloud masking
- `--cloud-mask-inference-resolution-m`: resolution used for OmniCloudMask inference input (default `10.0`)
- `--cloud-mask-command`: run external command template on full-scene image
- `--run-alignment` / `--no-run-alignment`: enable or skip final alignment stage
- `--run-radiometric-normalization` / `--no-run-radiometric-normalization`
- `--alignment-fixed-image`: required fixed/reference raster when alignment is enabled
- `--alignment-moving-band-index`: 0-based moving image band used for registration metric
- `--alignment-fixed-band-index`: 0-based fixed image band used for registration metric
- `--alignment-registration-mode`: `default` or `structural_wv3_lidar`
- `--alignment-clip-fixed-to-moving`: clip fixed domain to moving extent before alignment
- `--alignment-output-on-moving-grid`: write aligned output on moving-image grid/resolution (default: enabled)
- `--alignment-temp-dir`: temp directory for elastix/transformix artifacts
- `--alignment-keep-temp-dir`: retain alignment temp files for inspection/debugging

Alignment mode notes:
- `default`: raw-band elastix registration
- `structural_wv3_lidar`: common-grid structural registration path for optical-to-LiDAR alignment
  - despite the name, this is a generic structural edge-based mode, not a WorldView-only implementation

For large scenes, prefer `alignment_temp_dir` on a native Linux filesystem rather
than a mounted external path.

## Py6S Output Units

- By default, Py6S writes **surface reflectance scaled by 10,000** as `int16`.
- Example conversion:
  - pixel value `3560` means reflectance `0.356`
  - pixel value `10000` means reflectance `1.0`
- This is controlled by:
  - `--py6s-output-scale-factor` (default `10000`)
  - `--py6s-output-dtype` (default `int16`)

## WorldView IMD Calibration Notes

- By default, the workflow reads per-band `absCalFactor` and `effectiveBandwidth`
  from the scene IMD and computes radiance input for Py6S.
- For WV02/WV03, a built-in Maxar gain/offset table is applied by default:
  - `radiance = gain * (DN * absCalFactor / effectiveBandwidth) + offset`
- By default, workflow uses `radiance = DN * (absCalFactor / effectiveBandwidth)`.

## WorldView Geometry Notes (Py6S)

- Py6S `view_zenith` is populated from IMD `meanOffNadirViewAngle` (degrees from nadir).
- FLAASH-specific LOS conversion (`180 - meanOffNadirViewAngle`) is **not** used for Py6S.

## Calibration References

- WorldView-3 DN->Radiance gain/offset table used in this workflow:
  - https://docs.mcube.terradue.com/pre-processing/opt/worldview-3/
- Additional WorldView radiometric technical notes:
  - https://nctec.co/wp-content/uploads/2018/01/2018_03_28_DigitalGlobe_Radiometric_Use_of_WorldView-2_Imagery_RevA.pdf
  - https://dg-cms-uploads-production.s3.amazonaws.com/uploads/document/file/95/Radiometric_Use_of_WorldView-3_Imagery.pdf

## Skip Existing

`--skip-existing` works in two layers:

- per step, expected output filenames are planned first and only missing outputs are processed
- before running the scene at all, the workflow checks whether the final output TIFF already exists and skips the whole scene if it does

This means intermediate temp outputs can be safely removed without breaking scene-level skip behavior.

## Output Model

Each raster step uses the same storage rule:

- `save_<step>=temp` writes to `<temp_dir>/<step_name>`
- `save_<step>=output` writes to `<output_dir>/<step_name>`
- any other value is treated as a custom folder, and the step name is appended to it

If `temp_dir` is unset, `vhr-worldview` creates a real temporary directory with Python `tempfile`.

There is no extra final copy stage anymore. The scene's final raster is the output from the last enabled raster step, and the metadata JSON is written beside that raster.

Default save behavior:

- `save_cloud_mask=output`
- all other `save_*` settings default to `temp`

## Defaults

- Atmospheric method default: `py6s`
- If you switch to `flaash`, ENVI Task Engine settings/requirements apply.
- Alignment default: disabled (`--no-run-alignment`).
