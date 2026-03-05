# `vhr-worldview`

Runs full-scene WorldView preprocessing:

1. Scene discovery from input directory (uses `Root*.txt` overrides if present)
2. Atmospheric correction (`flaash`, `py6s`, or `none`)
3. Orthorectification (default RPC)
4. Pansharpening
5. Optional cloud masking (inferred from orthorectified MS and applied to current workflow output)
6. Write final full-scene output

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

- `--input-dir` (repeatable) unless set in config
- `--dem-file-path` unless set in config
- `--envi-engine-path` only when `--atmospheric-method=flaash` and not using `--skip-flaash`

## Common Options

- `--output-dir`: final scene output folder
- `--output-suffix`: filename suffix for scene outputs
- `--atmospheric-method`: choose `flaash`, `py6s`, or `none`
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
- `--cloud-mask-command`: run external command template on full-scene image

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

## Validation Rules

The command enforces:

- valid DEM percentile range
- scratch directory existence
- required path pairing for existing ortho inputs
- mutual exclusion for `--cloud-mask-method` and `--cloud-mask-command`
- `/mnt/<drive>/...` scratch only when FLAASH pathing requires Windows ENVI interop

## Defaults

- Atmospheric method default: `py6s`
- If you switch to `flaash`, ENVI Task Engine settings/requirements apply.
