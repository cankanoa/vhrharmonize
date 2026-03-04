# Atmospheric Correction with Py6S

This page describes exactly how `vhr-worldview` uses Py6S, which config values are actually consumed, and how to choose/supply them.

## Processing Context

In the WorldView workflow, atmospheric correction runs before orthorectification/pansharpening:

1. Parse scene metadata (IMD + optional root overrides)
2. Run atmospheric correction (`py6s`, `flaash`, or `none`)
3. Orthorectify MS and PAN
4. Pansharpen
5. Optional cloud masking

When `atmospheric_method: py6s`, the Py6S implementation writes **surface reflectance** by default as scaled integer values:

- `py6s_output_scale_factor: 10000.0`
- `py6s_output_dtype: int16`

Interpretation:

- `3560` means reflectance `0.356`
- `10000` means reflectance `1.0`

## Which Inputs Are Auto-Derived

The following are auto-derived per scene and do not usually need manual edits:

- solar/view geometry from IMD
- acquisition day/month from scene basename
- terrain elevation percentile from DEM + footprint
- WorldView per-band `absCalFactor/effectiveBandwidth` from IMD when IMD radiance calibration is enabled

## Py6S Fields in `configs/worldview.yml`

These fields are used when `atmospheric_method: py6s`:

- `py6s_atmosphere_profile`
- `py6s_aerosol_profile`
- `py6s_aot550`
- `py6s_visibility` (if set, this overrides `py6s_aot550`)
- `py6s_water_vapor`
- `py6s_ozone`
- `py6s_output_scale_factor`
- `py6s_output_dtype`
- `py6s_use_imd_radiance_calibration`
- `py6s_use_worldview_gain_offset_adjustment`
- `py6s_auto_atmos_source`
- `py6s_auto_atmos_grid_size`
- `py6s_auto_atmos_search_days`
- `py6s_auto_atmos_timeout_s`
- `py6s_executable` (optional explicit path)

## Parameter Precedence Rules

1. CLI flags override YAML values.
2. If `py6s_visibility` is set, Py6S uses visibility mode and ignores `py6s_aot550`.
3. IMD radiance calibration:
   - if `py6s_use_imd_radiance_calibration: true`, workflow uses IMD `absCalFactor/effectiveBandwidth`
   - if `py6s_use_worldview_gain_offset_adjustment: true`, an additional WV2/WV3 gain/offset table is applied
4. Atmosphere profile behavior:
   - if `py6s_atmosphere_profile: user`, user `py6s_water_vapor` and `py6s_ozone` are used directly
   - if profile is not `user` (for example `midlatitude_summer`), the named profile is used
5. Auto-source behavior:
   - if `py6s_auto_atmos_source: nasa_power` and profile is `user`,
     scene-date + bbox values auto-update `aot550`, `water_vapor`, `ozone`

## Recommended Defaults (Current Project Practice)

```yaml
atmospheric_method: py6s
py6s_atmosphere_profile: user
py6s_aerosol_profile: maritime
py6s_aot550: 0.2
# py6s_visibility: null
py6s_water_vapor: 2.5
py6s_ozone: 0.3
py6s_output_scale_factor: 10000.0
py6s_output_dtype: int16
py6s_use_imd_radiance_calibration: true
py6s_use_worldview_gain_offset_adjustment: true
py6s_auto_atmos_source: nasa_power
py6s_auto_atmos_grid_size: 3
py6s_auto_atmos_search_days: 1
```

Use `maritime` for island/coastal scenes by default. For inland scenes, test `continental`.

## Aerosol Profile Options

Supported aerosol profiles in this implementation:

- `continental`
- `maritime`
- `urban`
- `desert`
- `biomass_burning`
- `stratospheric`

## How to Source Atmospheric Inputs

### Water Vapor

- Can be estimated scene-wise from MODIS water vapor products (manual workflow or your standalone fetch script).
- Set result via `py6s_water_vapor`.

### AOT550

- Can be sourced from aerosol products (including MODIS aerosol workflows).
- Set via `py6s_aot550`, or provide `py6s_visibility` instead.

### Ozone

- Typically sourced from ozone products (not usually MODIS water-vapor workflow).
- Set via `py6s_ozone`.

### Automatic NASA API Mode

You can automatically source all three fields from NASA POWER daily API:

- `py6s_auto_atmos_source: nasa_power`
- requires `py6s_atmosphere_profile: user`
- computes bbox sample-grid mean from:
  - `AOD_55` -> `py6s_aot550`
  - `PW` (cm) -> `py6s_water_vapor` (g/cm^2)
  - `TO3` (Dobson Units) -> `py6s_ozone` (cm-atm via `DU / 1000`)

CLI example:

```bash
vhr-worldview \
  --config-yaml configs/worldview.yml \
  --py6s-atmosphere-profile user \
  --py6s-auto-atmos-source nasa_power \
  --py6s-auto-atmos-grid-size 3 \
  --py6s-auto-atmos-search-days 1
```

## Should You Change Values?

Start simple:

1. Keep defaults and confirm radiometric behavior.
2. Switch to `py6s_atmosphere_profile: user` if you want explicit control over water vapor/ozone.
3. Tune only one variable at a time (`aot550` or aerosol profile first).
4. Compare outputs with stable visualization settings (same band combo/stretch).

## Example Run

```bash
vhr-worldview \
  --config-yaml configs/worldview.yml \
  --atmospheric-method py6s \
  --py6s-atmosphere-profile user \
  --py6s-aot550 0.2 \
  --py6s-water-vapor 2.5 \
  --py6s-ozone 0.3
```

## Calibration Reference Links

- WorldView-3 DN->Radiance gain/offset reference used in this project:
  - https://docs.mcube.terradue.com/pre-processing/opt/worldview-3/
- Additional WV radiometric references:
  - https://nctec.co/wp-content/uploads/2018/01/2018_03_28_DigitalGlobe_Radiometric_Use_of_WorldView-2_Imagery_RevA.pdf
  - https://dg-cms-uploads-production.s3.amazonaws.com/uploads/document/file/95/Radiometric_Use_of_WorldView-3_Imagery.pdf
