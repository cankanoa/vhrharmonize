# WorldView Config (`configs/worldview.yml`)

This is the primary config template for `vhr-worldview`.

## Example

```yaml
input_dir:
  - /mnt/d/20171019_35cm_WV03_BAB_050311750010
dem_file_path: /mnt/d/dem/Hawaii_SRTM_GL1Ellip.tif
atmospheric_method: py6s
py6s_atmosphere_profile: midlatitude_summer
py6s_aerosol_profile: maritime
py6s_aot550: 0.2
py6s_water_vapor: 2.5
py6s_ozone: 0.3
py6s_output_scale_factor: 10000.0
py6s_output_dtype: int16
py6s_use_imd_radiance_calibration: true
py6s_use_worldview_gain_offset_adjustment: false
# py6s_executable: /usr/local/bin/sixsV1.1

# Required only if atmospheric_method: flaash
# envi_engine_path: "/mnt/c/Program Files/NV5/ENVI60/IDL90/bin/bin.x86_64/taskengine.exe"

epsg: 6635
nodata_value: -9999
scratch_dir: /mnt/d/test_data/scratch

output_dir: /mnt/d/test_data/out/scenes
output_suffix: _final

cloud_mask_method: omnicloudmask
cloud_buffer_pixels: 10
cloud_mask_omnicloud_kwargs_json:
  inference_device: cuda
  mosaic_device: cpu
  batch_size: 1
  patch_size: 768
  patch_overlap: 192
```

## Field Notes

- `input_dir` can be a list; command-line flags still override config values.
- `atmospheric_method` defaults to `py6s`.
- `envi_engine_path` is only required when `atmospheric_method: flaash`.
- Py6S output defaults to scaled reflectance (`int16` with factor `10000`).
  - `reflectance = pixel_value / py6s_output_scale_factor`
- WorldView IMD calibration is enabled by default (`py6s_use_imd_radiance_calibration: true`).
  - Uses IMD per-band `absCalFactor` and `effectiveBandwidth` for DN->radiance input to Py6S.
  - Optional WV02/WV03 gain/offset table is controlled by `py6s_use_worldview_gain_offset_adjustment`.
  - Default is `false` because these table adjustments can over-darken some scenes.
  - Set `py6s_use_imd_radiance_calibration: false` to disable.
- `scratch_dir` is used for scene intermediates.
- `output_dir` controls final scene output location.
- `cloud_mask_method` currently supports `omnicloudmask`.
- `cloud_mask_omnicloud_kwargs_json` accepts either:
  - a YAML mapping/dictionary (recommended in YAML config), or
  - a JSON string (same format used by CLI flag).

## Run

```bash
vhr-worldview --config-yaml configs/worldview.yml
```

## Calibration References

- WorldView-3 DN->Radiance gain/offset table used in this workflow:
  - https://docs.mcube.terradue.com/pre-processing/opt/worldview-3/
- Additional WorldView radiometric technical notes:
  - https://nctec.co/wp-content/uploads/2018/01/2018_03_28_DigitalGlobe_Radiometric_Use_of_WorldView-2_Imagery_RevA.pdf
  - https://dg-cms-uploads-production.s3.amazonaws.com/uploads/document/file/95/Radiometric_Use_of_WorldView-3_Imagery.pdf
