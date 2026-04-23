# WorldView Config (`configs/worldview.example.yml`)

`configs/worldview.example.yml` is the tracked template for `vhr-worldview`.

Create your local config first:

```bash
cp configs/worldview.example.yml configs/worldview.yml
```

## Example

```yaml
shared:
  input_file_glob:
    - /data/worldview/**/*.tif
  dem_file_path: /data/dem.tif
  epsg: 6635
  nodata_value: -9999
  dtype: int16
  log_to_console: false
  skip_existing: true
  # output_dir: ../Processed
  # temp_dir: ./temp

workflow:
  last_run_step: raw
  run_fetch_atmosphere: false
  save_fetch_atmosphere: temp
  run_atmospheric_correction: true
  save_atmospheric_correction: temp
  run_orthorectification: true
  save_orthorectification: temp
  run_pansharpen: true
  save_pansharpen: temp
  run_cloud_mask: true
  save_cloud_mask: output
  run_alignment: false
  save_alignment: temp
  run_radiometric_normalization: false
  save_radiometric_normalization: temp

atmospheric_correction:
  atmospheric_method: py6s
  py6s_atmosphere_profile: user
  py6s_aerosol_profile: maritime
  py6s_aot550: 0.2
  py6s_water_vapor: 2.5
  py6s_ozone: 0.3
  py6s_output_scale_factor: 10000.0
  py6s_output_dtype: int16
  py6s_use_imd_radiance_calibration: true
  py6s_use_worldview_gain_offset_adjustment: true
  fetch_atmosphere:
    fetch_atmosphere_source: nasa_power
    fetch_atmosphere_grid_size: 3
    fetch_atmosphere_search_days: 1
    fetch_atmosphere_timeout_s: 30.0
  # envi_engine_path: "/path/to/taskengine.exe"
  # py6s_executable: /usr/local/bin/sixsV1.1
  # flaash_param_MODTRAN_RES: 15.0

orthorectification:
  orthorectification_output_suffix: _ortho

pansharpen:
  pansharpen_output_suffix: _pansharpen

cloud_mask:
  cloud_mask_method: omnicloudmask
  cloud_buffer_pixels: 10
  cloud_mask_inference_resolution_m: 10.0

alignment:
  alignment_fixed_image: /data/fixed/reference.tif
  alignment_registration_mode: structural_wv3_lidar

radiometric_normalization:
  radiometric_normalization_method: spectralmatch
  # See: https://spectralmatch.github.io/spectralmatch/api/pipeline/
  # match_shared_window_size: 1024
```

## Field Notes

- `input_file_glob` is the main discovery input and is searched recursively for `.tif` files.
- YAML sections are flattened by the CLI, so nested sections are for readability and organization.
- `output_dir` is the base output location used when a step has `save_<step>: output`.
  - if unset, the default output base is `../Processed` from the MUL image folder
  - relative output paths resolve from the MUL image folder
- `temp_dir` is optional.
  - if unset, `vhr-worldview` creates a real temporary directory with Python `tempfile`
  - if set, `save_<step>: temp` writes under that directory
- every raster step has a matching `save_<step>` setting.
  - `temp` writes to `<temp_dir>/<step_name>`
  - `output` writes to `<output_dir>/<step_name>`
  - any other value is treated as a custom folder path, and the step name is appended to it
- `skip_existing: true` skips a whole scene when the final output TIFF already exists.
  - per-step output planning still skips only the missing outputs within each step
- `last_run_step` lets the workflow resume from an already completed raster stage.
- `run_fetch_atmosphere` controls the optional prefetch step for external atmosphere values.
- `run_radiometric_normalization` is the last step in the chain when enabled.
- by default, only `save_cloud_mask` uses `output`.
  - the other `save_*` values default to `temp`
- there is no separate final-copy stage.
  - the scene's final raster is whatever the last enabled raster step produced
  - the scene metadata JSON is written beside that raster
- `envi_engine_path` is only required when `atmospheric_method: flaash`.
- Py6S output defaults to scaled reflectance (`int16` with factor `10000`).
  - `reflectance = pixel_value / py6s_output_scale_factor`
- Built-in cloud masking infers the mask from orthorectified multispectral imagery and applies it to the current workflow output.
- `cloud_mask_omnicloud_kwargs_json` can be set in YAML as a mapping or from the CLI as a JSON string.
- `alignment_registration_mode: structural_wv3_lidar` mirrors the WV/LiDAR structural workflow.
  - despite the name, this mode is a generic structural edge-based registration path
- any `match_*` config keys are passed directly into the spectralmatch pipeline after stripping the `match_` prefix
- any `flaash_param_*` config keys are passed directly into FLAASH after stripping the prefix and uppercasing the key

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
