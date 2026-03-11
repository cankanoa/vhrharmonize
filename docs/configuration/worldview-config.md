# WorldView Config (`configs/worldview.example.yml`)

`configs/worldview.example.yml` is the tracked template for `vhr-worldview`.

Create your local config first:

```bash
cp configs/worldview.example.yml configs/worldview.yml
```

## Example

```yaml
input_dir:
  - /mnt/d/20171019_35cm_WV03_BAB_050311750010
dem_file_path: /mnt/d/dem/Hawaii_SRTM_GL1Ellip.tif
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
py6s_auto_atmos_source: nasa_power
py6s_auto_atmos_grid_size: 3
py6s_auto_atmos_search_days: 1
py6s_auto_atmos_timeout_s: 30.0
# py6s_executable: /usr/local/bin/sixsV1.1

# Required only if atmospheric_method: flaash
# envi_engine_path: "/mnt/c/Program Files/NV5/ENVI60/IDL90/bin/bin.x86_64/taskengine.exe"

epsg: 6635
nodata_value: -9999
scratch_dir: /mnt/d/test_data/scratch

output_dir: /mnt/d/test_data/out/scenes
output_suffix: _final
skip_existing: false

run_atmospheric_correction: true
run_pansharpen: true
run_cloud_mask: true
cloud_mask_method: omnicloudmask
cloud_buffer_pixels: 10
cloud_mask_inference_resolution_m: 10.0
cloud_mask_omnicloud_kwargs_json:
  inference_device: cuda
  mosaic_device: cpu
  batch_size: 1
  patch_size: 768
  patch_overlap: 192

# Optional final alignment
run_alignment: false
alignment_fixed_image: /mnt/d/lidar/mean_intensity_mosaic.tif
alignment_output_suffix: _aligned
alignment_moving_band_index: 0
alignment_fixed_band_index: 0
alignment_registration_mode: structural_wv3_lidar
alignment_parameter_map: rigid
alignment_no_tiling: true
alignment_tile_size: 1000
alignment_tile_buffer: 100
alignment_clip_fixed_to_moving: true
alignment_enforce_mutual_valid_mask: true
alignment_output_on_moving_grid: true
alignment_moving_nodata: -9999
alignment_fixed_nodata: -9999
# alignment_output_nodata: -9999
alignment_min_valid_fraction: 0.002
alignment_temp_dir: /tmp
alignment_keep_temp_dir: false
alignment_log_to_console: false
```

## Field Notes

- `input_dir` can be a list; command-line flags still override config values.
- `atmospheric_method` defaults to `py6s`.
- `envi_engine_path` is only required when `atmospheric_method: flaash`.
- Py6S output defaults to scaled reflectance (`int16` with factor `10000`).
  - `reflectance = pixel_value / py6s_output_scale_factor`
- WorldView IMD calibration is enabled by default (`py6s_use_imd_radiance_calibration: true`).
  - Uses IMD per-band `absCalFactor` and `effectiveBandwidth` for DN->radiance input to Py6S.
  - WV02/WV03 gain/offset table is controlled by `py6s_use_worldview_gain_offset_adjustment`.
  - Default is `true`.
  - Set `py6s_use_imd_radiance_calibration: false` to disable.
- Optional NASA auto-atmosphere mode:
  - `py6s_auto_atmos_source: nasa_power`
  - intended for `py6s_atmosphere_profile: user`
  - when `py6s_atmosphere_profile: user`, this auto-updates:
    - `py6s_aot550`
    - `py6s_water_vapor`
    - `py6s_ozone`
  - values are computed as bbox sample-grid mean on scene date (with configurable +/- day search).
- `scratch_dir` is used for scene intermediates.
- `output_dir` controls final scene output location.
- `skip_existing: true` skips scenes only when the final metadata JSON and final output TIFF already exist.
  - if alignment is enabled, the aligned TIFF must also exist
- `run_atmospheric_correction` defaults to `true`; set to `false` to skip atmospheric correction.
- `run_pansharpen` defaults to `true`; set to `false` to stop at orthorectified multispectral output.
- `run_cloud_mask` defaults to `true`; set to `false` to skip cloud masking without removing cloud settings.
- `cloud_mask_method` currently supports `omnicloudmask`.
- Built-in cloud masking always infers mask from orthorectified multispectral and applies it to current workflow output.
- `cloud_mask_inference_resolution_m` defaults to `10.0` so OmniCloudMask runs at coarser resolution suited to the model, then the mask is reprojected back to the workflow output grid.
- `cloud_mask_omnicloud_kwargs_json` accepts either:
  - a YAML mapping/dictionary (recommended in YAML config), or
  - a JSON string (same format used by CLI flag).
- `run_alignment` defaults to `false`; set `true` to align final scene output to a fixed raster.
- `alignment_fixed_image` is required when `run_alignment: true`.
- `alignment_moving_band_index` selects the moving (scene output) band for registration metric. Default is `0`; for many WV/LiDAR runs you may want to set this explicitly to `6`.
- `alignment_fixed_band_index` selects the fixed raster band for registration metric.
- `alignment_registration_mode: structural_wv3_lidar` mirrors the WV/LiDAR structural workflow.
- `alignment_clip_fixed_to_moving: true` allows using very large fixed rasters by clipping to overlap.
- `alignment_output_on_moving_grid: true` writes aligned output back to the original moving-image grid/resolution.

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
