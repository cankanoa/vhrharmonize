# Full-Scene Workflow

This describes the execution behavior of `vhr-worldview`.

## High-Level Flow

1. Parse YAML + CLI overrides
2. Validate required inputs and compatibility rules
3. Discover scenes from input directory structure (optional `Root*.txt` overrides)
4. Read scene/photo metadata and overrides
5. For each scene:
   - run atmospheric correction (`py6s` default, `flaash` optional) when `run_atmospheric_correction: true`
   - orthorectify multispectral and panchromatic imagery
   - pansharpen when `run_pansharpen: true`
   - optionally infer cloud mask from orthorectified MS (default inference at 10 m), then apply it to current workflow output
   - optionally run alignment (`run_alignment: true`) against a fixed/reference raster
   - write final output in `output_dir`

## Output Model

- Intermediates live in temporary scratch directories.
- Final deliverables are full-scene rasters.
- When alignment is enabled, an additional aligned raster is written per scene
  using `alignment_output_suffix`.
- Each final raster also writes a sidecar metadata report:
  - `<scene_basename><output_suffix>_metadata.json`
  - contains input file paths, workflow settings, and effective atmospheric values used.
- For Py6S runs, atmospheric output is surface reflectance.
  - Default storage is scaled reflectance (`int16`, factor `10000`).

## Py6S Metadata Behavior (WorldView)

- Py6S geometry uses IMD `meanOffNadirViewAngle` directly for view zenith.
- Optional IMD calibration (`py6s_use_imd_radiance_calibration`) applies per-band
  `absCalFactor/effectiveBandwidth` before Py6S correction.

## Alignment Notes

- Alignment runs after the final scene image is produced (after cloud masking if enabled).
- The moving image for alignment is the workflow final scene output.
- `alignment_clip_fixed_to_moving: true` allows large fixed rasters by using only overlap extent.
- For WV/LiDAR intensity use, prefer:
  - `alignment_registration_mode: structural_wv3_lidar`
  - `alignment_moving_band_index: 6`
  - `alignment_fixed_band_index: 0`
  - `alignment_no_tiling: true`

## Operational Recommendations

- Keep one config YAML per dataset family.
- Use output suffixes to separate processing runs.
- Store generated metrics/reports in `outputs/`.
