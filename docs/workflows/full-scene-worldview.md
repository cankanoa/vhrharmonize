# Full-Scene Workflow

This describes the execution behavior of `vhr-worldview`.

## High-Level Flow

1. Parse YAML + CLI overrides
2. Validate required inputs and compatibility rules
3. Discover scenes from input directory structure (optional `Root*.txt` overrides)
4. Read scene/photo metadata and overrides
5. For each scene:
   - run atmospheric correction (`py6s` default, `flaash` optional) or reuse existing atmospheric/ortho inputs
   - orthorectify multispectral and panchromatic imagery
   - pansharpen
   - optionally cloud-mask final scene
   - write final output in `output_dir`

## Output Model

- Intermediates live in temporary scratch directories.
- Final deliverables are full-scene rasters.
- For Py6S runs, atmospheric output is surface reflectance.
  - Default storage is scaled reflectance (`int16`, factor `10000`).

## Py6S Metadata Behavior (WorldView)

- Py6S geometry uses IMD `meanOffNadirViewAngle` directly for view zenith.
- Optional IMD calibration (`py6s_use_imd_radiance_calibration`) applies per-band
  `absCalFactor/effectiveBandwidth` before Py6S correction.

## Operational Recommendations

- Keep one config YAML per dataset family.
- Use output suffixes to separate processing runs.
- Store generated metrics/reports in `outputs/`.
