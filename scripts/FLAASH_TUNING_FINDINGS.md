# FLAASH Tuning Findings (Current Project)

## Scope
These findings summarize what was tested and fixed while running the WorldView wrapper workflow on:

- Scene: `17OCT19211717-M1BS-050311750010_01_P001`
- Area: Hawaiʻi (same/near-identical footprint comparisons)

## Key Fixes Implemented
1. Corrected line-of-sight zenith handling to ENVI/FLAASH convention:
- `LOS_ZENITH = 180 - meanOffNadirViewAngle`
- Required because ENVI expects LOS zenith in `[90, 270]`.

2. Fixed WSL-to-Windows path handling for ENVI FLAASH:
- FLAASH-facing paths now convert from `/mnt/<drive>/...` to `D:\...` style before `run_flaash`.
- This resolved ENVI hydration/path errors.

3. Added per-photo override support in practice via `Root_WV.txt`:
- `ParamsOverridesPerPhoto` is applied after base defaults, so it takes precedence.
- Confirmed working for this scene.

4. Added tile-only pipeline with final-output-only write behavior:
- Intermediate rasters stay in scratch temp space.
- Final product is tile-clipped output only.

## Current Cloud Mask Position
Cloud masking in wrapper is currently:

1. FLAASH
2. Orthorectify
3. Pansharpen
4. Cloud mask (optional, `omnicloudmask`)
5. Tile clip/write

So cloud mask does **not** currently influence FLAASH parameter estimation.

## What Was Tested for Negative Reflectance
Matrix results were written to:

- `scripts/flaash_matrix_results.csv`

Metric:

- `neg_avg` = mean of per-band sampled negative pixel percentages.

Scenarios tested (same WV override `0.9541`):

1. `baseline_initial`: `25.2845`
2. `maritime_auto_30km`: `36.7676`
3. `highvis_rural_disabled`: `23.9478` (best)
4. `lowvis_rural_disabled`: `23.9478` (best, tie)
5. `tropical_maritime_disabled`: `25.7064`
6. `tropical_tropo_auto_20km`: `39.8744`

## Observed Improvement
Best tested result improved from:

- Baseline `25.2845` to best `23.9478`
- Absolute gain: `1.3367` percentage points

This is a real but modest gain.

## Most Likely Drivers of Scene-to-Scene Differences
Even for near-identical footprint:

1. Cloud and cloud-adjacency effects (especially large bright cloud presence).
2. Changes in water vapor/aerosol between acquisition times.
3. Geometry differences (sun and off-nadir view angle).
4. Dark-content fraction differences (shadow/water/low-albedo areas).
5. Atmospheric model mismatch (fixed settings can be wrong scene-by-scene).

## Practical Guidance Going Forward
1. Keep per-photo overrides in `Root_WV.txt` (already working).
2. Prefer non-auto aerosol settings for this scene family unless strong evidence says otherwise.
3. Continue scene-level WV estimation (monthly/daily MODIS where available).
4. Track negatives as a diagnostic, but expect some negatives in difficult scenes.
5. If product specs allow, clamp final reflectance `< 0` to `0` after QA.

## Important Clarification
`USE_AEROSOL = "Automatic Selection"` does **not** come from WorldView metadata.
It is an internal FLAASH retrieval from image radiance/model assumptions.
