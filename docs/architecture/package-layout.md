# Package Layout

The package now follows a layered structure to support multi-provider VHR preprocessing.

## Layers

- `vhrharmonize/preprocess`: step-focused preprocessing modules
  - atmospheric correction
  - orthorectification
  - pansharpening
  - cloudmasking
  - registration
  - tiling
  - radiometric normalization (RRN/ARN)
- `vhrharmonize/providers`: provider adapters (`worldview`, `planet` scaffold)
- `vhrharmonize/pipelines`: orchestration layer for scene workflows
- `vhrharmonize/cli`: command entrypoints

## New Extension Points

- Atmospheric correction backends:
  - `preprocess.atmospheric_correction.atmospheric_correction(...)`
  - `preprocess.atmospheric_correction.FLAASHCorrector`
  - `preprocess.atmospheric_correction.Py6SCorrector`
- Relative radiometric normalization:
  - `preprocess.radiometric_normalization.apply_rrn(...)`
  - `preprocess.radiometric_normalization.spectralmatch_rrn(...)` (scaffold backend)
- Absolute radiometric normalization:
  - `preprocess.radiometric_normalization.apply_arn(...)` (scaffold)

## Packaging

- Type marker: `vhrharmonize/py.typed`
- Extras:
  - `.[py6s]` (backward-compat alias; Py6S is in base dependencies)
  - `.[spectralmatch]`
  - `.[cloud]`
  - `.[docs]`
  - `.[all]`
