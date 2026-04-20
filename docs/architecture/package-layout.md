# Package Layout

The package now follows a layered structure to support multi-provider VHR preprocessing.

## Layers

- `vhrharmonize/preprocess`: step-focused preprocessing modules
  - atmospheric correction
  - alignment
  - orthorectification
  - pansharpening
  - cloudmasking
  - coregistration
  - radiometric normalization
- `vhrharmonize/providers`: provider adapters
  - `worldview.py`: WorldView file discovery and IMD parsing
  - `planet.py`: placeholder provider module
  - `standardized.py`: provider-neutral metadata model
- `vhrharmonize/io`: geospatial IO plus output-planning helpers
- `vhrharmonize/cli`: command entrypoints

## New Extension Points

- Atmospheric correction backends:
  - `preprocess.atmospheric_correction.atmospheric_correction(...)`
  - `preprocess.atmospheric_correction.FLAASHCorrector`
  - `preprocess.atmospheric_correction.Py6SCorrector`
- Relative radiometric normalization:
  - `preprocess.radiometric_normalization.radiometric_normalization(...)`

## Packaging

- Type marker: `vhrharmonize/py.typed`
- Extras:
  - `.[fetch-atmosphere]`
  - `.[cloud]`
  - `.[py6s]`
  - `.[flaash]`
  - `.[orthorectification]`
  - `.[pansharpen]`
  - `.[align]`
  - `.[radiometric-normalization]`
  - `.[docs]`
  - `.[all]`
