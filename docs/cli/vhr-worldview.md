# vhr-worldview

## Overview

`vhr-worldview` is the main end-to-end workflow for WorldView scenes. It discovers MUL and PAN pairs, reads IMD metadata, and runs the selected preprocessing steps in sequence.

## Installation

```bash
pip install "vhrharmonize[all]"
```

## Usage

```bash
cp configs/example.worldview.yml configs/worldview.yml

vhr-worldview --config-yaml configs/worldview.yml
```

```bash
vhr-worldview \
  --input-file-glob "/data/worldview/**/*.TIF" \
  --dem-file-path /data/dem.tif \
  --output-dir ../../processed \
  --atmospheric-method py6s \
  --run-alignment \
  --alignment-fixed-image /data/reference.tif
```

- Config and discovery: `--config-yaml`, `--input-file-glob`, `--filter-basename`
- DEM and projection: `--dem-file-path`, `--dem-online-api-key`, `--dem-online-source`, `--dem-online-api-endpoint`, `--dem-online-timeout-s`, `--epsg`
- Shared output/runtime: `--output-dir`, `--temp-dir`, `--keep-temp-dir`, `--nodata-value`, `--dtype`, `--log-to-console`, `--concurrent-processing`, `--overview-scales`
- Scene reuse and filtering: `--skip-existing/--no-skip-existing`, `--run-from-existing/--no-run-from-existing`, `--max-cloud-cover-to-process`
- Per-step enable/save/overview flags: `--run-fetch-atmosphere`, `--save-fetch-atmosphere`, `--run-atmospheric-correction`, `--save-atmospheric-correction`, `--calculate-overviews-atmospheric-correction`, `--run-orthorectification`, `--save-orthorectification`, `--calculate-overviews-orthorectification`, `--run-pansharpen`, `--save-pansharpen`, `--calculate-overviews-pansharpen`, `--run-cloud-mask`, `--save-cloud-mask`, `--calculate-overviews-cloud-mask`, `--run-alignment`, `--save-alignment`, `--calculate-overviews-alignment`, `--run-seamline-metadata`, `--save-seamline-metadata`, `--run-radiometric-normalization`, `--save-radiometric-normalization`, `--calculate-overviews-radiometric-normalization`
- Per-scene step save paths support: `$temp`, `$temp/...`, `$output`, `$output/...`, `./relative/to/mul`, `/absolute/path`, `relative/to/current/working/directory`
- Aggregate single-output save paths support: `$temp/file`, `/absolute/file`, `relative/to/current/working/directory/file`
- Output suffixes: `--fetch-atmosphere-output-suffix`, `--atmospheric-correction-output-suffix`, `--orthorectification-output-suffix`, `--orthorectification-pan-output-suffix`, `--pansharpen-output-suffix`, `--cloud-mask-output-suffix`, `--cloud-mask-mask-suffix`, `--alignment-output-suffix`
- Atmospheric method and FLAASH controls: `--atmospheric-method`, `--envi-engine-path`, `--flaash-dem-ground-percentile`, `--flaash-modtran-atm`, `--flaash-modtran-aer`, `--flaash-use-aerosol`, `--flaash-default-visibility`, `--skip-flaash`, `--existing-flaash-input`, plus any `--flaash-param-*` passthrough argument
- Py6S controls: `--py6s-atmosphere-profile`, `--py6s-aerosol-profile`, `--py6s-aot550`, `--py6s-visibility`, `--py6s-water-vapor`, `--py6s-ozone`, `--py6s-output-scale-factor`, `--py6s-output-dtype`, `--py6s-executable`, `--py6s-use-imd-radiance-calibration/--no-py6s-use-imd-radiance-calibration`, `--py6s-use-worldview-gain-offset-adjustment/--no-py6s-use-worldview-gain-offset-adjustment`, `--py6s-auto-atmos-source`, `--py6s-auto-atmos-grid-size`, `--py6s-auto-atmos-search-days`, `--py6s-auto-atmos-timeout-s`, `--py6s-auto-atmos-power-endpoint`
- Fetch atmosphere controls: `--fetch-atmosphere-source`, `--fetch-atmosphere-grid-size`, `--fetch-atmosphere-search-days`, `--fetch-atmosphere-timeout-s`, `--fetch-atmosphere-power-endpoint`, `--fetch-atmosphere-ee-project`, `--fetch-atmosphere-authenticate/--no-fetch-atmosphere-authenticate`, `--fetch-atmosphere-env-file`, `--fetch-atmosphere-hours-window`
- Orthorectification controls: `--existing-mul-ortho-input`, `--existing-pan-ortho-input`, `--orthorectification-rpc-refinement-geojson`
- Cloud mask controls: `--cloud-mask-command`, `--cloud-mask-method`, `--cloud-mask-red-band-index`, `--cloud-mask-green-band-index`, `--cloud-mask-nir-band-index`, `--cloud-mask-classes`, `--cloud-buffer-pixels`, `--cloud-mask-inference-resolution-m`, `--cloud-mask-omnicloud-kwargs-json`
- Alignment controls: `--alignment-fixed-image`, `--alignment-band-index`, `--alignment-moving-band-index`, `--alignment-fixed-band-index`, `--alignment-moving-nodata`, `--alignment-fixed-nodata`, `--alignment-output-nodata`, `--alignment-min-valid-fraction`, `--alignment-split-factor`, `--alignment-clip-fixed-to-moving/--no-alignment-clip-fixed-to-moving`, `--alignment-output-on-moving-grid/--no-alignment-output-on-moving-grid`, `--alignment-trim-edge-invalid/--no-alignment-trim-edge-invalid`, `--alignment-edge-trim-depth`, `--alignment-edge-trim-detection-band-index`, `--alignment-edge-trim-invalid-below`, `--alignment-edge-trim-invalid-above`, `--alignment-enforce-mutual-valid-mask/--no-alignment-enforce-mutual-valid-mask`, `--alignment-use-edge-proxies/--no-alignment-use-edge-proxies`, `--alignment-solve-resolution`
- Seamline metadata controls: `--seamline-metadata-layer`, `--seamline-metadata-image-field-name`, `--seamline-metadata-footprint-source`, `--seamline-metadata-calculate-bounds-eight-connected/--no-seamline-metadata-calculate-bounds-eight-connected`
- Radiometric normalization controls: `--radiometric-normalization-method`, `--radiometric-normalization-kwargs-json`, `--group-by-basename`, plus current `--match-*` SpectralMatch passthrough arguments such as `--match-steps`

When `log_to_console` is enabled, processing uses this progress format:

```text
START CLOUD_MASK:
[17DEC08211836-M1BS-016445318010_01_P012] Start | in=image.tif | out=image_cloudmasked.tif, image_cloudmask.tif
[17DEC08211836-M1BS-016445318010_01_P012] Completed 3/123
```

Glob discovery announces `start glob_matches:` before searching; scene discovery uses `start load_worldview_scenes:`. Only processing-step headings are uppercase. Each scene step announces its start before processing, and footprint metadata reports each image as it is calculated or reused. `Completed X/Y` uses the scene's position in the input list and total scene count; parallel workers can finish out of order. Scene steps run within each worker, so step headings can repeat and interleave. An exception does not emit a successful completion message. Standalone preprocessing calls use the same format with `1/1`; SpectralMatch controls its own internal progress messages. Footprint workers return geometry and image basename; the main process writes each record as it arrives. Existing matching records count as already processed.
