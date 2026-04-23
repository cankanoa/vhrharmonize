# Full-Scene Workflow

This describes the execution behavior of `vhr-worldview`.

## High-Level Flow

1. Parse YAML + CLI overrides
2. Validate required inputs and compatibility rules
3. Discover WorldView TIFFs from `input_file_glob`
4. Group TIFFs into scenes with one MUL image and one PAN image
5. Parse IMD metadata and map it into standardized metadata
6. For each scene:
   - optionally skip immediately when the final output TIFF already exists
   - resolve the DEM path
     - use the configured local DEM, or
     - download an OpenTopography SRTM GL1 ellipsoidal DEM into the temp dir when `dem_file_path=online`
   - optionally fetch atmosphere values
   - optionally run atmospheric correction (`py6s`, `flaash`, or `none`)
   - optionally orthorectify multispectral and panchromatic imagery
   - optionally pansharpen
   - optionally cloud-mask the current raster
   - optionally align to a fixed/reference raster
   - optionally run radiometric normalization
   - treat the last enabled raster step output as the final raster
   - write a scene metadata report JSON

## Scene Model

Each discovered WorldView scene is now treated as a single acquisition bundle:

- one multispectral image bundle
- one panchromatic image bundle
- provider-specific IMD metadata for each image
- standardized metadata for each image
- step outputs registered back onto the scene object

The scene root is determined from the image path:

- if the TIFF is under a folder ending in `_MUL` or `_PAN`, the parent of that folder becomes the scene root
- otherwise the TIFF directory itself becomes the scene root

## Output Model

- intermediate and final step outputs live under per-step folders
- when `temp_dir` is unset, a real temporary directory is created automatically
- when `temp_dir` is set, any step with `save_<step>=temp` writes under that directory
- `output_dir` is the base used when a step has `save_<step>=output`
  - if unset, it defaults to `../Processed` from the MUL image folder
  - relative paths resolve from the MUL image folder
- each step follows the same save rule:
  - `temp` -> `<temp_dir>/<step_name>`
  - `output` -> `<output_dir>/<step_name>`
  - custom path -> `<custom_path>/<step_name>`
- by default, only cloud-mask saves to `output`
- all other raster steps default to `temp`
- the final raster is the output of the last enabled raster step
- each final raster also writes a sidecar metadata report in the same folder

## Skip Existing

`skip_existing` works in two layers:

- scene-level skip: if the final output TIFF already exists, the whole scene is skipped
- step-level skip: within each step, only missing outputs are processed and existing outputs are reused

This allows old intermediate temp folders to be deleted without breaking completed scenes.

## Alignment Notes

- alignment runs after cloud masking and before radiometric normalization
- the moving image for alignment is the current workflow raster at that point
- `alignment_registration_mode: structural_wv3_lidar` is a structural edge-based mode intended for optical-to-LiDAR style alignment, but the implementation itself is not WorldView-specific

## Operational Recommendations

- keep one config YAML per dataset family
- use output suffixes to separate processing runs
- install only the extras needed for the steps you plan to run
