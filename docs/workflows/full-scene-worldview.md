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
   - optionally fetch atmosphere values
   - optionally run atmospheric correction (`py6s`, `flaash`, or `none`)
   - optionally orthorectify multispectral and panchromatic imagery
   - optionally pansharpen
   - optionally cloud-mask the current raster
   - optionally align to a fixed/reference raster
   - optionally run radiometric normalization
   - copy/write the final raster to the final output folder
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

- intermediate step outputs live under per-step folders
- when `temp_dir` is unset, a real temporary directory is created automatically
- when `temp_dir` is set, step folders are created under that directory unless a step-specific output dir overrides them
- final outputs default to `<scene root>/Processed` when `output_dir` is unset
- each final raster also writes a sidecar metadata report:
  - `<scene_basename><output_suffix>_metadata.json`

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
