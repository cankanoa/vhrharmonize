# WorldView Pipeline

## B1 Scene Structure

The streamlined workflow currently targets WorldView Basic 1B scenes. A scene is treated as one matched multispectral and panchromatic acquisition pair.

The workflow expects:

- A multispectral TIF
- A panchromatic TIF
- Matching IMD metadata files containing the scene corner coordinates
- RPC companions (such as RPB) when the raster needs them for orthorectification

Scene discovery groups files by basename and product identifiers so the MUL and PAN files for the same acquisition are processed together.

Bounds remain as IMD corner values in `standardized_metadata.source_metadata`. `materialize_scene_bounds(source_metadata, epsg=4326)` constructs a Shapely polygon on demand, connecting UL, UR, LR, and LL. DEM sampling reprojects that geometry into the DEM CRS. `package_bounds` uses the same polygon, while `calculate_bounds` continues to polygonize the processed raster mask. Missing or inconsistent IMD corners raise an error only when bounds are requested.

Local staging and HPC uploads share the same file selection from each image directory: the raster, IMD, RPC and other supported image metadata/sidecars. No separate footprint directory is needed.

## File Handling

The workflow plans outputs before each step runs. Every raster step builds its output name from the previous raster output name, so filenames reflect the full chain of steps that produced them.

There are two shared roots:

- `temp_dir`
- `output_dir`

Each step save location can point to:

- `$temp`
- `$temp/...`
- `$output`
- `$output/...`
- `./relative/to/mul`
- `/absolute/path`
- `relative/to/current/working/directory`

If `temp_dir` is not set, the workflow creates a real temporary directory. If `keep_temp_dir` is false, temp-saved files are deleted after a scene finishes.

The workflow also supports:

- step-level reuse with `run_from_existing`
- scene-level skipping with `skip_existing`
- GDAL readability checks for reused rasters and JSON-object validation for cached metadata with `run_from_existing_check_validity`
- GDAL readability checks before whole-scene skips with `skip_existing_check_validity`
- grouped per-scene multiprocessing with `concurrent_processing`

Atmosphere caches must contain readable JSON objects before reuse. Empty, truncated, or non-object caches are logged as invalid and fetched again, including when raster validity checks are disabled. Valid cached atmosphere data and completed raster outputs remain reusable.

Workflow JSON files are written to a temporary file in the destination directory, flushed and closed, then atomically replaced. A failed write preserves any previous destination and removes the temporary file. After freeing space following a disk-full failure, rerun the same configuration with resume enabled to recover invalid atmosphere caches.

## Processing Steps

The scene pipeline works in this order:

1. Scene discovery and IMD parsing
2. Optional atmosphere fetch
3. Atmospheric correction with Py6S, FLAASH, or no correction
4. Orthorectification of the multispectral raster
5. Orthorectification of the panchromatic raster when pansharpening is enabled
6. Pansharpening
7. Optional cloud masking
8. Optional alignment to a fixed raster
9. Optional seamline metadata vector after scene-level raster steps complete
10. Optional radiometric normalization after scene-level raster steps complete

At the end of the scene run, the workflow writes a scene metadata JSON beside the final raster output. If the scene is skipped because of cloud cover filtering, it writes a short metadata JSON explaining the skip.
