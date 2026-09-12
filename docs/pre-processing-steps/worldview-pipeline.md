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

If `temp_dir` is not set, the workflow creates a real temporary directory.

- `delete_temp_dir`: Delete workflow temp directories after all scene and aggregate processing succeeds. Defaults to true; false retains them. Directories containing source inputs or saved outputs are preserved.
- `run_from_existing`: Reuse VHRHarmonize step outputs. Seamline metadata retains incremental resume behavior; SpectralMatch handles its own reuse.
- `skip_existing`: Skip completed scenes.
- `run_from_existing_check_validity`: Check VHRHarmonize step outputs for GDAL raster readability, TIFF strip/tile bounds (including internal overviews), or atmosphere JSON syntax, even if others are missing. SpectralMatch outputs are excluded. JSON checks exclude required fields and atmospheric values. At each step, automatically remove invalid primary outputs before rerunning and refetch invalid atmosphere JSON. Preserve valid files, sidecars, and source inputs. Startup inspection and upload/download planning never delete files.
- `skip_existing_check_validity`: Apply the same checks before skipping scenes, without deleting files.
- `run_spectralmatch`: Call SpectralMatch after scene processing. VHRHarmonize does not inspect, validate, delete invalid outputs, or skip its outputs; use `match_shared_resume_from_steps` for reuse.
- `save_spectralmatch`: Output file or folder matching the last `match_steps` entry. Use a file for a merged mosaic, a folder for tiled merges or multi-raster steps, and a vector file for seamline steps.
- `calculate_overviews_spectralmatch`: Enable the last overview-capable entry in `match_steps`: joint coregistration, global regression, local block adjustment, or merge. Requires overview scales and errors if any explicit `match_*_build_overviews` is true. Uses SpectralMatch's default step order when `match_steps` is omitted.
- `delete_temp_steps_proactively`: Delete individual temp step files and sidecars once non-temp outputs pass the scene-skip validity checks, including skipped scenes with incomplete temp steps. Leaves directories in place and operates independently of `delete_temp_dir`. Defaults to false. Preserves source/configured inputs and rasters needed by pending aggregate steps. Temp-only scenes require a saved non-temp aggregate output before per-file cleanup.
- `concurrent_processing`: Group multiprocessing tasks by scene.

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
10. Optional SpectralMatch after scene-level raster steps complete

At the end of the scene run, the workflow writes a scene metadata JSON beside the final raster output. If the scene is skipped because of cloud cover filtering, it writes a short metadata JSON explaining the skip.
