# WorldView Pipeline

## B1 Scene Structure

The streamlined workflow currently targets WorldView Basic 1B scenes. A scene is treated as one matched multispectral and panchromatic acquisition pair.

The workflow expects:

- A multispectral TIF
- A panchromatic TIF
- The matching IMD metadata files
- The matching footprint files such as SHP

Scene discovery groups files by basename and product identifiers so the MUL and PAN files for the same acquisition are processed together.

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
- grouped per-scene multiprocessing with `concurrent_processing`

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
9. Optional radiometric normalization after scene-level raster steps complete

At the end of the scene run, the workflow writes a scene metadata JSON beside the final raster output. If the scene is skipped because of cloud cover filtering, it writes a short metadata JSON explaining the skip.
