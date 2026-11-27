# Individual SpectralMatch functions

These adapters import the installed SpectralMatch functions directly. Processing algorithms remain in that package. YAML uses the native parameter names, with explicit connections in `shared` and the processing steps in [the WorldView recipe](https://github.com/cankanoa/vhrharmonize/blob/main/configs/example.worldview.yml). Statistics functions and the old pipeline adapter are excluded.

All functions are available from `vhrharmonize` in Python and as `plugin: function_name` in YAML. Each has a generated `vhr function_name --config recipe.yml` command. Use explicit file lists in recipes for core dependency tracking, reuse, validation and HPC transfers. The full native parameter catalog, defaults and concise option comments are in the example recipe.

```python
from vhrharmonize import global_regression, merge_rasters

matched = global_regression(input_images=["a.tif", "b.tif"],
                            output_images=["a_global.tif", "b_global.tif"],
                            custom_nodata_value=-9999, output_dtype="float32")
merge_rasters(input_images=matched, output_image_path="mosaic.tif")
```

## joint_coregistration

::: vhrharmonize.plugins.joint_coregistration

## global_regression

::: vhrharmonize.plugins.global_regression

## local_block_adjustment

::: vhrharmonize.plugins.local_block_adjustment

## align_rasters

::: vhrharmonize.plugins.align_rasters

## create_footprints

SpectralMatch owns both raw footprint creation and optional initial mosaic coverage.
The VHR adapter passes its native arguments through. `footprint_source: metadata`
uses the GeoJSON `geometry` in each ordered `metadata_records` entry without
reading raster pixels; the default `calculate_bounds` polygonizes validity masks.
Scalar metadata attributes are copied into the output. `metadata_crs` defaults to
EPSG:4326, and `epsg` optionally selects the output CRS.

Set `outline_path` (and optionally `outline_layer`) to define the goal area.
`image_rank_function` and `image_rank_descending` assign overlapping area to the
preferred image. `extrude_distance` buffers each owned area in projected CRS units,
clipped to its source footprint, preserving overlap for matching and movement.
Only connected components containing original owned area are retained: growth
cannot pick up detached source islands, but may join existing owned components.
Buffering also extends beyond the outline. Empty owners are omitted, and VHR
retains only the corresponding workflow scenes in their original order.

The WorldView recipe invokes `plugin: create_footprints` three times:

- `polygonize_images` runs per scene after cloud masking. It writes an individual
  `_footprint.gpkg` into the processed folder and validates cached vectors on reruns.
- `standardize_mosaic_coverage` uses `footprint_source: vector_bounds` with
  `input_images: collect:polygonized_footprint_path`. It reads those saved geometries,
  applies ranking and a margin, then `mask_mosaic_coverage` crops the rasters.
- `create_seamline_footprints` recalculates valid-data polygons after alignment and
  radiometric matching, directly before `postprocess_footprints` consumes them.

In vector mode `input_images` selects saved vectors, one original image per file.
`image_field_name` selects the ID column, while `image_path` contains the original
raster path. These identifiers and saved scalar attributes are retained; raster
files are not required or opened. `vector_layer` selects a layer, or defaults to
the first. Multiple features with the same image ID are combined before ranking.
Explicit `metadata_records` or CSV values override saved attributes with matching
names. VHR tracks and stages the vector files as dependencies and removes scenes
using their preserved image IDs.

The initial coverage recipe ranks each connected polygon component using:

```yaml
param:image_rank_function: >-
  5 * log(max({stat:area_crs}, 1) / 1000000) / log(2)
  + 200 * {stat:gap_coverage_fraction}
  - ({var:view_zenith} + {var:solar_zenith})
```

`stat:area_crs` is calculated from each component's geometry after outline clipping
and before ownership subtraction/buffering. Holes are subtracted. It uses squared
CRS units; the example uses metres, making `/1000000` an area in km². Doubling a
component's area offsets five degrees of worse combined view/solar zenith.
MultiPolygons are ranked as separate components with the same image ID and metadata,
so a small island cannot borrow a larger component's area bonus. The statistic is
not a table field; `var:area_crs` remains an ordinary arbitrary column reference.
Ranking without an area-statistic reference continues to operate per image.

`stat:gap_coverage_fraction` enables iterative gap filling in `create_footprints`.
It is zero for the first selection, which therefore uses the existing area/angle
preference. After each owner is chosen, the algorithm finds the largest connected
gap in the remaining available coverage and recomputes candidate scores. The
fraction is the candidate's intersection area with that gap divided by the gap's
area, from 0 to 1. Candidates that cover none of the gap are excluded from that
selection. A selected candidate owns all of its still-uncovered geometry, which
can also fill other gaps. The goal excludes outline areas with no source coverage.

At weight 200, covering 20 percentage points more of a gap adds 40 score points,
offsetting 40 degrees of worse combined zenith when the other terms are equal.
The area term remains each original clipped component's area, and image summaries
remain fixed across selections. Buffering happens after ownership and cannot
change gap scores. `weighted_score` records the score at selection and
`weighted_rank` records selection order; scores need not decrease across ranks.
This is a weighted tradeoff, not a guarantee of the fewest images. Expressions
without the gap statistic retain the existing one-time ranking. The gap statistic
requires a SpectralMatch version supporting it and is unavailable to weighted
seamlines and Markov image-quality ranking, which lack this coverage context.

Raw defaults preserve overlaps; outline/ranking/buffering recomputes coverage as
a group and replaces the polygon file only after success.

Set `param:do_postprocess_footprints: true` to call `postprocess_footprints` on
each owned piece before claiming its area. Optional
`param:postprocess_footprints_params: {}` forwards keyword arguments to that
function; omitted options inherit its current defaults. Input/output paths are
managed internally. Removed areas remain available to later candidates, and gap
ranking uses the cleaned ownership. With `stat:area_crs`, each ranked component
is processed separately. A candidate removed entirely is skipped. The overlap
buffer and source clipping still happen afterward, retaining only components
containing owned area. This prevents new detached islands, but does not guarantee
that every retained component contains a raster pixel.

::: vhrharmonize.plugins.create_footprints

## postprocess_footprints

::: vhrharmonize.plugins.postprocess_footprints

## voronoi_center_seamline

::: vhrharmonize.plugins.voronoi_center_seamline

## weighted_seamline

::: vhrharmonize.plugins.weighted_seamline

## markov_triangles

The adapter uses bounded raster sampling inside SpectralMatch's node workers.
Each requested window is average-resampled with GDAL's chunked warper using a
32 MiB warp budget and a separate 32 MiB sample cache. This avoids coarse warped
VRT reads expanding to nearly an entire full-resolution, multiband source in RAM.
The native grid, band selection, NoData handling, fallback sampling and both
default cost terms are retained. The budgets exclude GDAL's separate block cache
and the mesh/solver arrays. The override runs in serial, process-pool and supplied
executor workers without changing the installed SpectralMatch files.

::: vhrharmonize.plugins.markov_triangles

## mask_rasters

`param:seperate_multipolygons: true` writes one cropped raster for each disconnected
retained polygon. For an output base `scene_masked.tif`, files are
`scene_masked_1.tif`, `scene_masked_2.tif`, etc. Numbering restarts for each input
and follows polygon bounds in increasing coordinate order. Overlapping/touching
mask features are dissolved before splitting; interior holes remain masked.
An exclude mask splits the remaining raster area. The default is `false`.

The plugin reports the resulting scene changes directly to core. Each returned
numbered file becomes an independent scene with its own `current_image_paths`,
`file_path`, scene ID and basename, and an empty naming suffix. Parent metadata
and protected source paths are inherited. No additional import step is needed.

```yaml
mask_mosaic_coverage:
  plugin: mask_rasters
  core:run: true
  param:input_images: collect:current_image_paths
  var:suffix: expr:var.suffix & '_coverage'
  var:current_image_paths: expr:const.temp_dir & '/coverage/' & var.basename & var.suffix & '.tif'
  param:output_images: collect:current_image_paths
  core:requires: const:coverage_mask
  param:vector_mask: [include, const:coverage_mask, const:spectralmatch.image_field_name]
  param:seperate_multipolygons: true
  param:resume_from_outputs: "yes"
```

Downstream steps continue using `collect:current_image_paths` and `var:basename`.
`returned:$` still exposes the ordinary flat output-path list if requested.
Only actual returned outputs are used, so stale numbered files are ignored.
Native reuse restores the same scene mapping without rewriting completed files.
An empty component list raises a clear no-coverage error.

Core tracks numbered files for validation, collisions, overviews, cleanup and
transfers. Scene replacement happens after successful execution, never during
planning. HPC preparation of a downstream suffix must first run through the
splitting step with `run_to_step_before_prepare`; a terminal masking-only job can
be staged directly. Connected narrow arms are not split off.

::: vhrharmonize.plugins.mask_rasters

## merge_rasters

::: vhrharmonize.plugins.merge_rasters

## merge_vectors

::: vhrharmonize.plugins.merge_vectors

## band_math

::: vhrharmonize.plugins.band_math

## create_cloud_mask_with_omnicloudmask

::: vhrharmonize.plugins.create_cloud_mask_with_omnicloudmask

## process_raster_values_to_vector_polygons

::: vhrharmonize.plugins.process_raster_values_to_vector_polygons

## compute_overviews

::: vhrharmonize.plugins.compute_overviews

## search_paths

::: vhrharmonize.plugins.search_paths

## create_paths

::: vhrharmonize.plugins.create_paths

## match_paths

::: vhrharmonize.plugins.match_paths
