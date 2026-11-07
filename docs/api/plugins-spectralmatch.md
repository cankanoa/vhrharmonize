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
  - ({var:view_zenith} + {var:solar_zenith})
```

`stat:area_crs` is calculated from each component's geometry after outline clipping
and before ownership subtraction/buffering. Holes are subtracted. It uses squared
CRS units; the example uses metres, making `/1000000` an area in km². Doubling a
component's area offsets five degrees of worse combined view/solar zenith.
MultiPolygons are ranked as separate components with the same image ID and metadata,
so a small island cannot borrow a larger component's area bonus. The statistic is
not a table field; `var:area_crs` remains an ordinary arbitrary column reference.
Ranking without a geometry-statistic reference continues to operate per image.

Raw defaults preserve overlaps; outline/ranking/buffering recomputes coverage as
a group and replaces the polygon file only after success.

::: vhrharmonize.plugins.create_footprints

## postprocess_footprints

::: vhrharmonize.plugins.postprocess_footprints

## voronoi_center_seamline

::: vhrharmonize.plugins.voronoi_center_seamline

## weighted_seamline

::: vhrharmonize.plugins.weighted_seamline

## markov_triangles

::: vhrharmonize.plugins.markov_triangles

## mask_rasters

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
