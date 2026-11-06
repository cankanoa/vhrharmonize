# Individual SpectralMatch functions

These adapters import the installed SpectralMatch functions directly. Processing algorithms remain in that package. YAML uses the native parameter names, with explicit connections in `setup_spectralmatch` and subsequent steps in [the WorldView recipe](https://github.com/cankanoa/vhrharmonize/blob/main/configs/example.worldview.yml). Statistics functions and the old pipeline adapter are excluded.

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

The WorldView recipe invokes `plugin: create_footprints` twice: the early
`standardize_mosaic_coverage` call uses IMD geometry, ranking and a margin, then
`mask_mosaic_coverage` crops the rasters once before coregistration and matching.
The later `create_footprints` call retains raw valid-data footprints for final
seamlines. Raw defaults preserve overlaps; outline/ranking/buffering recomputes
coverage as a group and replaces the polygon file only after success.

The initial coverage recipe favors the lowest combined raw `view_zenith` and
`solar_zenith`, with equal weight per degree and `image_rank_descending: true`.
**Higher expression values rank on top by default.** Negating the raw sum makes
the lowest combined angle produce the highest ranking value:

```text
-({var:view_zenith} + {var:solar_zenith})
```

There is no standardization or cloud-cover term in this recipe. Both angles at
zero produce the best score, zero; larger combined angles produce lower scores.
The expression also works for a single image or constant-angle sets.

SpectralMatch's summary placeholders remain available for other ranking
expressions. `{var:field}` selects an image attribute. `{summary:field}` selects an
image-level summary across the candidates entering that ranking call, after
outline clipping and before ownership assignment. Supported summaries are
`mean`, `std`, `variance`, `min`, `max`, `median`, `sum`, `count`, and `range`.
`std` and `variance` use the population (`ddof=0`); each image counts equally.
Referenced summary fields must contain finite numeric values for every candidate.
The native `image_rank_function` default remains `None`, and the later raw
footprint call leaves it unset.

The former VHR `seamline_metadata` plugin has been removed. Use `input_images`,
`output_polygons`, `output_layer`, and `eight_connected` in place of its
`image_paths`, `output_path`, `layer`, and `calculate_bounds_eight_connected`.
Native raw-output reuse uses `resume_from_outputs`, replacing `reuse`. Install
the updated SpectralMatch checkout in both local and HPC environments.

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
