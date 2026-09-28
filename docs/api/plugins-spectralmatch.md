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
