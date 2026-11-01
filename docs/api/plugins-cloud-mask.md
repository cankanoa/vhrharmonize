Use `inference_image_path` to predict clouds from a separate raster, such as
orthorectified multispectral imagery, while `input_image_path` selects the raster
to mask, such as the pansharpened product. Omitting `inference_image_path` keeps
prediction and application on the same input. Band indices refer to the prediction
raster. The saved binary mask uses the inference grid; application resamples it
with nearest neighbor onto the target grid when `allow_mask_reprojection` is true
(the default). Output resolution and bands follow `input_image_path`.

::: vhrharmonize.plugins.cloud_mask
