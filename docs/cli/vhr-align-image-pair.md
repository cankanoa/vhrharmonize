# `vhr-align-image-pair`

Align a moving image (A) to a fixed/reference image (B) using elastix registration.

Default behavior is tile-wise registration for large rasters:

- tiling enabled
- tile size: `1000` pixels
- tile buffer/overlap: `100` pixels
- registration band index: `0` (0-based)

The transform is estimated per tile on a single selected band, then applied to all bands for output.

## Example

```bash
vhr-align-image-pair \
  --moving-image /data/image_a_cloudmasked.tif \
  --fixed-image /data/image_b_lidar.tif \
  --output-image /data/out/image_a_aligned_to_b.tif
```

## Important options

- `--band-index`: 0-based band index used for registration metric
- `--tile-size`: tile size in pixels
- `--tile-buffer`: overlap/buffer in pixels around each tile
- `--no-tiling`: disable tiling and run on full extent
- `--parameter-map`: elastix default map (`rigid`, `affine`, `bspline`, ...)
- `--parameter-file`: custom elastix parameter file (repeatable)
- `--moving-nodata`, `--fixed-nodata`: nodata overrides for mask generation
- `--output-nodata`: output nodata override
- `--keep-temp-dir`: preserve temporary tile artifacts for debugging

## Masking behavior

For robust matching with cloud-masked imagery, registration uses masks:

- moving mask from valid data in moving tile
- fixed mask from valid fixed pixels intersected with moving valid footprint

This prevents cloud/no-data regions from dominating the registration metric.
