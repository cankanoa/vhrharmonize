# `vhr-pansharpen-orthos`

Pansharpenes existing orthorectified multispectral and panchromatic rasters.

## Usage

```bash
vhr-pansharpen-orthos \
  --mul-ortho /data/mul_ortho.tif \
  --pan-ortho /data/pan_ortho.tif \
  --output /data/out/pansharp.tif
```

## Notes

- Input paths must exist.
- Optional `--nodata-value` sets output NoData.
