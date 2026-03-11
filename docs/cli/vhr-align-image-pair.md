# `vhr-align-image-pair`

Align a moving image (A) to a fixed/reference image (B) using elastix registration.

Default behavior is the WV/LiDAR structural workflow on full extent:

- `registration_mode: structural_wv3_lidar`
- tiling disabled
- fixed image clipped to moving extent
- mutual valid mask enforced
- output written on moving-image grid
- registration band index: `0` (0-based) unless you set `--moving-band-index`

The transform is estimated per tile on a single selected band, then applied to all bands for output.

## Requirements

- `itk-elastix` installed (for example: `pip install -e ".[elastix]"`)
- both rasters must have a defined CRS
- both rasters must currently use the same CRS
- both rasters should spatially overlap in that CRS

## Example

```bash
vhr-align-image-pair \
  --moving-image /data/image_a_cloudmasked.tif \
  --fixed-image /data/image_b_lidar.tif \
  --output-image /data/out/image_a_aligned_to_b.tif
```

## Important options

- `--band-index`: 0-based band index used for registration metric
- `--moving-band-index`: 0-based moving-image band index used for registration metric
- `--fixed-band-index`: 0-based fixed-image band index used for registration metric
- `--tile-size`: tile size in pixels
- `--tile-buffer`: overlap/buffer in pixels around each tile
- `--no-tiling` / `--tiling`: run on full extent or per-tile
- `--parameter-map`: elastix default map (`rigid`, `affine`, `bspline`, ...)
- `--parameter-file`: custom elastix parameter file (repeatable)
- `--moving-nodata`, `--fixed-nodata`: nodata overrides for mask generation
- `--output-nodata`: output nodata override
- `--registration-mode`: `default` or `structural_wv3_lidar`
- `--clip-fixed-to-moving`: limit fixed domain to overlap with moving image
- `--enforce-mutual-valid-mask`: constrain both masks to shared valid area
- `--output-on-moving-grid`: write result on moving-image grid/resolution
- `--keep-temp-dir`: preserve temporary tile artifacts for debugging

## Registration Modes

- `default`: raw-band elastix registration. Use this when both rasters are comparable radiometrically and you want a generic registration path.
- `structural_wv3_lidar`: intended for optical-to-LiDAR alignment. It builds structural proxy images on a common fixed-grid ROI and runs a chained `translation -> rigid` transform estimate, then applies that geometric transform to the original moving bands.

## Recommended WV/LiDAR Usage

```bash
vhr-align-image-pair \
  --moving-image /data/worldview_cloudmasked.tif \
  --fixed-image /data/mean_intensity_mosaic.tif \
  --output-image /data/worldview_aligned.tif \
  --registration-mode structural_wv3_lidar \
  --moving-band-index 6 \
  --fixed-band-index 0 \
  --no-tiling \
  --clip-fixed-to-moving \
  --enforce-mutual-valid-mask \
  --output-on-moving-grid
```

## Real-World Example

With current defaults, the following WV/LiDAR alignment only needs the non-default
moving band selection:

```bash
vhr-align-image-pair \
  --moving-image /mnt/s/Satellite_Imagery/Big_Island/Processed_cloudmasked/17SEP06212820-M1BS-200011893447_01_P001_cloud_masked.tif \
  --fixed-image /mnt/x/PROJECTS_2/Big_Island/ChangeHI_Trees/Dry_Forest/Data/Raster/mean_intensity/mean_intensity_mosaic.tif \
  --output-image /mnt/s/Satellite_Imagery/Big_Island/Processed_cloudmasked/17SEP06212820-M1BS-200011893447_01_P001_cloud_masked_aligned.tif \
  --moving-band-index 6
```

Expected completion output looks like:

```json
{
  "output_image_path": "/mnt/s/Satellite_Imagery/Big_Island/Processed_cloudmasked/17SEP06212820-M1BS-200011893447_01_P001_cloud_masked_aligned.tif",
  "total_tiles": 1,
  "successful_tiles": 1,
  "skipped_tiles": 0,
  "temp_dir": null
}
```

`Dataset has no geotransform, gcps, or rpcs. The identity matrix will be returned.`
messages can appear from transformix intermediate files and are not by themselves a failure.

## Masking behavior

For robust matching with cloud-masked imagery, registration uses masks:

- moving mask from valid data in moving tile
- fixed mask from valid fixed pixels intersected with moving valid footprint

This prevents cloud/no-data regions from dominating the registration metric.
