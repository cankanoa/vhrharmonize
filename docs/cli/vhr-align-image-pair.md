# `vhr-align-image-pair`

Align a moving image (A) to a fixed/reference image (B) using `coregix`.

## Requirements

- `coregix` installed (for example: `pip install -e ".[align]"`)

## Example

```bash
vhr-align-image-pair \
  --moving-image /data/image_a_cloudmasked.tif \
  --fixed-image /data/image_b_lidar.tif \
  --output-image /data/out/image_a_aligned_to_b.tif \
  --split-factor 2 \
  --trim-edge-invalid \
  --edge-trim-depth 8 \
  --edge-trim-invalid-below -3000
```

## Important options

- `--band-index`: 0-based band index used for registration metric
- `--moving-band-index`: optional moving-image registration band
- `--fixed-band-index`: optional fixed-image registration band
- `--moving-nodata`, `--fixed-nodata`: optional nodata overrides
- `--output-nodata`: optional output nodata override
- `--min-valid-fraction`: required valid overlap fraction in the registration ROI
- `--temp-dir`: optional parent directory for coregix working files
- `--keep-temp-dir`: keep the coregix working directory
- `--use-edge-proxies` / `--no-use-edge-proxies`
- `--split-factor`: coregix split factor
- `--clip-fixed-to-moving` / `--no-clip-fixed-to-moving`
- `--output-on-moving-grid` / `--no-output-on-moving-grid`
- `--trim-edge-invalid` / `--no-trim-edge-invalid`
- `--edge-trim-depth`: number of edge pixels to trim
- `--edge-trim-detection-band-index`: band used to detect edge artifacts
- `--edge-trim-invalid-below`: threshold used to identify invalid edge values
- `--edge-trim-invalid-above`: optional upper invalid threshold
- `--enforce-mutual-valid-mask` / `--no-enforce-mutual-valid-mask`
- `--solve-resolution`: optional target pixel size for the registration solve

Expected completion output looks like:

```json
{
  "output_image_path": "/data/out/image_a_aligned_to_b.tif"
}
```
