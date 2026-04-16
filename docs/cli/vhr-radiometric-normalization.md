# `vhr-radiometric-normalization`

Run the shared radiometric normalization step directly from the CLI.

Install extras with:

```bash
pip install -e ".[radiometric-normalization]"
```

This command wraps the upstream `spectralmatch.pipeline(...)` function and passes
through `match_*` parameters after stripping the prefix.

## Example

```bash
vhr-radiometric-normalization \
  rrn \
  --input-image /data/image_a.tif \
  --input-image /data/image_b.tif \
  --output-image /data/out/normalized.tif
```

## Notes

- any `match_*` argument is forwarded into spectralmatch after stripping `match_`
- shared defaults such as temp dir, nodata, and output dtype can still be supplied from the vhrharmonize CLI layer
