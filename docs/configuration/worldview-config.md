# WorldView Config

The main workflow config is the YAML file used by `vhr-worldview`.

Use the tracked template:

- [configs/example.worldview.yml](https://github.com/cankanoa/vhrharmonize/blob/main/configs/example.worldview.yml)

The file is organized into these sections:

- `shared`: input discovery, temp/output roots, DEM settings, concurrency, cloud-cover filtering, reuse controls
- `workflow`: which steps run, where each step saves, whether overviews are built, and filename suffixes
- `atmospheric_correction`: atmospheric method-specific settings
- `orthorectification`: orthorectification-specific settings
- `pansharpen`: pansharpening-specific settings
- `cloud_mask`: cloud masking settings
- `alignment`: pairwise alignment settings
- `radiometric_normalization`: SpectralMatch settings

Use the YAML template as the source of truth for current keys and examples. The CLI parser can still override any of these values at runtime.
