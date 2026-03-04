# vhrharmonize

`vhrharmonize` is a Python library and CLI suite for preprocessing very high resolution (VHR) satellite imagery.
Today, the tested and supported end-to-end workflow is WorldView; support for additional providers and sensors (for example, Planet) is being added.

Core capabilities:

- Atmospheric correction
- RPC orthorectification
- Pansharpening
- Optional cloud masking with OmniCloudMask
- Metadata processing tools for automated workflows, including scene discovery

Primary execution path is the `vhr-worldview` CLI, which runs the full-scene workflow and writes final scene outputs.

## Documentation Map

- Installation and environment setup: `Getting Started`
- Operational command usage: `CLI`
- YAML templates and settings: `Configuration`
- End-to-end pipeline behavior: `Workflows`
- Function-level module index: `API Reference`

## Quick Command Preview

```bash
vhr-worldview --help
vhr-fetch-modis-water-vapor --help
vhr-cloudmask-raster --help
vhr-pansharpen-orthos --help
```

## Project Structure

- `vhrharmonize/`: reusable Python package
- `vhrharmonize/cli/`: packaged CLI implementations
- `configs/`: current config templates
- `docs/`: MkDocs source
