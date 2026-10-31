# vhr workflow

This command is generated from the public `run_workflow` function's signature and docstring. `--config`, `--config-dir` and `--dry-run` correspond directly to its Python arguments. Defaults and validation live in the Python API.

```python
from vhrharmonize import load_workflow, run_workflow

counts = run_workflow("configs/example.worldview.yml", dry_run=True)
run_workflow("configs/example.worldview.yml")

# An in-memory config works too; choose its relative-path base explicitly.
workflow = load_workflow(config_dict, config_dir="/data/project")
workflow.plan()
workflow.run()
context = workflow.records[0]["context"]  # {"const": {...}, "var": {...}}
scene_variables = context["var"]
constants = workflow.context["const"]
```

Run an ordered, sensor-neutral plugin recipe:

```bash
vhr workflow --config configs/example.worldview.yml --dry-run
vhr workflow --config configs/example.worldview.yml
vhr workflow --config configs/example.planet.yml
```

`--dry-run` prints existing output counts and actual pending work without running ordinary processing steps or cleanup. Scene-setting functions run to discover their records; their own import metadata saves may occur. Final processing JSON is only written by a real run. All step options live in the YAML. See [workflow configuration](../configuration/workflow-config.md) for paths, metadata mapping, plugin registration and resume behavior.

The former `vhr-worldview` command and flat configurations are not supported in version 3.

## Anchor images to a raster reference with joint registration

The `lidar_joint_coregistration` step in `configs/example.worldview.yml`
uses `joint_coregistration` with `core:scope: var`. It anchors each current
image to `const:reference_path` before the image-to-image joint registration
step. Set `core:run: false` to skip this per-scene anchoring. Workflow workers
can process the independent image/reference pairs concurrently.

Install the updated SpectralMatch checkout with its `elastix` extra, or use
VHRHarmonize's `elastix` extra once that SpectralMatch version is available.
The example environment includes the Elastix dependency.

`param:input_images` lists moving raster paths; `param:fixed_images` separately
lists reference raster paths. References are read-only and have no output
entries. `param:input_images_bands` and `param:fixed_images_bands` each default
to `[1]` and accept one-based band lists; multiple bands are averaged into
grayscale for matching. Output rasters retain all moving-image bands.

`param:moving_tie_method` defaults to `orb` for moving-to-moving pairs, and
`param:fixed_tie_method` defaults to `elastix` for moving-to-reference pairs.
Both accept `orb` or `elastix`. The fixed matcher is unused without references.
Each pair/window reads its selected bands once and reuses them for validity
and grayscale conversion. Matching methods leave the joint solver unchanged:
`param:local_model: none` gives a global correction, while `bilinear` or
`piecewise_affine` also applies a local correction.

For one combined reference and image-to-image solve, disable
`lidar_joint_coregistration` and uncomment `param:fixed_images` in the normal
`joint_coregistration` step. The default `partial_influence` coverage policy
uses reference ties only at valid reference pixels; moving-to-moving ties
still connect the network. Disconnected moving components may align without
a reference, while a solve with no usable ties errors.

## Py6S spectral bands and output validity

Py6S requires **exactly one** of these parameters (also checked after resolving
fields supplied through `variables`):

- `param:wavelengths_center_um`: band centers in micrometers, for example
  `[0.48, 0.55, 0.66]`. A single-band raster may pass a scalar.
- `param:wavelengths_bounds_um`: uniform-response `[lower, upper]` intervals
  in micrometers, for example `[[0.45, 0.51], [0.51, 0.58]]`. A single-band
  raster may pass `[0.45, 0.51]` directly.
- `param:wavelengths_response_csv`: a path to sampled spectral responses.

Each option must describe exactly as many bands as the raster contains, in
raster band order. The former `band_wavelengths_um` parameter and inline
response mappings are replaced by these explicit options.

### Spectral-response CSV convention

This project uses a UTF-8, comma-delimited wide table with a header and no
index or comment rows. This is a project convention, not a sensor-specific
interchange standard:

```csv
wavelength_nm,BAND_B,BAND_G
450,0.0,0.0
451,0.1,0.0
452,0.2,0.01
```

The first column must be named `wavelength_nm`. Subsequent columns have unique,
nonempty band labels and contain dimensionless relative responses; their order
must match the raster's band order. Wavelengths must be finite, strictly
increasing, and within 200–4000 nm. At least two wavelength rows are required.
Responses must be finite and nonnegative, with a positive response in every
band. Blank cells, duplicate wavelengths and mismatched column counts are
rejected. Arbitrary increasing wavelength spacing is supported.

Keep the original measured samples in the CSV. The plugin converts nm to
micrometers and linearly interpolates onto the 2.5 nm grid required by 6S only
at runtime; it peak-normalizes each interpolated band. When bounds do not land
on that grid, the small extensions are filled with zero. It does not average
the source samples into bins or rewrite the CSV.

The two example CSVs contain the original **1 nm** VNIR response samples,
without resampling, averaging, or additional normalization. Columns are
Coastal, Blue, Green, Yellow, Red, Red Edge, NIR1, NIR2, matching the example's
`const:band_order`. The original spreadsheet response values are already
normalized by their authors.

| CSV (under `configs/`) | Wavelengths | Original workbook |
| --- | --- | --- |
| `example.worldview2.response.csv` | 350–1100 nm; 751 rows | [DigitalGlobe WV2 response workbook](https://wp-cdn.apollomapping.com/web_assets/user_uploads/2011/10/14124816/QB02_WV01_WV02_Spectral_Response_Final_Correct_Jan292013.xlsx), `WV02 Relative Radiance Response` sheet, source columns G, C, D, H, E, I, F, J |
| `example.worldview3.response.csv` | 350–2450 nm; 2101 rows | [DigitalGlobe WV3 response workbook](https://apollomapping.com/image_downloads/WV3_VNIR_SWIR_response.xlsx), `Data VNIR SWIR 20131101` sheet, source columns C–J; VNIR responses beyond 1100 nm are zero |

Source workbook SHA256 checksums:

```text
WV2 ac376e045e37e9e50eb5dd3b01c74598b130345af0c6b1be559652283f83f5f9
WV3 ab4b095eb8cbaa2b02436770573126cdf72cca2eec46f11b21698f42a1b29eb9
```

`example.worldview.yml` selects the CSV using `IMAGE_1.satId`. Its `path:./...`
values resolve relative to the YAML, so the CSV files should stay beside it.
For another sensor, supply a CSV whose response columns match that raster.

### HPC staging

`wavelengths_response_csv` is declared as an input for dependency planning,
existence checks, protection, and HPC upload. `example.hpc.yml` maps
`var:wavelengths_response_csv` to the remote references folder. Preparation
uploads each required CSV once and remaps the selected paths when correction
runs remotely. If correction has already run locally and its outputs are
reused, the CSV is not needed for those loaded steps and is not uploaded for
them. A custom HPC recipe must likewise map the selected CSV or its containing
directory; declaring a file parameter alone does not choose a remote folder.

Py6S samples DEM pixels within the footprint's bounding window. Its output
keeps invalid input pixels and nonfinite correction results invalid; an output
NoData override does not filter input DN values. When an inherited NoData value
would collide with valid reflectance, the plugin selects NaN for float
outputs or an available integer sentinel. Integer outputs require an explicit
scale greater than one and fail on overflow instead of wrapping. Output band
scales convert stored values back to reflectance, and band tags record the
6S coefficients. `clip_reflectance: false` retains finite values outside [0, 1]
when the output dtype can represent them.

Coefficient precision remains limited by the installed 6S executable's printed
output; this plugin does not modify or recompile 6S.
