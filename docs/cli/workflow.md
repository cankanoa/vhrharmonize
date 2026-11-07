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

Press **Ctrl+C** to cancel a workflow while the progress display is active.
Cancellation interrupts waits for worker results, stops the workflow's local
process workers, and restores the terminal before returning to the shell.

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

## Atmospheric correction plugins

`plugin: py6s_atmospheric_correction` runs Py6S; `plugin: flaash_atmospheric_correction`
runs ENVI FLAASH. Each has its own parameters and module. The former
`plugin: atmospheric_correction` and `param:method` selector are replaced by
these explicit plugins. The public Python functions `run_py6s` and `run_flaash`
remain available from `vhrharmonize`.

The WorldView example enables the Py6S step and includes a separate commented
FLAASH step immediately below it. To switch algorithms, set the Py6S step's
`core:run` to `false` and uncomment the FLAASH step. Configure its ENVI engine
path and required output as needed. Both steps use the preceding enabled
step's image path.

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

## Atmosphere lookup warnings

Both `fetch_atmosphere` sources (`nasa_power` and `modis_gee`) return and save
`warnings`, `lookup_attempts`, and `defaults_used` in the configured atmosphere
JSON. With `param:log_to_console: true`, every warning is also printed with
`WARNING |` and the same message. JSON diagnostics are retained even when
console logging is disabled.

Warnings identify missing NASA POWER values at each point/date, missing MODIS
images/bands/valid regional values, alternate MODIS bands or collections, and
use of values from another UTC date. POWER also warns before trying each
alternate date. Its lookup order remains acquisition day, following day,
previous day (and further offsets when `search_days` is increased); another
date is tried only when the preceding date supplied no atmospheric values.
MODIS still selects the closest valid Terra/Aqua match within `hours_window`;
water-vapor and aerosol dates are reported separately.

Configure missing-value substitutions explicitly for either source:

```yaml
param:fallback_values:
  aot550: 0.2
  water_vapor: 2.5  # g/cm^2
  ozone_cm_atm: 0.3 # cm-atm
```

Each substitution generates a warning naming the parameter and value and is
listed in `defaults_used`. Omit this mapping to retain nulls, with explicit
unavailable-parameter warnings. Both sources expose `aot550`, `water_vapor`,
and `ozone_cm_atm` for the following correction step. MODIS also retains its
original raw/scaled fields; it does not retrieve ozone, so ozone needs a
configured default when using these Py6S examples. `lookup_attempts` preserves
the actual retrieved samples even when effective values contain defaults.

The examples apply defaults in `fetch_atmosphere`, rather than in later YAML
expressions, so the saved JSON describes the actual inputs passed to Py6S.
Request/authentication errors remain errors: the adapter saves their lookup
diagnostics with `status: error`, then re-raises the original exception. They
do not trigger an automatic provider switch or default substitution.

Pansharpening uses `core:tile_threads` for both Orthority worker pools and GDAL
warping/compression within each image. `core:concurrent_processing` still controls
concurrent images. Native Orthority overviews are disabled; request overviews
through `core:calculate_overviews` on the desired output step.

NASA POWER retrieval preserves each parameter from the closest date that supplies
it, searching the configured alternate days for parameters still missing before
using configured defaults. The saved JSON's `parameter_dates` records the source
date for each parameter (`null` for a default or missing value). The legacy
`date_used` field records the first date that supplied any value; consult
`parameter_dates` when values come from multiple dates.

Initial coverage (`create_footprints`) accepts `param:min_polygon_area: 1.0` to remove individual owned polygon pieces smaller than 1 square metre before buffering. Use YAML `null` (Python `None`) to disable this filter. The default applies when outline, ranking or buffering controls are used; ordinary raw footprint generation is unchanged.

The WorldView recipe caches per-image polygons in `polygonize_images`, beside the processed rasters. Initial coverage reads those GeoPackages with `footprint_source: vector_bounds`, `input_images: collect:polygonized_footprint_path`, and `vector_layer: footprints`. After alignment and matching, `create_seamline_footprints` regenerates polygons just before seamline postprocessing.

Coverage ranking can use `{stat:area_crs}` to calculate and rank each connected polygon's area in squared CRS units, after outline clipping and before overlap subtraction. This statistic is separate from `{var:field}` table columns. Saved vectors supply their original configured image ID, `image_path`, and attributes; the separate `input_vectors` parameter is no longer used by `create_footprints`.
