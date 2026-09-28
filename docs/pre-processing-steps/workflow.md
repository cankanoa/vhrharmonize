# Plugin processing

Top-level step declaration order and prefixed settings define execution and connections. Each unique step name selects an implementation with `plugin:`; pluginless steps can prepare constants and variables. Core provides `plugin: shared` for defaults. Processing plugins include `import_files`, `file_source`, `fetch_atmosphere`, `fetch_dem`, `atmospheric_correction`, `orthorectification`, `pansharpen`, `cloud_mask`, `alignment`, `seamline_metadata`, and the [individual SpectralMatch functions](../api/plugins-spectralmatch.md).

The WorldView recipe imports primary MUL rasters, maps IMD data to JSON, derives PAN paths, and explicitly connects correction, two orthorectification branches, pansharpening, cloud masking and alignment. The Planet example begins with already orthorectified reflectance. Both use the same runner and plugins.

Adapters receive one `params` mapping for all function arguments plus shared settings. File-feature declarations identify which arguments core manages as paths. The engine resolves YAML against workflow-wide `const` and per-scene `var` objects; functions receive those objects only through explicit parameters such as `param:settings: const:$` or `param:variables: var:$`. Py6S requires mapped viewing/solar geometry, acquisition day/month and explicit band wavelengths. When `dem_file_path` is supplied, Py6S samples ground elevation within the mapped geometry unless `ground_elevation_km` is explicitly set. `fetch_dem` can download an elevation raster from mapped bounds for subsequent steps. No sensor calibration defaults are embedded in the preprocessing code. FLAASH accepts mapped metadata or explicit native parameters.

Any function can establish or replace scenes by declaring `scene_records_return`. `import_files` is one such plugin; generated HPC recipes use `restore_scenes` for saved snapshots. Core refreshes scene state and replans the remaining steps after each replacement.

See [configuration](../configuration/workflow-config.md) for the dependency and resume rules, and [HPC](../cli/hpc.md) for remote execution.
