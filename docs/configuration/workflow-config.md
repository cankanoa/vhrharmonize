# Workflow configuration

Recipes map **unique step names** to settings. Each step may select one implementation with `plugin: registered_name`; ordinary steps run in declaration order. Names are arbitrary: `orthorectify_mul` and `orthorectify_pan` may both select `plugin: orthorectification`. Each step must be a mapping. Duplicate keys, unknown explicit plugins, plugin lists and unprefixed setting keys are errors.

A step without `plugin:` only evaluates `var:`, `const:` and `core:` settings. It cannot use function `param:` or `returned:` values. Use these steps for explicit setup and context updates. Every step, including setup and shared defaults, requires `core:run: true`; the default is disabled.

`plugin: shared` is provided by core. Enabled shared blocks supply global defaults before ordinary steps execute, regardless of their position. Multiple shared blocks merge in declaration order; later values override earlier ones. Use an ordinary pluginless step for changes that must happen at a specific point in processing. The name `shared` alone has no special behavior.

```yaml
shared:
  plugin: shared
  core:run: true
  param:epsg: 6635
  core:output_metadata_path: expr:var.relative_output_dir & '/processing.json'
  # core:delete_final_json_first: true
  const:output_dir: path:./output
  const:temp_dir: path:./temp
  core:cleanup_dirs: [const:temp_dir]
  const:reference_path: /data/reference.tif
  const:band_wavelengths_um: [0.4273, 0.4779, 0.5462]
  # core:run_from_existing: true

discover_inputs:
  plugin: import_files
  core:run: true
  param:search_glob: /data/images/*.tif
  var:relative_output_dir: >-
    path:expr:$replace(var.file_path, /[^\/]+$/, '') & '../output'
  var:mul: returned:file_path
  var:basename: expr:$replace($split(var.file_path, '/')[-1], /\.[^.]*$/, '')
  var:suffix: ""
  param:create_metadata_json:
    metadata:
      to_json: literal:expr:$replace(var.file_path, '.tif', '.json')
  var:solar_zenith: var:metadata.properties.solar_zenith

alignment:
  plugin: alignment
  core:run: true
  core:require_outputs: param:output_image_path
  param:moving_image_path: var:mul
  param:fixed_image_path: const:reference_path
  var:suffix: expr:var.suffix & '_aligned'
  var:aligned: expr:var.relative_output_dir & '/' & var.basename & var.suffix & '.tif'
  param:output_image_path: var:aligned
```

The WorldView example preserves setup and tuning comments, lists supported function parameters (including SpectralMatch), identifies defaults and comments out default settings. Every plugin explicitly shows `core:run`.

## Keys and values

A key and its value have independent roles. In `param:output_path: var:atmosphere_json`, the key selects the Python argument and the value reads a scene variable.

| Key | Meaning |
| --- | --- |
| `param:name` | Pass a resolved value as the Python function argument `name`. |
| `var:name` | Assign a per-scene variable. Dotted names update nested objects. |
| `const:name` | Assign a workflow-wide variable. Dotted names update nested objects. |
| `core:name` | Configure the runner. These settings control the runner and are never forwarded as function parameters. |

| Value | Meaning |
| --- | --- |
| Native YAML scalar, list or object | Pass the value, recursively resolving typed strings inside collections. |
| `var:name` | Read a scene variable, including dotted fields. Preserve its JSON type. |
| `const:name` | Read a workflow-wide variable, including dotted fields. Preserve its JSON type. |
| `expr:expression` | Evaluate JSONata with the combined context as the root. |
| `returned:field` | Read a field of the current function's return value. `returned:$` selects the whole result. |
| `collect:name` | Collect a scene variable across all current records. `collect:$` selects all scene `var` objects. |
| `literal:text` | Return text without interpreting a reserved prefix. |

Unquoted `param:output_path:` keys and `expr:var.suffix & '_aligned'` values are valid YAML. The prefix colon is not followed by whitespace. An expression containing YAML punctuation such as `: ` may need a quoted scalar or a block scalar:

```yaml
var:summary: >-
  expr:{"angle": var.solar_zenith, "quality": var.cloud_cover < 20 ? "clear" : "cloudy"}
```

Prefixes apply to complete strings; there is no embedded `${name}` interpolation. Ordinary strings, URLs and file paths are literals. `$` belongs to JSONata syntax inside `expr:`; it is not a second workflow reference system. Nested data keys such as `WV03`, `BAND_C`, GeoJSON's `type` or an ENVI task's `SENSOR_TYPE` keep their native names.

## One context, two scopes

The engine maintains a JSON context for resolving each plugin’s YAML settings:

```json
{"const": {"band_wavelengths_um": [0.4273, 0.4779]}, "var": {"basename": "scene_01"}}
```

`const` holds values shared across all scenes: calibration tables, wavelengths, common paths and aggregate results. `var` holds image paths, acquisition geometry, per-image calibration factors, returned values and naming state. Names can exist in both scopes without colliding. `const` describes scope, not immutability: a later enabled step can update an earlier constant.

Any initialized scene or aggregate step can assign both `var:` and `const:`. The prefixes select separate namespaces: `var:gain` and `const:gain` never collide. `collect:gain` gathers every scene's gain as a list, so `const:gains: collect:gain` stores that complete list as one shared value. Collection is available to scene functions too.

Within a scene step, `collect:` reads all scenes as they stood at the start of that step, so each invocation receives the same list.

**Ordinary `var:` assignments always evaluate separately for each scene**, regardless of `core:scope`. Both `var:name` and JSONata `var.name` refer to that scene's value. A literal list or object is stored intact in each scene; it is not automatically split. Dotted assignments retain the scene's other fields. `core:scope` determines how often the function runs, not how these assignments evaluate.

Aggregate function parameters have no selected scene. Direct `var:` references, including inside JSONata expressions or nested parameters, raise `ValueError`; use `collect:name` or `const:name` instead. To transform a collection with JSONata, explicitly collect it into a constant first and reference that constant in the expression.

An aggregate **`returned:` value assigned to `var:`** is distributed: the selected return field must be a list matching scene count or a dictionary keyed by exactly the current scene IDs. Incorrect lengths, missing/extra IDs and scalar returns raise `ValueError`. Each mapped item may itself be any JSON value. This also applies to returned selectors nested inside a scene assignment; subsequent expressions read that scene's mapped result.

Shared constants are arbitrary JSON values and have no scene-count requirement. Scene-independent constant expressions resolve once per enabled step (including with zero scenes), preserving the existing behavior of `const:scale: expr:const.scale + 1`. Constants derived from individual scene variables or function returns are allowed; scene invocations must agree on their shared value or core raises `ValueError`. Use `collect:` to retain differing values as a shared list, or `var:` to keep them per scene. Disabled steps make no assignments.

Scene-setting steps apply their `var:` mappings to the newly returned records. Their `const:` mappings may then collect or reference those discovered scene values. Before scenes exist, ordinary `var:` reads/writes and `collect:` remain invalid; `plugin: shared` runs before discovery and initializes constants only.

For a batch step, capture inputs before replacing the current paths:

```yaml
match:
  plugin: global_regression
  core:run: true
  param:input_images: collect:current_image_paths
  var:current_image_paths: expr:const.temp_dir & '/matched/' & var.basename & '.tif'
  param:output_images: collect:current_image_paths
```

The first collection reads old paths and the second reads the replacements. Later scene steps use `param:input_path: var:current_image_paths` to receive their single path. Returned assignments use the same count/ID validation after execution; known output expressions allow planning and HPC staging before execution.

Expressions use [JSONata](https://docs.jsonata.org/overview.html), evaluated by [jsonata-python](https://github.com/rayokota/jsonata-python). JSONata sees both scopes: `expr:var.relative_output_dir & '/' & var.basename & '.tif'`. Use `&` for strings and JSONata functions such as `$map`, `$lookup`, `$replace` and `$substring`. Python comprehensions and function calls are not supported.

Assignments resolve in declaration order and may update an existing value, such as `var:suffix: expr:var.suffix & '_aligned'`. Function parameters see preceding assignments. `returned:` assignments and their dependents take effect after the function returns. A function parameter cannot depend on that same invocation's result; place its input reference before the return assignment. An explicitly passed scope snapshot contains the values at its parameter’s position in the settings. Missing fields or undefined expression results raise `ValueError`; JSON null remains `None`. JSONata's `??` supplies a default for a missing field.

### Nested objects

Tables stay ordinary nested YAML; only the setting itself needs a prefix:

```yaml
shared:
  plugin: shared
  core:run: true
  const:calibration:
    WV02:
      BAND_C: [0.977, -5.552]
    WV03:
      BAND_C: [0.905, -8.604]
```

`const:calibration.WV03.BAND_C` returns the entire array. `const:calibration.WV03.BAND_C.0` returns its first value. JSONata uses brackets: `expr:const.calibration.WV03.BAND_C[0]`. For a sensor chosen per scene, use `expr:$lookup(const.calibration, var.sensor_id).BAND_C[0]`. A key such as `const:calibration.WV03.BAND_C` replaces that nested field, retaining its siblings; dotted assignment segments must be identifiers. Arbitrary object keys can be supplied in a nested mapping and read with JSONata `$lookup`.

### Function defaults and returned values

Function argument precedence is **explicit plugin `param:` > shared `param:` > function default**. Runner `core:` controls are never passed as function arguments; use an explicit `param:` even when the name is the same. Defining a `var:` or `const:` field never supplies a same-named argument, including file paths. Connect it explicitly:

```yaml
param:solar_zenith: var:solar_zenith
param:band_wavelengths_um: const:band_wavelengths_um
```

Shared `param:` settings are offered to every processing adapter; `FunctionPlugin` forwards only accepted arguments, applying adapter aliases when necessary. Unsupported shared settings are ignored, while unsupported explicit plugin parameters are errors. A missing required argument raises an error even if a matching field exists in the context.

Functions and processing adapters receive no automatic context argument. Names such as `context` and `variables` have no special behavior. Whole scopes can be passed as ordinary arguments:

```yaml
param:variables: var:$   # All scene variables.
param:settings: const:$  # All workflow constants.
param:context: expr:$    # Both scopes, retaining their const/var namespaces.
```

Choose parameter names that match the function. `"var:"` and `"const:"` also select whole scopes but must be quoted in YAML. A merged object can be requested with `param:variables: expr:$merge([const, var])`. All these explicit references track runtime dependencies automatically; no adapter dependency declaration is needed. Passed objects are snapshots. Mutating them does not publish changes: return data and assign it explicitly in YAML.

Use explicit [context operations](#explicit-context-files) to persist and load selected returned values. Existing output files alone do not restore metadata. Value-only plugins run when needed. Unavailable values from unused branches are omitted from final context. Python callers can inspect `workflow.records[i]["context"]` and `workflow.context["const"]` after execution.

## Importing files and scene identity

`import_files` is an ordinary `FunctionPlugin` wrapping `import_files(...)`. It declares `var_records_return = "scenes"`; any plugin can declare a returned field containing a list of plain dictionaries to establish scenes. The adapter selects whether to replace or merge them. `import_files` uses merge mode. No plugin is required to be first. Before any plugin establishes scenes, enabled functions run once using ordinary parameters and constants. Reading or assigning `var` raises `ValueError`, including through JSONata or `collect:`. The context contains only `const` at this point, so `expr:$` still works. Once scenes are established, scene functions run per record and aggregate functions run once. An initially empty list establishes zero scenes. A later empty import preserves existing scenes. See [adding plugins](../getting-started/adding-plugins.md#creating-or-replacing-scenes) for the contract and optional declarations.

The import function returns `{"scenes": [...]}`. Each scene contains `scene_id`, `file_path`, `source_paths` and the fields named in `create_metadata_json`. Each dictionary becomes a scene's `var` object directly. The importer creates no directory roots or constants; define those with YAML assignments. For example:

```json
{
  "scenes": [{
    "scene_id": "/data/scene.TIF",
    "file_path": "/data/scene.TIF",
    "source_paths": ["/data/scene.TIF", "/data/scene.RPB", "/data/scene.IMD"],
    "rpc": ["/data/scene.RPB"],
    "metadata": {"IMAGE_1": {"cloudCover": 0.1}}
  }]
}
```

There are no automatic `filename`, `basename`, `directory` or `suffix` fields. Derive names from `var.file_path` and initialize naming state in YAML, as above. Within scene assignments, `returned:` selects from that individual dictionary: `returned:metadata.IMAGE_1` and `var:metadata.IMAGE_1` access the same imported data. Function parameters and `const:` assignments use the whole invocation's context/result.

`param:search_glob` accepts one pattern or a list, including recursive, brace and extended globs. The WorldView example selects primary MUL files and derives PAN paths from them. The ordinary function argument `param:scene_id` selects the merge identity. It defaults to the per-file reference `var:file_path`; YAML writes that default as `literal:var:file_path`. Supply a literal ID, metadata reference, or JSONata expression:

```yaml
# Default: identify by the complete file path.
# param:scene_id: literal:var:file_path
# Identify by a parsed metadata field instead:
param:scene_id: literal:var:metadata.properties.id
# Or derive an ID from the filename:
# param:scene_id: literal:expr:$replace($split(var.file_path, '/')[-1], /\.[^.]*$/, '')
```

The importer evaluates this after `create_metadata_json`, before `where` and the YAML scene mappings, and returns the result as `scene_id`. Per-file expressions therefore use imported fields, not later aliases such as `var.basename`. The ID must be a nonempty string; use `$string(...)` in JSONata for numeric metadata IDs. IDs must be unique within one import call. Separate imports merge records with the same ID even when their file paths differ; different IDs create separate scenes even for the same file. Configure identity through `param:scene_id`; assigning a later `var:scene_id` does not change core's record identity. Output paths must remain unique. In Python, use `import_files(pattern, scene_id="var:metadata.properties.id", ...)` without the YAML `literal:` delay.

`param:create_metadata_json` maps arbitrary field names to rules. Each rule must contain exactly one of:

- `path`: a glob or list of globs. Return a sorted, deduplicated list of existing files, or `[]` when nothing matches.
- `to_json`: one metadata file. Decode JSON/GeoJSON, IMD, XML or YAML into a JSON-compatible object. Missing files or invalid metadata raise an error. No JSON file is written.

```yaml
param:create_metadata_json:
  rpc:
    path: literal:expr:$replace(var.file_path, '.TIF', '.{RPB,rpb}')
  metadata:
    to_json: literal:expr:$replace(var.file_path, '.TIF', '.IMD')
param:where: literal:expr:100 * var.metadata.IMAGE_1.cloudCover <= 75
```

Rules resolve relative paths from the matched file's parent; absolute paths and `~` are also supported. Both rule types add discovered inputs to `source_paths`, which core uses to protect originals when `core:protect_source_files` is enabled. Field names `scene_id`, `file_path` and `source_paths` cannot be overwritten by rules.

`literal:` delays the expression until the importer evaluates it for each discovered file. Rules run in declaration order and may read `var.file_path` and fields produced by earlier rules. The filter runs after all rules, before YAML scene mappings; therefore it reads `var.metadata`, not a later YAML alias. Plain relative paths need no expression or `literal:` prefix. Block mappings keep expressions containing commas unquoted; YAML flow mappings such as `{path: '*.RPB'}` also work.

The importer receives no workflow constants object. Use ordinary `const:` or `expr:` arguments to resolve specific values before the call, for example `param:search_glob: const:input_pattern`. Per-file expressions see the imported scene fields. Parsing preserves structure; YAML defines sensor mappings and calibration. A later step can consume `var:rpc` directly through `core:requires` or `file_source`'s `param:companions` argument. Shared function parameters follow the usual precedence.

### Importing more files later

Use unique step names with `plugin: import_files` wherever another batch should join the workflow:

```yaml
initial_images:
  plugin: import_files
  core:run: true
  param:search_glob: /data/initial/*.tif
  param:scene_id: literal:expr:$replace($split(var.file_path, '/')[-1], /\.[^.]*$/, '')
  var:initial_image: returned:file_path
  var:relative_output_dir: path:/data/products

# Processing steps for the initial images can go here.

additional_images:
  plugin: import_files
  core:run: true
  param:search_glob: /data/additional/*.tif
  param:scene_id: literal:expr:$replace($split(var.file_path, '/')[-1], /\.[^.]*$/, '')
  var:additional_image: returned:file_path
  var:relative_output_dir: path:/data/products
```

Imports build on the existing context. New IDs add new scenes; matching IDs merge into the existing scene without duplicating it. The example above matches filenames across two folders. Existing values win, including lists and nulls; nested objects gain missing fields. Per-import `var:` mappings initialize new scenes and missing fields without resetting existing values such as a processed path or suffix. Previously completed or cached function returns remain available. Source-protection path lists accumulate all imported inputs. An empty import does not clear the scene list. When different paths share an ID, the original `file_path` stays intact. Use a new alias such as `var:additional_image: returned:file_path` to retain the newly imported path.

Directory variables follow the same assignment rules as other metadata. Per-import `var:` assignments fill new scenes and missing fields; existing scene values survive reimport. Define shared roots once with `const:`; an explicit later `const:` assignment changes the shared value.

All subsequent scene steps process the combined collection. New scenes do not run earlier steps, so downstream parameters must be available on both existing and new scenes. Imports that depend on unfinished work execute at their position and rebuild the remaining plan. HPC preparation can stage multiple imports when their files and preceding required results are already available; unresolved later discovery remains a preparation error. `file_source` only copies files and its companions; it does not establish or replace scenes.

### Directory roots

Directory roots are ordinary YAML values. For example:

```yaml
shared:
  plugin: shared
  core:run: true
  const:output_dir: path:./output
  const:temp_dir: path:./temp
  core:cleanup_dirs: [const:temp_dir]

discover_inputs:
  plugin: import_files
  core:run: true
  param:search_glob: /data/images/**/*.tif
  var:relative_output_dir: >-
    path:expr:$replace(var.file_path, /[^\/]+$/, '') & '../../processed'
```

`path:./output` and `path:./temp` resolve from the YAML directory. The per-scene expression derives an absolute output root from the original input's parent; `path:` normalizes `..` components. Use `path:sys` to explicitly allocate a system temporary directory instead. A stable path is preferable when saved context must survive between runs. No directory is allocated by `import_files` itself.

The names have no intrinsic meaning: use any `const:` or `var:` field in downstream path expressions. The WorldView example puts per-scene products under `var.relative_output_dir`, intermediates under `const.temp_dir`, and the final mosaic under `const.output_dir`. Output directories need no core declaration. HPC `path_mappings` can reference these values directly.

`core:cleanup_dirs` is a shared control, default `[]`. It accepts one `const:name` or `var:name` reference, or a list of these references (including nested fields). Each selected value can contain one root or a list of roots. It makes roots eligible for cleanup; it does not create values or request processing. Missing or not-yet-available fields are ignored until populated. Existing plugin `temporary_directory_context_paths` declarations remain supported and are combined with these explicit selections.

Cleanup still requires the function's `output_temporary_cleanup_paths` declaration and the existing cleanup controls. It removes only owned, consumed intermediate outputs, retaining requested products and protected inputs. `core:require_outputs` selects deliverables independently of directories.

**Migration:** remove `param:temp_dir`, `param:output_dir`, `param:temp_dir_scope`, and `param:output_dir_scope` from `import_files` steps. Replace them with explicit YAML assignments, update path expressions/HPC mappings, and select cleanup roots with `core:cleanup_dirs` if needed. The removed Python arguments are rejected rather than silently ignored. Processing functions that independently accept `param:temp_dir` still support it.

Scene-setting functions may run during planning, dry-run and HPC preparation without an opt-in. Ordinary processing functions are not run by dry-run. The planner stops at scene discovery whose preceding required work has not finished, reporting later enabled steps as `pending`; execution resumes planning after the scene update. HPC requires scenes and file paths to be known during preparation. HPC recipes preserve ordinary imports and explicit context loads; no scene snapshot is embedded in the generated YAML.

### Final metadata JSON

The final JSON destination is an explicit shared YAML control. Its value may be an absolute path, any `var:`/`const:` reference, or a JSONata expression:

```yaml
shared:
  plugin: shared
  core:run: true
  core:output_metadata_path: var:report_path # path | null; default: null (disabled).
  # core:delete_final_json_first: true      # true | false; default: true.
```

Define `var:report_path` in an enabled plugin, for example `expr:var.relative_output_dir & '/' & var.basename & '.json'`. A constant destination such as `/data/results/processing.json` combines scenes in one file. Relative destinations resolve from the YAML directory.

Each completion of the last needed step appends its combined `{const, var}` context to a JSON array. A scene step writes as each scene finishes, including when its output is reused from cache. A final aggregate step writes one entry per scene using the aggregate's final constants; a workflow without scenes writes its const context. Discovery-only runs export their imported scenes. Different resolved paths produce separate files naturally.

`delete_final_json_first: true` replaces old contents on the first successful write to each destination during that workflow run. Core remembers canonical absolute paths, so subsequent writes to that path append without clearing earlier completions. With `false`, previous entries are retained; an existing JSON object becomes the array's first entry. Replacement is atomic. Dry-run does not clear or write this final JSON.

Imported metadata stays in memory unless selected by an explicit context save or included in this final export. Implicit output checkpoints are not used.

## File paths and names

Selected file arguments accept a path or flat list of paths. A list registers each file individually for dependencies, collisions, reuse, validation and transfers. Existing directories are never considered complete cached products by core; the invoked function manages their contents and any per-file resume policy.

Destinations such as `param:output_path` and `param:output_image_path` are ordinary function arguments. Adapters independently select function parameter names for each core file feature. The engine does not infer behavior from argument names. All selections default to empty; an argument absent from every selection is an ordinary value passed to the function.

### Plugin file features

These are Python adapter attributes, not YAML parameters. Each contains a set or sequence of function parameter names. Features are independent: selecting a path for one job never opts it into another. Built-ins explicitly select the features they need.

| Feature declaration | Core behavior | YAML control |
| --- | --- | --- |
| `output_path_resolution_paths` | Normalize selected destinations and expand `~`. | Always for selected paths. |
| `input_dependency_paths` | Link selected inputs to preceding file producers. | Always; `core:requires` adds extra dependencies. |
| `output_dependency_paths` | Register selected destinations as file producers. | Always. |
| `output_target_paths` | Internal target selection populated by core from YAML; adapters leave it empty. | `core:require_outputs: false` by default; select `param:name` or a list of selectors. |
| `input_existence_check_paths` | Require selected inputs to exist before invocation. | Always. |
| `input_protection_paths` | Protect selected inputs/references from overwrite and cleanup, including aliases. | Always for selected inputs; imported source protection is separately controlled by `shared.core:protect_source_files: true`. |
| `output_parent_creation_paths` | Create parents of selected destinations before invocation. | Always. |
| `output_collision_check_paths` | Reject duplicate destinations when either declaration selects collision checking. | Always. |
| `output_reuse_paths` | Allow skipping a function when its requested products already exist. Every requested product must be selected for reuse. | Plugin `core:reuse`; default `shared.core:run_from_existing: true`. |
| `output_validation_paths` | Validate selected products during reuse and after invocation. Reuse always checks existence; format checks apply only to this selection. | Plugin/shared `core:check_validity: true`; `shared.core:validity_check_grid_size: 2048`. |
| `output_invalid_removal_paths` | Inspect and remove corrupt selected products before regeneration. No missing-file error or post-call check is implied. | Runs before processing, independently of reuse/validation selections. |
| `output_overview_calculation_paths` | Build overviews on selected TIFF outputs after processing. | Plugin `core:calculate_overviews: false`; requires `shared.param:window_scales` when enabled. |
| `output_temporary_cleanup_paths` | Remove selected regular files inside declared temporary roots, with associated sidecars, after required consumers succeed. Only applies under populated temporary roots; explicitly required products are retained. | Shared `core:delete_temp_steps_proactively: true` / `core:delete_temp_dir: false`. |
| `input_hpc_staging_paths` | Rewrite and upload selected inputs needed by processing steps. | HPC remote reference directory and upload settings. |
| `output_hpc_staging_paths` | Rewrite selected destinations and upload selected reusable products needed by the remote run. | HPC remote output/temp directories and upload settings. |
| `output_hpc_download_paths` | Rewrite selected destinations and include requested products in the download map. Does not enable upload. | `core:require_outputs`, HPC `path_mappings`, and download-conflict setting. |

For example, a plugin can validate an image and a JSON report while building overviews only for the image:

```python
output_validation_paths = {"output_image", "output_report"}
output_overview_calculation_paths = {"output_image"}
```

These selections do not imply path normalization, directory creation, reuse or cleanup; declare those separately when needed. Built-ins sometimes assign the same immutable `frozenset` to several features as an explicit convenience. Reassigning one feature leaves the others unchanged. The old generic declarations and `manages_reuse`/`manages_overviews` flags are rejected. A nested pipeline leaves the relevant core selections empty and performs those operations itself.

Selected paths remain subject to structural constraints: output destinations must be single paths known during planning, and the engine cannot overwrite protected sources. Paths not selected for normalization are passed unchanged; use absolute paths for other core file operations unless working-directory-relative paths are intentional. Cleanup also removes associated raster sidecars. Input protection keeps intermediate files alive until their required readers finish, even when those readers omit file-dependency tracking. `core:requires` remains an explicit YAML shorthand for extra file dependencies, including resolution, existence checks, protection and HPC staging.

Relative discovery globs resolve from the Python working directory. The importer resolves metadata files and companion globs from each discovered file's directory. Core passes input parameters unchanged. Selected relative output parameters resolve from the YAML directory. Registered output roots do not affect resolution. `const:` alone does not give a string path semantics: its use by a selected path-resolution parameter does. Use absolute common file paths when scene directories differ. `~` expands on the execution host. Output destinations must resolve during planning; HPC also requires inputs selected for staging to resolve then. Returned scalar/object values can remain deferred until execution.

Files embedded in compound options need an explicit dependency. For example, in an aggregate step:

```yaml
const:mask: /data/mask.gpkg
core:requires: const:mask
param:global_regression_vector_mask: [include, const:mask, image]
```

Filename accumulation is explicit scene state: `var:suffix: expr:var.suffix & '_ortho'`. The output expression combines `var.basename`, `var.suffix` and an extension. The WorldView recipe accumulates only the MUL suffix; its PAN output uses a fixed `_pan_ortho.tif` ending. Output paths are explicit function arguments. Core and plugins do not append naming suffixes; YAML expressions define the complete paths before planning.

## Controls, reuse and cleanup

All steps default to disabled. `core:run: false` does nothing: no adapter loading, value resolution, assignments, cache restoration, processing, transfers or cleanup. Downstream inputs must explicitly reference an enabled step or an imported file. A plugin-only CLI command keeps enabled upstream names and reuses their files, but does not compute upstream processing steps.

Use `core:require_outputs` to choose deliverables independently of file locations:

```yaml
alignment:
  plugin: alignment
  core:run: true
  core:require_outputs: param:output_image_path
  param:moving_image_path: var:current_image_paths
  param:fixed_image_path: const:reference_path
  var:current_image_paths: path:expr:var.destination & '/aligned.tif'
  param:output_image_path: var:current_image_paths
```

The default `false` leaves the step available as a dependency. A `param:name` selector or list such as `[param:output_images, param:save_adjustments]` requests those output parameters. Each selected parameter must be a declared output and resolve during planning to a nonempty path or nonempty flat list of paths. Unknown parameters, inputs, nulls, objects and nested lists are errors. `true` requests the step and all supplied output parameters, including functions with no file outputs. Reuse still applies: requesting a result does not force recomputation. `core:requires` remains separate and declares extra **input dependencies**.

### Explicit path values

Prefix a value with `path:` to normalize it before it reaches a function:

| Value | Resolution |
|---|---|
| `path:/data/image.tif` | Absolute path. Native Windows drive and UNC paths follow the host's normal path rules. |
| `path:~/data/image.tif` | Expand the current user's home directory. |
| `path:./results`, `path:../results`, `path:results` | Relative to the workflow YAML directory (the API's `config_dir` for in-memory configurations). |
| `path:var:images`, `path:const:reference`, `path:collect:images` | Resolve the reference, then normalize the resulting path or flat list. |
| `path:expr:const.folder & '/result.tif'` | Evaluate the expression, then normalize. |
| `path:sys` | Create a unique system temporary directory, stable for that assignment during planning, execution and staging. |

Declare a shared temporary path once, for example `const:scratch: path:sys`, then reference it with `const:scratch`. Scene assignments can create separate paths per scene. `path:` does not mark a path as an input, deliverable, upload or cleanup candidate. `var:` and `const:` alone simply read values. There are no `output:` or `temp:` resolution prefixes. Plugin-selected legacy output normalization remains available, using the YAML directory; other unprefixed values keep their existing function-specific semantics. The importer resolves metadata and companion rules relative to each input; `path:` values arrive already absolute. Write `path:sys` without a space after the colon, or quote the entire scalar.

Migrating older recipes: rename `core:process_for_paths` to `core:require_outputs`, keeping its value unchanged. Explicitly request desired products, including the final processing step. Moving a file outside temporary storage no longer requests it. Use `path:` or explicit expressions for destinations that previously depended on a registered output root.

Per-step controls are `core:skip_plugin_call`, `core:run`, `core:scope`, `core:reuse`, `core:check_validity`, `core:calculate_overviews`, `core:require_outputs`, `core:requires` and `core:processing_direction`. Scope defaults to the adapter's declaration. Seamline metadata and SpectralMatch declare `aggregate`; their names receive no special execution branch in the core.

Shared runner controls are:

| Name after `core:` | Default |
| --- | --- |
| `protect_source_files` | `true` |
| `output_metadata_path` | `null` (disabled) |
| `delete_final_json_first` | `true` |
| `run_from_existing` | `true` |
| `check_validity` | `true` |
| `validity_check_grid_size` | `2048` |
| `cleanup_dirs` | `[]` (explicit `const:`/`var:` cleanup-root references) |
| `delete_temp_dir` | `false` |
| `delete_temp_steps_proactively` | `true` |
| `log_to_console` | `true` |
| `show_progress` | `true` |
| `report_progress` | `false` |
| `save_statistics_path` | `statistics.jsonl` |
| `load_statistics_path` | `statistics.jsonl` |
| `processing_direction` | `vertical` |
| `concurrent_processing` | `1` |
| `concurrent_processing_backend` | `process_pool` |
| `dask_scheduler` | `null`; `[file, path]` or `[address, URL]` |

The engine implements these controls; plugin functions do **not** need to accept them all. There is no mandatory set of raster parameters such as `custom_nodata_value`, `output_dtype`, `epsg` or `window_scales`. Functions accept the settings they support. Raster options use SpectralMatch names such as `custom_nodata_value`, `output_dtype` and `window_scales`. Native utilities that use `custom_output_dtype` retain that name; set it directly when overriding their dtype. The engine additionally uses shared `param:window_scales` when `core:calculate_overviews` requests raster overviews. Per-step `core:reuse` and `core:check_validity` override shared reuse/validation. Raster validation checks readability/TIFF bounds; JSON validation checks syntax.

The planner follows file and runtime-variable dependencies backward from deliverables. A cached cloudmasked image can feed alignment even when correction, orthorectification and pansharpen temporary files are absent. Static suffixes/constants do not force upstream recomputation. Only explicit `core:require_outputs` targets request processing products, regardless of directory. Terminal processing steps are not requested automatically. `core:require_outputs: true` requests the step and all supplied declared outputs; use this for steps with side effects but no output files. Scene discovery still runs to establish the plan, and explicitly selecting a plugin through `run_plugin` requests its enabled steps. `loaded` counts reusable outputs, `processing` counts required work and `unused` counts bypassed nodes; a loaded node can also be unused. `pending: 1` marks a step awaiting scenes from a runtime scene-setting function, so its invocation count is not yet known.

With `core:log_to_console: true`, core announces workflow startup, discovery, graph construction, planning, execution and enabled final cleanup using `[core:<stage>] Start`. The first line, `[core:workflow] Start`, appears as soon as core settings are resolved, before input discovery or processing. Cached calls to `plan()` do not repeat the announcement; a new graph and plan after scene discovery do. Early stages print before the dashboard starts, while execution-time messages enter the progress API's recent messages and terminal display. Setting `core:log_to_console: false` suppresses these core messages.

With the dashboard disabled, the scheduler also logs `[scene_id core:step_name] Start current/processing/total` immediately before each function call or worker submission. `current` is the dispatch number for this step, `processing` is the number of calls needed in this run, and `total` includes cached and unused calls. Aggregate steps count function calls, so a single call handling many images has a total of one. These are start counters; concurrent calls can finish in a different order. Single-image operation logs end with `Completed` without a separate counter.

Shared `core:save_statistics_path` appends raw core timings across executions, and
`core:load_statistics_path` uses successful past timings to seed console and dashboard
ETAs. Both default to `statistics.jsonl` beside this YAML. Set either to `null` to
disable it independently. A missing load file is allowed; malformed existing history
raises a file/line error. Choose filenames yourself for different input sizes/types.
The old `statistics_path` setting is renamed to `save_statistics_path`.
Run `vhr statistics --input-path statistics.jsonl --output-path summary.json` to
generate a separate report. Recording enables the reporting backend independently
of the display and console logging. See
[Statistics](../api/statistics.md) for the format and event API.

### Live progress dashboard

The dashboard is enabled by default. Its setting can be overridden in a shared block:

```yaml
shared:
  plugin: shared
  core:run: true
  core:show_progress: true
```

Set `core:show_progress: false` to disable the `prompt_toolkit` frontend. `core:report_progress: true` collects progress and writes snapshots without a display; `show_progress: true` implies reporting. With both controls disabled and no Python reporting callback, explicit snapshot path or statistics file, ordinary console logging is used. Logging controls remain independent: `core:log_to_console` controls core messages, while plugin messages follow `param:log_to_console`.

In an interactive terminal, logs stay in the original terminal session and scrollback. `prompt_toolkit` redraws one live box below new output, titled **VHRHarmonize Workflow Progress**, with a single divider between the workflow table and active operations. Both tables place Progress in the rightmost column. The display uses no alternate screen or mouse capture: use the terminal’s normal scrolling, selection and copy controls. Scrolling back can move the progress box out of view. Rows use the **configured step name**, including repeated uses of the same plugin. Core-generated overviews have a separate `<step name> / overviews` row. `(no cb)` beside a name means detailed worker callbacks are unavailable; core still tracks that step's starts, completions and elapsed time. On smaller terminals, active steps take priority and an omitted-row count appears; the total always includes every row. Very short terminals show a compact total and first active operation. Widen the terminal to see all columns. On exit, pending log lines are flushed and one final summary remains without repeating messages already printed. Redirected output (including Slurm log files), noninteractive input and `TERM=dumb` receive a single final summary without terminal animation. Redirected output has no escape sequences, and ASCII-only streams use simple borders.

The columns are:

- **Unused, Loaded, Done, Run, All:** five separate columns, each displaying `count(percentage of All)`. Unused means scenes the planner does not need for this step. Loaded is the existing cached/reused-output count (`reused` in the API). Done means successful completions **during this run**, excluding previously completed outputs. Run is all work selected for this run, including Done; it does not decrease as work finishes. All is the full planned scene count, including loaded and unused scenes. These columns overlap: Done is a subset of Run, and loaded outputs can also be unused. Aggregate calls represent the current scene batch; their Done count advances by the batch size only after successful validation. Overview rows count raster files. Steps after runtime scene discovery show `TBD` in ETA until an estimate is available.
- **Progress:** segments run left to right as unused (blue), loaded (purple), done (green), and remaining work from Run (light gray). The four column titles use the same fixed colors, chosen for visibility on both light and dark terminals. The full bar represents All. Cached scenes that are also unused appear once, in the Loaded segment; the blue segment covers other unused scenes. Forced reprocessing does not count old outputs as loaded. Percentages appear only in the count columns; the bars have no percentage labels or bottom legend.
- **Active:** function calls currently executing, excluding queued submissions. The top `total` row sums the counts and active calls below it.
- **ETA:** approximate remaining work based on completed-call durations. Per-step estimates use observed concurrency; the total combines remaining work and divides by available workflow workers (Dask worker threads for Dask). Unknown estimates show `TBD`. Dependencies, unequal scene sizes and changing concurrency can change the estimate. Cached scenes do not enter timing samples.

The active-operation panel uses one row per operation, with **Step, ID, Status, Elapsed, ETA, Progress** columns. Elapsed is shown in seconds; unknown ETAs show `TBD`. Status is the callback's status when provided, otherwise `starting` or `working`. The bar uses tqdm-style tile, band, point or image progress; backends without measurable progress use an indeterminate bar. Static/HPC displays include every active operation; small live terminals show a count of additional operations when space runs out. Worker processes send events to core rather than creating their own terminal displays. See [function progress callbacks](../api/io-logging.md) for the Python interface.

When running a YAML file with reporting enabled, core writes `<workflow.yml>.progress.json` beside that file. The snapshot is replaced atomically at most once per second, with initial and final updates on success or failure. `vhr hpc-status --config <staged.hpc.yml>` fetches the snapshot through `get_slurm_progress()` and gives it to the same `prompt_toolkit` frontend for a static display. This works while Slurm output is redirected. Snapshots from a different Slurm job ID are ignored, and the displayed timestamp identifies the last update. The staged HPC YAML retains the fetched data under `workflow_progress`. Apps can receive the same versioned data through a Python callback, `Workflow.get_progress()`, the JSON file, or the HPC accessor; see the [progress API](../api/progress.md).

Proactive cleanup applies to **`output_temporary_cleanup_paths` regular-file outputs under selected temporary roots** (`core:cleanup_dirs` or plugin declarations). It waits for all required consumers to succeed and protects imported inputs, reference files and their aliases. Explicitly required products and terminal temporary results are retained. `delete_temp_dir` performs this cleanup at the end and removes emptied subdirectories; it does not recursively erase the root. Earlier temporary outputs are retained across runtime scene replacements because future consumers were not yet known. Undeclared plugin caches and directory outputs are not automatically removed. Private diagnostic files remain the function’s responsibility. The individual SpectralMatch steps expose their intermediate raster paths to core cleanup.

### Horizontal and vertical execution

`core:processing_direction` accepts `horizontal` or `vertical` (default) in shared settings and on individual steps. A step setting overrides the shared setting. Horizontal execution completes a step across all required scenes before moving on. Consecutive vertical scene steps let each scene advance as its own preceding work finishes; with one worker, one scene completes that section before the next starts. With multiple workers, a fast scene can advance while another scene is still upstream.

```yaml
shared:
  plugin: shared
  core:run: true
  core:processing_direction: vertical
  core:concurrent_processing: 4

# A particular step can request a full synchronization point:
check_images:
  core:run: true
  core:processing_direction: horizontal
```

Core synchronizes automatically at `collect:` reads, aggregate calls, scene-setting plugins, and steps returning or deriving shared constants that are unavailable until execution. Scene invocations must agree on shared values before later steps consume them. Explicit horizontal steps also synchronize. Disabled steps do not interrupt a vertical section. Within a vertical section, declared file and value dependencies are still enforced, including dependencies on another scene. Static filename/constant assignments do not force an extra synchronization point.

`core:scope` remains independent: it selects one function call per scene or one aggregate call. `processing_direction` only changes scheduling. The process-pool worker limit applies across the whole vertical section; Dask uses cluster capacity and gives downstream tasks higher priority. Cached outputs still skip unnecessary work, and temporary files are removed only after their required consumers finish and explicitly loaded values have been published. HPC staging preserves the selected direction.

## Aggregates and concurrency

An aggregate uses `param:input_images: collect:current_image_paths` or `param:metadata_records: collect:footprint_metadata`. Collection reads scene variables at that position in the workflow. Aggregate `const:` assignments become available to subsequent scene and aggregate steps. Repeated invocations use separate unique step names selecting the same `plugin:`. Counts, explicit context operations and execution order identify those names, such as `orthorectify_mul` and `orthorectify_pan`.

Scene steps can process records concurrently. `core:concurrent_processing` accepts a positive integer or `num_cpu`. Dask requires `core:concurrent_processing: 1` plus a scheduler address/file when work must run. Horizontal processing completes each step before the next begins; vertical processing advances scenes independently between synchronization points. Aggregate functions manage their internal parallelism. SpectralMatch inherits supported shared `param:` settings; no adapter translates core worker counts or scheduler settings into native parameters.

See [Adding plugins](../getting-started/adding-plugins.md) for registration and [HPC execution](../cli/hpc.md) for staging.


## Individual SpectralMatch functions

`setup_spectralmatch` in the WorldView example is an ordinary pluginless aggregate step. It defines shared polygon identifiers. Each image keeps its path and filename suffix in `var`; batch functions receive paths through `collect:`. Ordinary settings are direct `param:` entries under `shared`, inherited automatically by functions that accept them. The enabled example runs matching, footprint generation, Markov seamlines and masking, with `merge_rasters` last. Alternative stages remain fully commented with their parameter catalogs. Disabled steps leave the current image list and suffix unchanged. Copy a step under another name to invoke the same plugin again.

```yaml
shared:
  plugin: shared
  core:run: true
  param:output_dtype: float32
  param:window_scales: [2, 4, 8, 16, 32]

setup_spectralmatch:
  core:run: true
  core:scope: aggregate
  const:spectralmatch.images: collect:aligned

match_global:
  plugin: global_regression
  core:run: true
  core:calculate_overviews: true
  param:input_images: const:spectralmatch.images
  param:output_images: [/work/a_global.tif, /work/b_global.tif]
  const:spectralmatch.images: returned:$
```

The adapters import algorithms from the installed `spectralmatch` package. They expose native function parameters and convert YAML lists to native tuple options where required. There is no nested SpectralMatch pipeline or implicit forwarding of stage outputs. Statistics functions are not registered. See the [function catalog](../api/plugins-spectralmatch.md).

For workflow file management, use explicit input/output file lists, as in the example. Native folder/glob/template forms remain available to direct Python callers; use `search_paths` or `create_paths` to calculate lists when needed. A value returned at runtime cannot determine that invocation's own output destinations; write planning expressions when HPC needs paths before execution. Compound options such as `vector_mask: [include, path, field]` need a matching `core:requires` path for dependency tracking and staging.

`resume_from_outputs` controls native reuse inside a function. Core handles reuse for explicit files separately. For tiled `merge_rasters`, set `output_tiles: true` and choose a directory destination; core invokes the function even if that directory exists, and transfers the directory recursively. `image_threads`, `concurrent_processing_backend` and `dask_scheduler` apply only to tiled merge; keep them `null` for a single mosaic. `compute_overviews` always invokes its function when requested because an existing raster alone does not prove it contains the requested overview levels.

The WorldView example uses VHR `core:calculate_overviews: true` on selected raster-producing steps, with levels from `shared.param:window_scales`. Native `build_overviews` remains unset (its default is false), and the separate `compute_overviews` alternative stays commented out. Footprints are calculated after radiometric matching so their image identifiers exactly match the seamline inputs.

## Explicit context files

Context persistence is opt-in on individual enabled steps, including steps without
a plugin. Each of the four controls maps filenames to one selector or a list:

```yaml
core:save_context:
  "path:./context/import_metadata.json":
    - dependencies
    - defined
    - var.metadata
    - const.band_order
```

Use dotted selectors for fields, including nested fields such as
`var.metadata.IMAGE_1`. `defined` selects that step's `var:`/`const:` assignments.
`dependencies` selects the context fields referenced recursively by that step's
parameters, including applicable shared parameters. The only spelling is
`dependencies`. Combined selections are deduplicated; selecting a parent includes
its children. The file contains selected resolved context values, not parameter
history or expressions. Save separate files before and after an assignment if
both states are needed.

Use `all` to save every current `const` and per-scene `var`, or load every field
present in the saved JSON:

```yaml
core:save_context:
  "path:./context/complete.json": all
core:load_context:
  "path:./context/complete.json": all
```

`all` also works with both upsert controls and in a list with other selectors.
Each scene keeps its own fields, including fields supplied by imports. Loading
matches scenes by ID and retains current fields absent from the file; fields
present in the file follow the replacement or merge rules below. Saving includes
all context available to that invocation, so selecting `all` can require earlier
steps to compute metadata that would otherwise be unused.

| Control | Action |
| --- | --- |
| `core:save_context` | Replace the destination with selected current context. |
| `core:save_upsert_context` | Insert/update selected data in the destination, retaining other entries. |
| `core:load_context` | Replace selected context values from the file. |
| `core:load_upsert_context` | Merge selected file data into the current context. |

Upserts merge objects recursively; incoming values win, lists and scalars replace.
Load operations are applied in workflow order during planning, before deciding
which processing is needed. They can establish scenes from the saved IDs. A missing
file logs a warning only when `core:log_to_console` is true, then continues without
loading or scheduling a producer. A present file with invalid structure or a missing
selected field is an error. Saves happen when their enabled step runs successfully;
they never force an otherwise unused function to run. A disabled step performs no
context I/O. Parent-owned writes combine individual scenes safely, including when
processing workers run concurrently.

An explicitly loaded import can supply its saved scene assignments without
rediscovering raw files. Include every needed field and the same scene-ID scheme.
Context values still needed by downstream functions must be supplied explicitly or
computed by the producing function. There is no implicit `<output>.context.json`
lookup and no automatic call to a plugin's `restore()` method.

See [saved file formats](../saved-file-formats.md#explicit-context-json) and
[HPC staging](../cli/hpc.md).

## Importing existing outputs

A scene-discovery step whose plugin declares `var_records_return` and
`var_path_return` may declare:

```yaml
core:satisfies: {cloud_mask: output_raster_path}
```

Keys identify named workflow steps, and values identify declared output parameters.
Each scene's discovered file supplies that output. Multiple imports merge by scene
ID. The binding updates the corresponding direct `var:`/`const:` output reference
so downstream steps receive the imported path. Duplicate conflicting bindings,
unknown step names and input/non-file parameter targets are rejected. Imported
products must pass the configured reuse checks. Without `core:satisfies`, imports
provide ordinary variables consumed through explicit parameter references.

`core:skip_plugin_call: true` keeps a step in the dependency plan without invoking
its function. Skipped calls form a suffix; returned values remain pending rather
than being replaced with fabricated results. HPC preparation injects this control
after its local cutoff. `core:scope: var` selects per-record execution. Without initialized records, default-scope functions can run
once with constants, while explicit `var` scope and `var:` access are rejected.
