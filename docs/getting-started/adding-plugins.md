# Adding plugins

Keep the Python function and its adapter in one module, such as `my_package/copy_file.py`:

```python
from pathlib import Path
from shutil import copy2
from vhrharmonize.plugins.base import FunctionPlugin


def copy_file(input_path: str, output_path: str):
    """Copy a file to an explicit destination, returning its size."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    copy2(input_path, output_path)
    return {"bytes": Path(output_path).stat().st_size}


class CopyFile(FunctionPlugin):
    target = "my_package.copy_file:copy_file"
    input_dependency_paths = {"input_path"}
    input_existence_check_paths = {"input_path"}
    input_protection_paths = {"input_path"}
    input_hpc_staging_paths = {"input_path"}

    output_path_resolution_paths = {"output_path"}
    output_dependency_paths = {"output_path"}
    output_parent_creation_paths = {"output_path"}
    output_collision_check_paths = {"output_path"}
    output_reuse_paths = {"output_path"}
    output_validation_paths = {"output_path"}
    output_invalid_removal_paths = {"output_path"}
    output_temporary_cleanup_paths = {"output_path"}
    output_hpc_staging_paths = {"output_path"}
    output_hpc_download_paths = {"output_path"}
    # Optional for raster products:
    # output_overview_calculation_paths = {"output_path"}
```

Each feature independently selects **function parameter names**, not literal filenames. Selected arguments can contain one path or a flat list of paths; each list item is tracked separately. All selections default to empty: keep only the features the plugin needs. Selecting path resolution does not enable validation, reuse, cleanup or transfers. For example, `output_overview_calculation_paths = {"output_image"}` restricts core overviews to that argument even if the function writes other files. See the [complete feature table](../configuration/workflow-config.md#plugin-file-features) for effects and YAML controls. Explicit `var:`, `const:`, `collect:` and JSONata references in parameters are tracked automatically. Context fields are never matched to function arguments by name. Keep optional imports inside functions so registration and recipe CLI help remain lightweight. Thin wrappers around another package can set `wrapper.__parameter_sources__ = (factory,)`, where `factory()` returns the delegated function. This exposes its native header and docstrings to argument validation and generated function options without copying processing code.

Register the adapter in your package's `pyproject.toml` and install with `pip install -e .`:

```toml
[project.entry-points."vhrharmonize.plugins"]
my_copy = "my_package.copy_file:CopyFile"
```

Its registered name selects the implementation through `plugin:` and becomes a CLI command. The top-level step name is independent and must be unique:

```yaml
shared:
  plugin: shared
  core:run: true
  # core:protect_source_files: true
  core:output_metadata_path: expr:var.output_dir & '/processing.json'
  # core:delete_final_json_first: true

discover_inputs:
  plugin: import_files
  core:run: true
  param:search_glob: /data/files/*.tif
  param:output_dir: ../output
  # param:temp_dir: sys
  var:image: returned:file_path
  var:basename: expr:$replace($split(var.file_path, '/')[-1], /\.[^.]*$/, '')
  var:suffix: ""

copy_products:
  plugin: my_copy
  core:run: true
  core:require_outputs: param:output_path
  const:copy_suffix: _copied
  param:input_path: var:image
  var:suffix: expr:var.suffix & const.copy_suffix
  var:copied: expr:var.output_dir & '/' & var.basename & var.suffix & '.tif'
  param:output_path: var:copied
  var:copied_bytes: returned:bytes
```

Later steps can read `var:copied` and `var:copied_bytes`. JSONata sees `{"const": {...}, "var": {...}}`: constants belong to the workflow, variables to a scene. Scene and aggregate plugins can assign either scope, including from returned values. An adapter with `scope = "aggregate"` runs once; its ordinary `var:` assignments still evaluate separately for each scene, and its parameters use `collect:` for scene values. Both scopes support `collect:image` or `collect:$` after discovery. Scene-independent constants resolve once per step; constants derived separately by scene invocations must agree. See [assigning shared and scene values](#assigning-shared-and-scene-values) below. Disabled steps do nothing and export no values; all steps default to disabled.

```bash
vhr my_copy --config recipe.yml
vhr workflow --config recipe.yml --dry-run
```

Python callers can call `copy_file(...)` directly or use the recipe API:

```python
from vhrharmonize import run_plugin, run_workflow

run_plugin("my_copy", "recipe.yml")
run_workflow("recipe.yml")
```

Validation belongs in the function/API so Python and CLI callers get the same errors. Return JSON-compatible data. Argument precedence is **explicit plugin `param:` > shared `param:` > function default**. Shared `param:` settings are offered to all adapters; only supported arguments reach the function. Use `FunctionPlugin.aliases` to translate external API names. Pass stored fields explicitly, for example `param:gain: const:gain`.

No context is injected into functions or processing adapters. `context` and `variables` are ordinary argument names. To pass an entire object, link it explicitly:

```yaml
param:scene: var:$       # Entire scene variable object.
param:settings: const:$  # Entire workflow constant object.
param:context: expr:$    # Both scopes: {const: {...}, var: {...}}.
```

Use the parameter names your function accepts. Dependencies are tracked from these references automatically, including returned values supplied by explicit context loads. `"var:"` and `"const:"` also select whole scopes, but need YAML quotes; the `$` forms above do not. Passed objects are snapshots. Publish changes through returned values and YAML assignments, not by mutating arguments.

No function must implement every shared parameter. The engine owns scheduling and the file features selected by the adapter. Proactive cleanup covers selected temporary regular files after their consumers finish; private caches and directory outputs remain the plugin's responsibility. Add `custom_nodata_value`, `epsg`, `output_dtype` or other options when meaningful for the function.

There is one adapter class: `FunctionPlugin`. If custom execution is needed, override `run(self, *, params, shared)`. All resolved arguments—including input paths and output paths—are in `params`. The independent file-feature declarations determine what core does with each argument. Leave core reuse, validation, invalid-removal or overview selections empty for operations the function handles itself. Core never skips an invocation merely because its output directory exists.

## Creating or replacing scenes

Any ordinary function can return a list of plain dictionaries and declare which returned field contains it:

```python
def find_items():
    return {"items": [{"name": "first"}, {"name": "second"}]}


class FindItems(FunctionPlugin):
    target = "my_package.discovery:find_items"
    var_records_return = "items"  # Use "$" if the function returns the list itself.
```

After registering this plugin as `find_items`:

```yaml
find_items:
  plugin: find_items
  core:run: true
  var:label: returned:name

my_processor:
  plugin: my_processor
  core:run: true
  core:require_outputs: true
  param:name: var:label
```

Before a scene-setting function returns, `var` is unavailable: reads, assignments, JSONata access and `collect:` raise `ValueError`. Other enabled functions can run once with constants and ordinary parameters. `shared.var:` is therefore invalid. Whole-context expressions (`expr:$`) contain only `const` before initialization.

The scene-setting function runs once with aggregate scope. Each dictionary becomes one scene's `var` object; optional `var:` assignments map each dictionary afterward, with `returned:` selecting from that dictionary. Plugin `const:` assignments refer to the whole function result. Existing constants remain available unless explicitly reassigned. By default, a later scene-setting function replaces the list, and core rebuilds scene indexes, contexts and remaining dependencies. Functions need no scheduling bookkeeping. In this mode an empty list removes all scenes; subsequent aggregate steps still run. Set `var_records_mode = "merge"` with a stable `var_id_return` to add new scenes and fill missing fields in existing scenes. Existing values win recursively, including nulls and lists. An empty merge leaves existing scenes in place; source-protection paths are combined. Previously completed and cached values survive the update.

`import_files` uses this same contract with `var_records_mode = "merge"` and `var_id_return = "scene_id"`; the function computes each returned ID from its ordinary `scene_id` parameter (configured as `param:scene_id` in YAML). The adapter only selects the returned field; it does not decide which input field or expression supplies the ID. Scene-setting functions run during planning, including dry-run and HPC preparation, once their required inputs and preceding work are available. HPC preparation limits calls to the explicitly selected local prefix; `core:skip_plugin_call: true` prevents discovery calls as well as ordinary processing calls. A scene setter depending on unfinished processing remains pending until execution supplies its inputs. Keep discovery functions lightweight: dry-run can invoke them, including their own file-writing behavior.

Optional declarations:

| Declaration | Default | Purpose |
|---|---|---|
| `var_records_mode` | `"replace"` | `"replace"` resets scenes; `"merge"` adds scenes and missing fields while preserving existing values. Merge requires `var_id_return`. |
| `var_id_return` | `None` | Field within each returned scene dictionary for its unique ID; otherwise use the list index. |
| `var_path_return` | `None` | Primary file path within each returned scene, used by `core:satisfies` and per-scene HPC assignments. Nested fields such as `asset.uri` are supported. |
| `source_file_protection_paths_return` | `None` | Field within each scene dictionary containing original paths to protect. Controlled by `shared.core:protect_source_files`, default `true`. |
| `constant_values_return` | `None` | Returned constant dictionary. Replace mode updates fields; merge mode recursively fills only missing fields. Explicit YAML `const:` assignments take precedence. |
| `temporary_directory_context_paths` | `()` | Ordered JSON locations containing temporary roots, for example `("var.paths.work", "const.paths.work")`. |
| `output_directory_context_paths` | `()` | Ordered JSON locations containing output roots, for example `("var.paths.products", "const.paths.products")`. |
| `directory_parameters` | `{}` | Map a declared directory location to the function parameter that supplies it, for example `{"var.paths.products": "destination"}`. HPC staging updates that parameter when relocating an implicit directory. |
| `discovery_input_parameter` | `None` | Parameter name used to label discovery uploads, for example `"catalog"`. Without a declaration, the label is `discovery`. |

Returned-field selectors support dotted fields and `$` for the whole result. Directory selectors start with `const.` or `var.` and may point to a path or list of paths. Core resolves existing declared locations; it does not create context fields or invent directory defaults. The workflow YAML directory is the base for selected relative output arguments. Cleanup applies only under populated temporary roots and retains explicitly required outputs. Roots do not request processing products.

`import_files` supplies directory defaults through its ordinary `temp_dir`, `output_dir`, `temp_dir_scope` (default `const`) and `output_dir_scope` (default `var`) arguments. Its adapter registers `var.temp_dir`/`const.temp_dir` and `var.output_dir`/`const.output_dir`. Other plugins can publish completely different field names and declare those locations. Input arguments are passed unchanged by core: the importer resolves metadata and companions relative to each found file and returns absolute paths for later functions.

Final JSON saving is independent of plugin return selectors: `shared.core:output_metadata_path` resolves against the completed context and appends `{const, var}` to a JSON array. The same destination collects multiple completions; scene-specific destinations naturally produce separate files. `core:delete_final_json_first` defaults to `true`, replacing old contents on the first write to each resolved destination during a workflow run. See [final metadata](../configuration/workflow-config.md#final-metadata-json).

Context persistence is configured explicitly with the four `core:*context` controls,
not a plugin checkpoint declaration or automatic `restore()` hook. Runtime scene
setters execute when their scene results are needed. HPC preparation requires
known scenes and paths; use the named local preparation cutoff or explicit context
loads to supply them. Earlier temporary files are retained across a scene reset
until their future readers are known.

### Custom discovery staging

Core uses the declarations above without checking plugin names. `core:satisfies`
requires `var_records_return` and `var_path_return`; saved context can supply
the declared primary field or a YAML variable assigned directly from that returned
field. Per-scene HPC assignments use a declared primary path, or a declared scene
ID when no suitable primary path is available. No `file_path` field is required.

An optional `stage_settings(...)` method handles discovery parameters that need
more than ordinary root substitution, such as globs or metadata lookup rules:

```python
def stage_settings(self, *, settings, params, returned, path_mappings,
                   file_paths, discovery_paths, config_dir):
    from vhrharmonize.workflow.values import remap_paths

    return {"param:catalog": [
        remap_paths(item["asset"]["uri"], path_mappings)
        for item in returned["items"]
    ]}
```

`settings` contains the selectively rewritten YAML block; `params` contains its
resolved local function parameters; `returned` is the original discovery result.
`path_mappings` maps local roots/files to remote destinations. `file_paths` is the
set of individually mapped local files, and `discovery_paths` contains all source
files discovered during planning. `config_dir` is the local workflow directory.
All mutable arguments are detached copies. Return only `param:name` overrides;
the default implementation returns `{}`. Core applies these updates through the
same comment-preserving YAML writer. A discovery step restored from explicit
context does not call this hook, since its discovery function will not rerun.

For a directory parameter that varies by scene, core supplies a `literal:expr:`
lookup using the declared scene fields. The discovery function must resolve that
expression within each scene, or its staging method can translate it into the
function's own input format. `ImportFiles` implements its glob and companion-rule
rewrites in this hook; those rules are not part of core.

See [configuration](../configuration/workflow-config.md) for controls and binding rules.


A pluginless step can update context before the next function:

```yaml
setup_values:
  core:run: true
  const:gain: 2
```

Only `plugin: shared` is a core-provided special implementation. It defines workflow defaults, not a callable processing function or CLI command. Pluginless steps cannot have `param:` or `returned:` values. Ordinary registered plugins require no core changes.

### Assigning shared and scene values

Plugins need no extra declaration to write either scope through YAML. `const:` stores a shared JSON value; `var:` stores one value per scene. Ordinary `var:` assignments evaluate per scene in every step; literal lists and dictionaries remain intact per scene. Aggregate function parameters cannot directly reference `var`, including through JSONata; use `collect:`. Aggregate `returned:` assignments to `var:` must select lists matching scene count or dictionaries keyed by exactly the current scene IDs. Core distributes those results and rejects mismatches. `collect:name` reads all current scene values and can also be passed by scene steps. Scene invocations writing the same shared constant must agree; use `const:values: collect:name` to preserve different scene values as a shared list.

### User-selected deliverables

Adapters do not populate `output_target_paths`. Core fills that selection from each step's `core:require_outputs` setting (`false` by default, `param:name`, a list of parameter selectors, or `true`). Continue declaring output dependencies, reuse, validation and transfer capabilities independently. Outputs inside temporary roots can be requested explicitly and are then protected from cleanup. Paths can be normalized in YAML with `path:var:name`, `path:./file` or `path:sys`; directory roles do not select the relative-path base.

`core:scope: var` selects per-record execution; `aggregate` selects one call for
the collection. Plugins declare `var_records_return`, `var_records_mode`,
`var_id_return`, and `var_path_return`. The former `scene` scope and `scene_*`
declarations are rejected. Return-field names are selected by those declarations
and belong to the function's data schema.
