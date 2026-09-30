# Saved file formats

This page describes files written by VHRHarmonize core and links to processing
libraries' own format documentation. Paths shown below are defaults or examples;
processing output names come from the workflow recipe.

## Timing history: `statistics.jsonl`

`core:save_statistics_path` appends one completed OpenTelemetry Python SDK span per
line. `core:load_statistics_path` reads the same format for runtime estimates. Both
default to `statistics.jsonl` beside the YAML. A file can contain many runs; group
them by `attributes["vhr.run_id"]`. There is no enclosing JSON array.

The example below shows a complete task record, pretty-printed for readability.
The SDK defines the outer structure. All `vhr.*` attribute names and their meanings
are VHR-specific. `service.name`, `host.name`, `process.pid` and the optional
`error.type` use OpenTelemetry attribute conventions; VHR supplies their values.
This is SDK span JSON, not OTLP wire JSON.

```json
{
  "name": "alignment",
  "context": {
    "trace_id": "0x11111111111111111111111111111111",
    "span_id": "0x2222222222222222",
    "trace_state": "[]"
  },
  "kind": "SpanKind.INTERNAL",
  "parent_id": "0x3333333333333333",
  "start_time": "2026-09-29T10:00:00.000000Z",
  "end_time": "2026-09-29T10:00:45.000000Z",
  "status": {"status_code": "OK"},
  "attributes": {
    "vhr.schema_version": 1,
    "vhr.run_id": "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
    "vhr.kind": "task",
    "vhr.name": "alignment",
    "vhr.status": "completed",
    "vhr.duration_seconds": 45.0,
    "vhr.backend": "process_pool",
    "vhr.processing_direction": "vertical",
    "vhr.concurrent_processing": 4,
    "vhr.config": "/project/workflow.yml",
    "vhr.job_id": "123456",
    "vhr.plugin": "alignment",
    "vhr.scene": "scene_P004",
    "vhr.scene_units": 1,
    "vhr.task_id": "1:0",
    "vhr.measured": true
  },
  "events": [],
  "links": [],
  "resource": {
    "attributes": {
      "service.name": "vhrharmonize",
      "host.name": "compute-node-01",
      "process.pid": 12345
    },
    "schema_url": ""
  }
}
```

| Record / optional fields | Structure and meaning |
| --- | --- |
| `vhr.kind: workflow` | One root span per run, with `parent_id: null`; overall elapsed duration. |
| `vhr.kind: core` | Core-stage measurement such as discovery, planning or execution. |
| `vhr.kind: task` | One plugin call; scene/unit fields describe its work. Aggregate calls can represent multiple scene units. |
| `vhr.kind: step_summary` | Final integer counts: `vhr.unused`, `vhr.reused` (Loaded), `vhr.done`, `vhr.run`, `vhr.all`. Its raw duration is zero; analysis excludes it as a timing sample. |
| `vhr.phase` | `preflight` for calls during scene discovery; these may omit `vhr.task_id`. |
| `vhr.config`, `vhr.job_id` | Omitted when the workflow has no filename or no Slurm job ID. |
| `error.type` | Exception class when known; status becomes `ERROR`. |
| `status.description` | Optional SDK explanation; used when recording ends without a workflow completion event. |
| Incomplete workflow | May omit duration and execution attributes; retained for identifying unfinished runs. |

`vhr.status` is `completed`, `failed`, `incomplete`, `reused`, `skipped` or `waiting`.
`vhr.kind` is separate from the SDK's `kind`. The resource identifies the recording
parent process, even when task durations were measured by remote workers.

`validate_statistics_record(record)` validates a single record. Loading and analysis
share this validator and add filename/line numbers to errors. Duplicate trace/span
IDs are counted once. The adjacent `.lock` file coordinates writers and contains no
statistics. See [statistics and timing](api/statistics.md) for the API, estimate
calculation and append behavior, and the
[SDK JSON serializer](https://opentelemetry-python.readthedocs.io/en/latest/sdk/trace.html#opentelemetry.sdk.trace.ReadableSpan.to_json).

## Statistics report: a separate summary JSON

`vhr statistics` / `summarize_statistics()` writes pandas' `orient="table"` format:

```jsonc
{
  "schema": {
    "fields": [
      // One {"name": "column_name", "type": "string|integer|number"}
      // descriptor per column below; types reflect the generated table.
    ],
    "pandas_version": "1.4.0"
  },
  "data": [
    {
      "kind": "task",
      "step": "alignment",
      "plugin": "alignment",
      "backend": "process_pool",
      "status": "completed",
      "runs": 1,
      "unfinished_runs": 0,
      "records": 1,
      "samples": 1,
      "total_seconds": 45.0,
      "mean_seconds": 45.0,
      "std_seconds": null,
      "min_seconds": 45.0,
      "p50_seconds": 45.0,
      "p90_seconds": 45.0,
      "p95_seconds": 45.0,
      "max_seconds": 45.0,
      "scene_units": 1,
      "done": 0,
      "run": 0,
      "all": 0,
      "reused": 0,
      "unused": 0
    }
  ]
}
```

`--per-run` adds `run_id` to the row and schema. Pandas defines the container/schema
format; VHR defines the columns and grouping. `pandas_version` describes pandas'
table-schema revision. Missing measurements become `null`, including standard
deviation for one sample. Progress counts belong to `step_summary` rows. The report
is replaced atomically; raw history stays unchanged. Load it with
`pandas.read_json(path, orient="table")`. Use the raw JSONL, not this report, as
`load_statistics_path`.

## Progress snapshot: `<workflow.yml>.progress.json`

This VHR-specific JSON object is replaced atomically as execution progresses:

| Field | Structure |
| --- | --- |
| `version` | `2` |
| `run_id`, `updated_at`, `job_id`, `status` | Run ID, UTC timestamp, nullable Slurm job ID, run status. |
| `total` | One progress row for the whole workflow. |
| `rows` | Array of progress rows for configured steps and core overview steps. |
| `active` | Array of `{task_id, step, scene, stats, eta_seconds}` objects. |
| `messages` | Last five message strings. |
| `message_history` | `{sequence, messages}`: cumulative sequence number and last 1,000 messages; optional in older snapshots. |

Each row contains `name`, `unused`, `reused`, `done`, `run`, `all`, `percentages`,
`fraction_done`, `active`, `pending`, `worker_progress`, `status` and `eta_seconds`.
`percentages` maps the five count keys to percentages of `all`; `fraction_done`
is Done/Run or `null`. ETAs are seconds or `null`. An active operation's `stats`
contains tqdm-style `n`, `total`, `prefix`, `unit`, `elapsed`, `rate`, with optional
`operation`, `status` and `scene`. Historical estimates do not change this schema.
See [the progress API](api/progress.md) and `validate_progress_snapshot()`.

## Final metadata JSON

`core:output_metadata_path` writes an array of completed contexts:

```json
[
  {
    "const": {"output_dir": "/project/output"},
    "var": {"scene_id": "P004", "result": "/project/output/P004.tif"}
  }
]
```

Keys inside `const` and `var` come from the recipe and plugins. A workflow without
scenes has only its constant scope. `core:delete_final_json_first` controls whether
the first write replaces previous contents. See
[final metadata](configuration/workflow-config.md#final-metadata-json).

## Explicit context JSON

The four step controls `core:save_context`, `core:save_upsert_context`,
`core:load_context`, and `core:load_upsert_context` use the same simple format:

```json
{
  "const": {"band_order": ["red", "green", "blue"]},
  "scenes": {
    "P004": {"metadata": {"cloud_cover": 0.1}, "basename": "P004"},
    "P005": {"metadata": {"cloud_cover": 0.2}, "basename": "P005"}
  }
}
```

`const` contains selected workflow-wide values once. `scenes` maps stable scene IDs
to selected per-scene variable fields. Either object can be empty or absent on
load. These are resolved values, without expression definitions, parameter dumps,
assignment histories, or automatic cache fingerprints. Use different filenames to
save different points in the workflow. A filename-keyed mapping selects fields
using `var.metadata`, `const.band_order`, `dependencies`, `defined`, or `all`.
`all` saves every current constant and scene variable, or loads every field present
in the file. It uses this same JSON structure and supports scenes with different
sets of fields, including empty variable objects.

`save_context` replaces the old file; saves from separate scenes in the same step
are collected into that run's file. `save_upsert_context` preserves other entries.
`load_context` replaces selected values; `load_upsert_context` merges dictionaries
recursively. Incoming leaves win, lists/scalars replace, and scenes match by ID.
Writes are atomic and use an adjacent lock file. Invalid JSON, malformed structure,
or missing selected fields are errors. A missing file only emits a warning when
console logging is enabled and is otherwise ignored. See
[explicit context files](configuration/workflow-config.md#explicit-context-files).

The former implicit `<output>.context.json` mechanism is removed. Old sidecars are
neither automatically loaded nor deleted.

## Processing outputs and atmosphere JSON

Raster products (including fetched DEMs) use their selected GDAL format, typically
GeoTIFF with bands, nodata, a geotransform and a CRS. Vector outputs use the selected
GIS format and plugin-defined attribute fields. TIFF overviews are raster resolution
pyramids, not separate JSON records.

`fetch_atmosphere` saves its function's returned mapping as JSON. Fields and units
depend on the selected source (`nasa_power` or `modis_gee`); see the return structures
in [fetch_atmosphere](api/plugins-fetch-atmosphere.md). Converted sensor metadata
preserves the imported document's structure rather than imposing a new core schema;
see [metadata I/O](api/io-metadata.md).

## Generated HPC files and logs

The staged workflow YAML is a selectively edited copy of the original recipe,
with comments and expressions retained. Explicit context JSON files may have
separate generated copies with paths rebased for the remote host; scene metadata
is not embedded in the YAML.
File mappings may replace import globs with remote filename lists and add path-only
lookup expressions for companion rules whose directory layout changed. Directory
root mappings continue to preserve relative layouts and ordinary expressions.

The staged HPC YAML preserves user settings/comments and adds resolved control
paths, SSH settings, job/status fields, and local-to-remote transfer maps:
`uploaded_input_paths`, `uploaded_reference_paths`, `download_output_paths`, and
`download_log_paths`. `upload_groups` is a list of `{step, variable, files}` objects:
the originating step, the mapped path reference (or parameter/control label), and
local upload paths. Each path occurs in only one group, independently of scene
count. It labels the upload dashboard without changing transfer selection.
The YAML may also contain `workflow_progress` using the snapshot
structure above. A `.prepare.yml` copy records an optional named local cutoff.
The generated `.sbatch` file is a shell script. Slurm logs are ordinary text,
including the final static progress report. See [HPC commands](cli/hpc.md).

## SpectralMatch files

See **[SpectralMatch's File Formats and Input Requirements](https://spectralmatch.github.io/spectralmatch/formats_and_requirements/)**
for its saved-file structures. Those definitions are maintained there and are not
duplicated here:

- [Regression parameters](https://spectralmatch.github.io/spectralmatch/formats_and_requirements/#regression-parameters-file)
- [Tie-point adjustments](https://spectralmatch.github.io/spectralmatch/formats_and_requirements/#tie-point-adjustments-file)
- [Tie-point GIS export](https://spectralmatch.github.io/spectralmatch/formats_and_requirements/#tie-point-gis-export)
- [Block maps](https://spectralmatch.github.io/spectralmatch/formats_and_requirements/#block-maps-file)
