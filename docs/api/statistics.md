# Statistics and timing API

Core measures workflow stages and plugin invocations, including aggregate steps,
overviews and scene discovery. Plugins need no timing code or progress callback.
The recorder consumes unthrottled events from the same core reporting backend as
the progress API; it does not scrape logs or sample dashboard refreshes.

## Keep raw history across runs

```yaml
shared:
  plugin: shared
  core:run: true
  core:save_statistics_path: statistics.jsonl
  core:load_statistics_path: statistics.jsonl
```

Both defaults are `statistics.jsonl` beside the workflow YAML, or in `config_dir`
(the working directory by default) for in-memory recipes. Set either to `null` to
disable that operation independently. `statistics_path` has been renamed to
`save_statistics_path`; update existing recipes. Relative paths use the workflow YAML's
directory, or `config_dir` for in-memory recipes. Parent directories for saving are created.
Every real run appends; existing history is never cleared. Recording works with
both `core:show_progress: false` and `core:log_to_console: false`. It also enables
the usual progress snapshot backend.

The file contains one OpenTelemetry span per line using the Python SDK's
[`ReadableSpan.to_json()` format](https://opentelemetry-python.readthedocs.io/en/latest/sdk/trace.html),
written by its `ConsoleSpanExporter` to a locked append stream. This is **SDK span
JSON, not OTLP wire JSON**. No collector or network service is needed, and VHR does
not replace an application's global OpenTelemetry provider. Install VHR's updated
dependencies on the execution host; the recorder stays in the parent process.

Each run has a new `vhr.run_id` (also used by progress snapshots) and a separate
OpenTelemetry trace. The workflow span is the root; its children describe core
stages, task calls and final step counts. Records contain timestamps, duration,
outcome, step/plugin names, scene identifiers, configured concurrency/backend,
host/process and, when available, config path and Slurm job ID.

[`filelock`](https://py-filelock.readthedocs.io/en/latest/) protects each complete
line against concurrent VHR writers using the same file. The adjacent `.lock`
file is normal. On HPC, use storage that supports file locking. Records are
flushed as spans finish. An unwritable path fails before execution; a later
consumer/write error warns and disables that consumer without stopping processing.

## Estimate runtime from history

Before appending the current run, core reads `load_statistics_path`. A missing file
is an empty history, allowing the first run to create it. An existing malformed file
raises an error with its filename and line number before new statistics are appended.

Successful measured task calls are pooled by configured step name, plugin and core
backend. The estimate is total measured seconds divided by total scene units;
aggregate calls retain their scene weight. Duplicate spans, failed/incomplete calls,
cached outputs and count-only records do not contribute. Completed current-run calls
join those samples as execution proceeds. Core uses the resulting means for remaining
step/workflow time and for active operations lacking a callback rate. Measurable
callback ETAs take precedence for the active operation. No history or current samples
means `TBD`; pending scene discovery also remains `TBD`.

The same estimates appear in console execution logs and the progress API consumed by
the dashboard and HPC status. Step estimates account for configured/observed concurrent
calls; total time is approximate because dependency barriers and shared resources can
limit concurrency. These are task-runtime estimates, not estimates of startup or cleanup.
Pick different filenames yourself for different input sizes/types; core does not classify
inputs or choose history files automatically. Save and load paths may be different.

`validate_statistics_record(record)` validates one raw SDK span and its VHR attributes,
returning it unchanged or raising `ValueError`. Both history loading and summary generation
use this validator. `load_statistics(path)` returns `(seconds, scene_units)` totals keyed
by `(step, plugin, backend)`. See [saved file formats](../saved-file-formats.md) for
record structures and which fields belong to OpenTelemetry versus VHR.

## Generate a statistics file

```bash
vhr statistics --input-path statistics/history.jsonl --output-path statistics/summary.json
vhr statistics --input-path statistics/history.jsonl --output-path statistics/by-run.json --per-run
vhr statistics --input-path statistics/history.jsonl --output-path statistics/one-run.json --run-id RUN_ID
```

Command paths are relative to the current directory. It writes a separate summary
file atomically, replacing the previous summary, and prints the same JSON. The
raw input remains unchanged; using the raw file as the output is rejected.

```python
import pandas as pd
from vhrharmonize import summarize_statistics

report = summarize_statistics("history.jsonl", "summary.json", per_run=True)
table = pd.read_json("summary.json", orient="table")
step_timings = table[table["kind"] == "task"]
```

The report uses pandas' standard
[`orient="table"` JSON format](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_json.html),
with `schema` and `data`. Rows group by `kind`, configured `step` name, `plugin`,
`backend` and `status`; `--per-run` adds `run_id`. Runs with the same grouping
values are pooled, including differing worker counts, unless `--per-run` is used.

| Columns | Meaning |
| --- | --- |
| `runs`, `unfinished_runs` | Contributing runs, and those with no completed/failed workflow record. |
| `records`, `samples` | Raw spans and usable duration measurements. |
| `total_seconds`, `mean_seconds`, `std_seconds` | Sum, arithmetic mean and sample standard deviation. |
| `min_seconds`, `p50_seconds`, `p90_seconds`, `p95_seconds`, `max_seconds` | Range and pandas' interpolated duration percentiles. |
| `scene_units` | Scene weights represented by task invocations; aggregates can cover several scenes. |
| `unused`, `done`, `run`, `all`, `reused` | Final progress counts on `step_summary` rows, summed across runs. |

`done` means completed during that run and overlaps `run`. Cached/skipped work
contributes counts, not zero-duration timing samples. Discovery-only steps have
task timings but do not appear in processing count rows. An aggregate call is one
timing sample, regardless of its scene weight. Standard deviation is `null` for a
single sample; all timing statistics are `null` when no measurements exist.
Duplicate trace/span IDs are ignored. Empty input produces an empty table.
Malformed/truncated records cause a line-numbered error without replacing the
previous summary.

## Timing boundaries

* `task`: worker duration around the plugin call, excluding scheduler queue time.
  A call is completed only after core validates its outputs. Validation failures
  become failed calls. Preflight calls used for scene discovery are measured in
  core, including their immediate finalization where applicable. No plugin-internal
  phases are inferred.
* `core`: initialization, discovery, graph build, actual planning, execution and
  final cleanup. These are inclusive, sometimes nested measurements; do not add
  them together for total runtime.
* `workflow`: elapsed time from workflow construction through execution and final
  reporting. With `load_workflow()` followed by a later `run()`, this includes the
  intervening time. Subsequent runs on that object get a new clock and run ID.
* `step_summary`: counts only, with no duration sample in the report.

Durations use monotonic clocks on the machine doing the work; UTC nanosecond start
times are for correlation. Concurrent task durations overlap, so their sum is
**processing time across tasks**, not workflow wall time. `workflow` rows measure
elapsed time. Failed calls are grouped separately from successful calls. Tasks
interrupted without a worker result are incomplete and excluded from timing
statistics. Completed spans in a killed run remain usable and contribute to
`unfinished_runs` until a final workflow span exists.

Construction/planning measurements are buffered until `run()` attaches consumers.
Dry runs and HPC preparation do not write statistics. Errors that prevent workflow
construction or occur before `run()` is entered cannot be recorded.

## Consume timing events in another application

```python
from queue import SimpleQueue
from vhrharmonize import run_workflow

events = SimpleQueue()
run_workflow("workflow.yml", event_callback=events.put)
```

`run_workflow()`, `run_plugin()` and `Workflow.run()` accept
`event_callback(event)`. It works without a statistics file or terminal display. Events use
the public `TimingEvent` / `TimingCallback` types in `vhrharmonize.progress`:

| Field | Meaning |
| --- | --- |
| `version` | Timing contract version `1`, independent of snapshot version. |
| `run_id`, `run_started_ns` | Run identity and UTC start time in nanoseconds. |
| `kind`, `name` | Measurement kind and configured step/core-stage name. |
| `start_time_ns`, `duration_seconds` | UTC start and monotonic elapsed duration. |
| `status` | `completed`, `failed`, `incomplete`, `reused`, `skipped`, or `waiting` for unresolved counts. |
| `attributes` | OpenTelemetry-compatible details using `vhr.*` keys plus `error.type` when known. |

Each ended measurement is delivered once, independently of UI refresh intervals.
Callbacks receive detached data and are serialized in the parent process, on the
calling or progress reader thread. Keep them quick; enqueue for UI/event loops.
Callbacks are never sent to multiprocessing workers. An exception disables that
callback with a warning; other consumers and processing continue.

For an application-managed destination, use the public consumer directly:

```python
from vhrharmonize import StatisticsRecorder, run_workflow

# Set core:save_statistics_path: null when supplying the recorder yourself.
with StatisticsRecorder("history.jsonl") as recorder:
    run_workflow("workflow.yml", event_callback=recorder)
```

## HPC files

HPC preparation rewrites `save_statistics_path` beneath the remote workspace products directory and
registers it for download. `load_statistics_path` is staged separately beneath the remote
inputs directory and uploaded when the local file exists. Preparation does not write
timings. It never uploads local history over the remote append destination, even when
the local load and save paths are identical. Each prepared job uses that input snapshot
for initial estimates; download updated history before preparing the next job to include
new measurements. A missing history input remains optional.
Executions using the same remote destination append there; a different HPC run
directory has its own history. Normal download conflict policy applies: use
`override_download_conflict: yes` to refresh a local mirror from the remote file.
Downloading replaces that mirror; it does not merge unrelated local and remote
histories. Keep separate destinations for separate remote runs when retaining all
their files locally. Run `vhr statistics` after download or on the execution host.

::: vhrharmonize.statistics
