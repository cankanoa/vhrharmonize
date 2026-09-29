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
  core:statistics_path: ./statistics/history.jsonl
```

The default is `null` (no statistics file). Relative paths use the workflow YAML's
directory, or `config_dir` for in-memory recipes. Parent directories are created.
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

# Leave core:statistics_path unset when supplying the recorder yourself.
with StatisticsRecorder("history.jsonl") as recorder:
    run_workflow("workflow.yml", event_callback=recorder)
```

## HPC files

HPC preparation rewrites `statistics_path` beneath the remote output directory and
registers it for download. It never uploads local history over remote history.
Executions using the same remote destination append there; a different HPC run
directory has its own history. Normal download conflict policy applies: use
`override_download_conflict: yes` to refresh a local mirror from the remote file.
Downloading replaces that mirror; it does not merge unrelated local and remote
histories. Keep separate destinations for separate remote runs when retaining all
their files locally. Run `vhr statistics` after download or on the execution host.

::: vhrharmonize.statistics
