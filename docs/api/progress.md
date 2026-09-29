# Workflow progress API

Workflow progress is a versioned JSON data contract. The backend collects events
from serial, process-pool and Dask workers; the `prompt_toolkit` dashboard consumes
these snapshots. Apps do not need to import the frontend, access scheduler internals, or parse
terminal output.

For durable task timings, use the related [statistics and timing API](statistics.md).
Its `event_callback` delivers every ended measurement without UI throttling;
`core:statistics_path` enables its append-only OpenTelemetry file consumer.

## Receive updates in Python

```python
from queue import SimpleQueue
from vhrharmonize import run_workflow

updates = SimpleQueue()
run_workflow("workflow.yml", progress_callback=updates.put)
# Another thread/event loop can consume updates and update an application UI.
```

`progress_callback(snapshot)` enables reporting independently of
`shared.core:show_progress`, which defaults to `true`. Set that control to `false`
to run without a display. `run_plugin()` and `Workflow.run()` accept the same
callback. A callback receives one positional snapshot dictionary. This workflow
API is distinct from a processing function's `progress_callback(**tqdm_fields)`;
the existing function callbacks and multiprocessing transport are unchanged.

Callbacks run in the parent process: the initial and final callbacks run on the
calling thread, and intermediate updates on the progress reader thread. Updates
are coalesced to at most four per second, plus immediate initial and final
snapshots. Short operations may finish between updates. Keep callbacks quick and
enqueue updates for GUI event loops. The backend gives each consumer an independent
copy. An exception disables that consumer and adds a recent message; processing
and other consumers continue. Processing exceptions still propagate normally,
after publishing a final `failed` snapshot.

## Request a snapshot

```python
from vhrharmonize import load_workflow, read_progress_snapshot

workflow = load_workflow("workflow.yml")
# In an application worker thread:
workflow.run(progress_callback=updates.put)
# From another thread during execution, or after it finishes:
snapshot = workflow.get_progress()

# From another process:
snapshot = read_progress_snapshot("workflow.yml.progress.json")
```

`get_progress()` returns an independent snapshot, or `None` before a reported run
or when reporting is disabled. It retains the final snapshot after success or
failure. Dry runs do not publish execution progress.

For reporting to a file without a display or Python callback:

```yaml
shared:
  plugin: shared
  core:run: true
  core:report_progress: true
  core:show_progress: false
```

Either control enables collection. When loading a YAML file, reporting writes
`<workflow.yml>.progress.json`. An explicit `progress_path="/path/status.json"`
on `run_workflow`, `run_plugin`, or `Workflow.run` also enables reporting and
overrides this destination. In-memory configurations have no default file.
The destination's parent directory must exist. Snapshots are replaced atomically
at most once per second, plus initial and final updates. File-write failures are
recorded as messages and do not stop processing. Messages are collected even with
the display off; usual core/plugin logging controls determine which messages
are emitted.

`read_progress_snapshot()` raises `FileNotFoundError` before a file exists, and
`ValueError` for invalid data or an unsupported schema version.

## HPC and rendering

```python
from prompt_toolkit import print_formatted_text
from vhrharmonize import get_slurm_progress, render_progress

snapshot = get_slurm_progress("configs/1.staged.hpc.yml")
if snapshot is not None:
    print_formatted_text(render_progress(snapshot), end="")
```

`get_slurm_progress()` fetches one snapshot over SSH, without printing, fetching
Slurm logs, or modifying the staged YAML. It accepts a staged YAML path or its
loaded mapping. It returns `None` when no snapshot is available, its data/version
is invalid, or its job ID belongs to another Slurm job. SSH execution errors
propagate. Apps can poll this function at their own interval.

`vhr hpc-progress --config configs/1.staged.hpc.yml` prints just the JSON snapshot
when available. `hpc-status` uses the same fetch API, retains the snapshot in its
staged YAML, and displays it through the same `prompt_toolkit` renderer as local execution.
`render_progress(snapshot, width=120, ascii_only=False)` returns `FormattedText`
and never starts a live display. Use `prompt_toolkit.formatted_text.to_plain_text()`
to obtain an unstyled string. `hpc-status` and redirected workflow output omit
terminal escape sequences; terminals without Unicode support use ASCII borders.
These consumers use numeric ETAs from the snapshot; they do not recalculate work.

Interactive workflow runs use a full-screen dashboard. Messages fill the available
space as a borderless, full-width main window without a heading. Mouse-wheel scrolling,
Up/Down, and PageUp/PageDown scroll messages while Workflow progress and
Active operation remain pinned below. Both tables put Progress in the rightmost
column. Home moves to the oldest retained message;
End resumes following new messages. Left/Right reveal long message lines.
The dashboard retains up to 10,000 message lines, restores the terminal on exit,
and prints one final summary. Processing stays on the calling thread, independently
of the UI thread. Noninteractive input/output and `TERM=dumb` use a static summary.

## Snapshot version 2

`ProgressSnapshot`, `ProgressRow`, `ActiveOperation`, and `ProgressCallback` are
public Python types in `vhrharmonize.progress`. `validate_progress_snapshot()`
validates and copies received data. Consumers should check `version`, tolerate
additional fields, and use the machine-readable values:

| Field | Meaning |
| --- | --- |
| `version` | Schema version, currently `2`; earlier experimental version 1 is unsupported. |
| `run_id` | Unique identifier for this execution, independent of the HPC run-directory name. |
| `job_id` | Slurm job ID, or `null` for a local run. |
| `updated_at` | UTC ISO timestamp of the snapshot. |
| `status` | `running`, `completed`, or `failed`. |
| `rows` | Ordered step rows, using configured step names. |
| `total` | Combined counts and active tasks across all rows, plus combined ETA. |
| `active` | Active tasks with `task_id`, `step`, `scene`, `stats`, and numeric `eta_seconds`. |
| `messages` | Up to five recent plain-text messages, without terminal formatting. |
| `message_history` | Optional bounded history: `sequence` counts messages appended during this run, and `messages` holds the most recent 1,000 messages. |

The last history entry has sequence number `message_history.sequence`. Consumers
can use this cursor to collect messages between screen refreshes, including
identical consecutive messages. A sequence gap larger than the history length
means earlier messages have expired. Reset the cursor when `run_id` changes.
Older version 2 snapshots without this optional field still render normally.

Each row contains `unused`, `done`, `run`, `all`, `reused`, `percentages`,
`fraction_done`, `active`, `pending`, `worker_progress`, `status`, and
`eta_seconds`. **Done means completed during this run** and overlaps Run.
`percentages` contains each count divided by All, multiplied by 100; these
overlapping columns do not form an exclusive partition. `fraction_done` is
Done / Run, or `null` when Run is zero. Total counts sum the step counts, so All
on the total row represents scene-step work rather than distinct scenes.

The dashboard labels the `reused` count **Loaded** and orders its count columns as Unused,
Loaded, Done, Run, All. Their headers and bar segments share fixed blue, purple,
green and light gray colors, respectively (All remains neutral). Cached
scenes that also count as unused occupy only the Loaded segment of the bar.
This display change preserves the version 2 snapshot fields.

Row statuses are `running`, `waiting`, `completed`, `reused`, `skipped`, or
`failed`. `eta_seconds` is a nonnegative number or `null` when not yet estimable.
Step ETAs use completed task durations and observed concurrency; the combined ETA
uses estimated remaining work and worker capacity. These are estimates, especially
across dependency barriers. `worker_progress: false` means detailed worker
callbacks are unavailable; lifecycle and completion accounting remain available.

Active `stats` retain tqdm's `n`, `total`, `prefix`, `unit`, `elapsed`, and `rate`,
plus any `operation`, `status`, and `scene` labels emitted by the function.
`total` and `rate` can be `null` for indeterminate operations. Elapsed and ETA are
seconds, rate is units per second. Queued tasks are not counted as active, and
Done advances only after core accepts the completed output.

The dashboard renders one row per active operation with Step, ID, Status, Elapsed,
ETA and Progress columns. It shows elapsed seconds and `TBD` for unknown ETAs. Detailed
phase descriptions and unit counts remain available in `stats` for API consumers.

::: vhrharmonize.progress
