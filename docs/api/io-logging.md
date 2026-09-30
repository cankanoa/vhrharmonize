# Logging and progress callbacks

When `core:log_to_console` is true, core prints `[core:workflow] Start` before input
discovery or processing, then announces its discovery, build, planning, execution
and enabled final cleanup stages with the same format. Planning announces only
actual work, including replanning after scene changes. With progress reporting
active, execution-time logs are captured in snapshot messages; early startup logs
print immediately before the dashboard is initialized.

For persistent timings, set `shared.core:save_statistics_path` to an append-only JSONL
file. It and `core:load_statistics_path` default to `statistics.jsonl` beside the YAML;
the latter seeds runtime estimates in execution logs and progress snapshots. Set either
to `null` to disable it. Core's unthrottled `event_callback` API reports measurements independently
of console logging and UI refreshes; the OpenTelemetry recorder is one consumer.
Use `vhr statistics` or `summarize_statistics()` to produce a separate report.
See [statistics and timing](statistics.md) for the formats and Python API.

Processing functions accept an optional Python-only `progress_callback`. The same callback contract works for standalone calls, workflow tasks and custom plugins. It receives keyword fields from tqdm's `format_dict`, including:

| Field | Meaning |
| --- | --- |
| `n` | Completed units in this operation |
| `total` | Expected units, or `None` for an indeterminate operation |
| `prefix` | Operation description (tqdm's `desc`) |
| `unit` | Unit label, such as `tiles`, `bands` or `images` |
| `elapsed` | Elapsed seconds |
| `rate` | Units per second, or `None` before a rate is available |
| `operation` | Identifier for a bar; changes when a new phase starts |
| `scene` | Optional input image identifier for work inside an aggregate call |

Accept `**stats` so additional tqdm fields remain compatible. Lifecycle failures also include `status="failed"` and never report successful completion.

```python
from vhrharmonize import align_image_pair

def report(**stats):
    print(stats["prefix"], stats["n"], stats["total"])

align_image_pair(
    "moving.tif", "reference.tif", "aligned.tif",
    progress_callback=report,
)
```

A callback is a Python callable, not a YAML expression or CLI option. With `core:report_progress: true`, `core:show_progress: true`, or an application progress callback, core supplies the reporting context automatically. It exposes a [public workflow snapshot API](progress.md); the optional `prompt_toolkit` frontend consumes those snapshots. Functions remain independent of the frontend. Python messages from each workflow worker are routed back to the snapshot's recent messages. Unmanaged backend console bars are disabled when the backend accepts `log_to_console`.

Messages use ordinary Python text output: `print()`, `sys.stdout` / `sys.stderr`
(`write`, `writelines`, and `flush`), and standard `logging.StreamHandler` output.
Console handlers created before a workflow starts are temporarily routed through
the same capture stream; their formatters, filters and levels remain in effect.
File handlers keep writing to their original destinations. Streams and console
handlers are restored when the workflow exits, including on failure. Worker
messages use the existing queue/Dask transport and reach the UI as plain text
through the progress API; plugins do not need to call the UI.

Blank lines, partial writes and multiline tracebacks are preserved, with CRLF
and carriage returns normalized to message boundaries. This captures Python text
streams; native file-descriptor writes and inherited subprocess output need to be
explicitly read from a pipe and forwarded to a captured text stream.

To report measurable work in a custom function, use the tqdm-compatible helper:

```python
from vhrharmonize.io.progress import progress

def process_images(images, progress_callback=None):
    with progress(total=len(images), desc="Processing", unit="images",
                  callback=progress_callback) as bar:
        for image in images:
            process_one(image)
            bar.update(1)
```

`progress` uses tqdm for counters, throttling and rate estimates. When a callback exists, it sends snapshots without rendering; otherwise it is silent by default. Use `disable=False` for an ordinary standalone tqdm display. Functions decorated with `reports_progress` inherit callbacks through nested calls and expose the optional `progress_callback` keyword. Workflow plugin calls also receive this lifecycle reporting automatically. Use `@reports_progress(worker_progress=True)` when the function forwards measurable progress from its workers. Lifecycle-only functions display `(no cb)` beside their step name; counts and core completion tracking still work.

SpectralMatch uses the same callback fields. Its shared image-task runner reports completed images after the parent has committed each result, and forwards worker GDAL progress through queues for local threads/processes or Dask events for remote workers. The application callback stays in the parent process and can be a closure; renderers are never pickled. VHR's own nested footprint and FLAASH tasks use the same transport pattern. Aggregate functions can therefore report both batch progress and work inside individual images. Install the updated SpectralMatch code on the execution host and all Dask workers to enable these callbacks; older versions remain usable but show `(no cb)`.

Callbacks provide only progress the function can measure. Raster loops, Py6S bands, atmosphere sample points and seamline results report increments; GDAL orthorectification reports its fractional callback. Opaque third-party calls and rasterio's overview building display their current phase until the call returns. Callback updates never mark an entire workflow task complete; core does that after output validation.

::: vhrharmonize.io.logging
