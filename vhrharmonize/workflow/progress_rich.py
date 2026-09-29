"""Rich frontend consuming the public progress snapshot, with no scheduler access."""

from dataclasses import dataclass
from threading import RLock

from rich.console import Console, Group
from rich.layout import Layout
from rich.live import Live
from rich.measure import Measurement
from rich.panel import Panel
from rich.progress_bar import ProgressBar
from rich.table import Table
from rich.text import Text

from vhrharmonize.progress import validate_progress_snapshot


def _duration(seconds):
    if seconds is None:
        return "estimating…"
    seconds = max(0, round(seconds))
    hours, seconds = divmod(seconds, 3600)
    minutes, seconds = divmod(seconds, 60)
    return f"~{hours:d}h {minutes:02d}m" if hours else f"~{minutes:d}m {seconds:02d}s"


def _bar(row, width=22):
    total = max(row["all"], 1)
    reused = min(width, round(width * row["reused"] / total))
    complete = min(width - reused, round(width * row["done"] / total))
    bar = Text()
    bar.append("━" * reused, style="grey50")
    bar.append("━" * complete, style="green")
    bar.append("━" * (width - reused - complete), style="grey23")
    return bar


def _eta(row):
    status = row["status"]
    return _duration(row["eta_seconds"]) if status == "running" else "done" if status == "completed" else status


@dataclass
class _StepLabel:
    row: dict

    def __rich_measure__(self, console, options):
        return Measurement(8, len(self.row["name"]) + (8 if not self.row["worker_progress"] else 0))

    def __rich_console__(self, console, options):
        suffix = " (no cb)" if not self.row["worker_progress"] else ""
        label = Text(self.row["name"])
        label.truncate(max(1, options.max_width - len(suffix)), overflow="ellipsis")
        label.append(suffix, style="dim")
        yield label


class RichProgressDisplay:
    """Consume snapshots from callbacks or HPC, then render live or once."""

    def __init__(self, *, console=None):
        self.console = console or Console(stderr=True)
        self.snapshot = None
        self.lock = RLock()
        self.live = None

    def update(self, snapshot):
        data = validate_progress_snapshot(snapshot)
        with self.lock:
            self.snapshot = data

    def table(self, data=None, max_rows=None):
        data = self.snapshot if data is None else data
        table = Table(expand=True, box=None, padding=(0, 1), pad_edge=False)
        bar_width = max(5, min(22, self.console.width - 94))
        table.add_column("Step", ratio=1, min_width=12, overflow="ellipsis", no_wrap=True)
        for name in ("Unused", "Done", "Run", "All"):
            table.add_column(name, justify="right")
        table.add_column("Progress", min_width=bar_width, no_wrap=True)
        table.add_column("Active", justify="right", min_width=6, no_wrap=True)
        table.add_column("ETA", justify="right", max_width=11, no_wrap=True)
        if data is None:
            return table
        rows = data["rows"]
        hidden = 0
        if max_rows is not None and len(rows) > max_rows:
            selected = sorted(rows, key=lambda r: (not r["active"], r["done"] >= r["run"], r["pending"]))[:max(0, max_rows - 1)]
            selected_names = {r["name"] for r in selected}
            hidden = len(rows) - len(selected)
            rows = [r for r in rows if r["name"] in selected_names]
        for row in [data["total"], *rows]:
            table.add_row(
                _StepLabel(row),
                *(f"{row[key]}({row['percentages'][key]:.0f}%)" for key in ("unused", "done", "run", "all")),
                _bar(row, width=bar_width), str(row["active"]), _eta(row),
            )
        if hidden:
            table.add_row(Text(f"… {hidden} more steps", style="dim"), "", "", "", "", "", "", "")
        return table

    def render(self, *, live=None):
        live = self.console.is_terminal if live is None else live
        with self.lock:
            data = self.snapshot
        if data is None:
            return Text("Waiting for progress", style="dim")
        active = data["active"]
        operations = []
        for task in active[:3]:
            stats = task["stats"]
            scene = str(stats.get("scene", task["scene"]))
            scene = scene if len(scene) <= 48 else "…" + scene[-47:]
            operations.append(Text(f"{scene} · {task['step']} · {stats['prefix']}", no_wrap=True, overflow="ellipsis"))
            n, total = stats["n"], stats["total"]
            count = f"{n:g} / {total:g} {stats['unit']}" if total is not None else "working"
            meter = Table.grid(padding=(0, 1))
            meter.add_row(ProgressBar(total=total, completed=n, width=25),
                          Text(f"{count} · {stats['elapsed']:.0f}s · ETA {_duration(task['eta_seconds'])}", style="dim"))
            operations.append(meter)
        if len(active) > 3:
            operations.append(Text(f"+ {len(active) - 3} other active operations", style="dim"))
        if not operations:
            operations = [Text("No active operations", style="dim")]
        active_height = min(len(operations) + 2, 10)
        available = max(6, self.console.height - active_height - 3)
        max_rows = len(data["rows"])
        table = self.table(data)
        if live:
            options = self.console.options.update(width=max(1, self.console.width - 4))
            while True:
                step_height = len(self.console.render_lines(table, options, pad=False)) + 2
                if step_height <= available or max_rows <= 1:
                    break
                max_rows -= 1
                table = self.table(data, max_rows=max_rows)
        else:
            step_height = available
        steps = Panel(table, title="Workflow progress", title_align="left")
        recent = Panel(Group(*(Text(m) for m in data["messages"])) if data["messages"] else Text("Waiting for messages", style="dim"),
                       title="Recent messages", title_align="left")
        active_panel = Panel(Group(*operations), title="Active operation", title_align="left")
        if live:
            layout = Layout()
            layout.split_column(Layout(recent, minimum_size=3), Layout(steps, size=step_height),
                                Layout(active_panel, size=active_height))
            return layout
        return Group(recent, steps, active_panel)

    def __enter__(self):
        if self.console.is_terminal:
            self.live = Live(console=self.console, screen=True, get_renderable=self.render,
                             refresh_per_second=4, redirect_stdout=False, redirect_stderr=False)
            self.live.start()
        return self

    def __exit__(self, *args):
        if self.live is not None:
            self.live.stop()
        if self.snapshot is not None:
            self.console.print(Panel(Group(self.table(), *(Text(m) for m in self.snapshot["messages"])),
                                     title="Workflow failed" if self.snapshot["status"] == "failed" else "Workflow progress"))
