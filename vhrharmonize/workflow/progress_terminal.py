"""Portable terminal frontend consuming progress snapshots, without scheduler access."""

from collections import deque
import os
import sys
from threading import Event, RLock, Thread
from time import monotonic
import warnings

from prompt_toolkit.application import Application, run_in_terminal
from prompt_toolkit.formatted_text import FormattedText, to_plain_text
from prompt_toolkit.input import create_input
from prompt_toolkit.key_binding import KeyBindings
from prompt_toolkit.layout import Layout, Window
from prompt_toolkit.layout.controls import FormattedTextControl
from prompt_toolkit.output import create_output
from prompt_toolkit.shortcuts import print_formatted_text
from prompt_toolkit.utils import get_cwidth

from vhrharmonize.progress import validate_progress_snapshot


COLORS = {"unused": "#4b87b9", "reused": "#956ac3", "done": "#298845", "run": "#858585"}
COUNT_COLUMNS = (("Unused", "unused"), ("Loaded", "reused"), ("Done", "done"), ("Run", "run"), ("All", "all"))
TITLE = "VHRHarmonize Workflow Progress"


def _duration(seconds):
    if seconds is None:
        return "TBD"
    hours, seconds = divmod(max(0, round(seconds)), 3600)
    minutes, seconds = divmod(seconds, 60)
    return f"~{hours}h {minutes:02d}m" if hours else f"~{minutes}m {seconds:02d}s"


def _fit(parts, width, *, right=False):
    """Fit styled cells by terminal character width, including non-ASCII IDs."""
    result, used = [], 0
    for style, text in parts:
        clipped = ""
        for character in str(text).replace("\n", " ").replace("\r", " ").expandtabs(4):
            size = get_cwidth(character)
            if used + size > width:
                result.append((style, clipped))
                padding = [("", " " * max(0, width - used))]
                return padding + result if right else result + padding
            clipped += character
            used += size
        result.append((style, clipped))
        if used >= width:
            break
    padding = [("", " " * max(0, width - used))]
    return padding + result if right else result + padding


def _bar(row, width=22, *, ascii_only=False):
    # Loaded includes cached-but-unused scenes; paint each scene just once.
    counts = (("unused", max(0, row["all"] - row["reused"] - row["run"])),
              ("reused", row["reused"]), ("done", row["done"]),
              ("run", max(0, row["run"] - row["done"])))
    total, units, end, parts = max(1, row["all"]), 0, 0, []
    for key, count in counts:
        units = min(total, units + count)
        boundary = round(width * units / total)
        parts.append((f"fg:{COLORS[key]}", ("#" if ascii_only else "▬") * (boundary - end)))
        end = boundary
    parts.append((f"fg:{COLORS['run']}", ("-" if ascii_only else "▬") * (width - end)))
    return FormattedText(parts)


def _line(cells, widths, *, right=(), gap=2):
    result = []
    for index, (cell, width) in enumerate(zip(cells, widths)):
        if index:
            result.append(("", " " * gap))
        result.extend(_fit(cell, width, right=index in right))
    return result


def _panel(title, lines, width, *, ascii_only=False):
    width = max(4, width)
    horizontal, vertical = ("-", "|") if ascii_only else ("─", "│")
    left, right, bottom_left, bottom_right = ("+",) * 4 if ascii_only else ("╭", "╮", "╰", "╯")
    label = to_plain_text(FormattedText(_fit([("", " " + title + " ")], width - 4))).rstrip()
    result = [("", left + horizontal + label + horizontal * max(0, width - 3 - get_cwidth(label)) + right + "\n")]
    for line in lines:
        if line is None:
            divider_left, divider_right = ("+", "+") if ascii_only else ("├", "┤")
            result.append(("", divider_left + horizontal * (width - 2) + divider_right + "\n"))
        else:
            result += [("", vertical + " "), *_fit(line, width - 4), ("", " " + vertical + "\n")]
    result.append(("", bottom_left + horizontal * (width - 2) + bottom_right + "\n"))
    return result


class TerminalProgressDisplay:
    """Display progress below logs in the terminal's original screen and scrollback.

    The application owns a UI thread only. Processing and application callbacks
    remain on their original threads; redirected output never creates an app.
    """

    def __init__(self, *, stream=None, width=None, ascii_only=None):
        self.stream = stream if stream is not None else sys.stderr
        self.width = width
        self.snapshot = None
        self.lock = RLock()
        self.messages = deque(maxlen=10000)
        self.message_sequence = 0
        self._message_task = None
        self.app = self.thread = None
        self.ready, self.stopping = Event(), Event()
        self.error = None
        self._pinned = FormattedText()
        self._pinned_height = 0
        if ascii_only is None:
            try:
                "╭─│▬…".encode(getattr(self.stream, "encoding", None) or "utf-8")
                ascii_only = False
            except UnicodeEncodeError:
                ascii_only = True
        self.ascii_only = ascii_only

    def update(self, snapshot):
        data = validate_progress_snapshot(snapshot)
        with self.lock:
            previous = self.snapshot
            if previous is None or data["run_id"] != previous["run_id"]:
                self.messages.clear()
                self.message_sequence = 0
                previous = None
            history = data.get("message_history")
            if history is not None:
                count = max(0, history["sequence"] - self.message_sequence)
                incoming = history["messages"][-count:] if count else []
                if count > len(history["messages"]):
                    incoming = ["[Earlier messages are no longer available]", *incoming]
                self.message_sequence = max(self.message_sequence, history["sequence"])
            else:
                # Older saved snapshots have only a short recent-message window.
                old = previous["messages"] if previous else []
                overlap = next((n for n in range(min(len(old), len(data["messages"])), 0, -1)
                                if old[-n:] == data["messages"][:n]), 0)
                incoming = data["messages"][overlap:]
            for message in incoming:
                for line in message.splitlines() or [""]:
                    self.messages.append(line)
            self.snapshot = data
        if self.app is not None:
            self.app.invalidate()

    def table(self, data=None, *, width=None, max_rows=None):
        data = self.snapshot if data is None else data
        if data is None:
            return [[("", "Waiting for progress")]]
        width = width or self.width or 116
        rows = data["rows"]
        hidden = 0
        if max_rows is not None and len(rows) > max_rows:
            chosen = sorted(rows, key=lambda r: (not r["active"], r["done"] >= r["run"], r["pending"]))[:max(0, max_rows - 1)]
            names = {r["name"] for r in chosen}
            hidden = len(rows) - len(chosen)
            rows = [r for r in rows if r["name"] in names]
        rows = [data["total"], *rows]
        counts = [[f"{r[key]}({r['percentages'][key]:.0f}%)" for _, key in COUNT_COLUMNS] for r in rows]
        etas = [_duration(r["eta_seconds"]) if r["status"] in {"running", "waiting"}
                else "done" if r["status"] == "completed" else r["status"] for r in rows]
        gap = 1 if width < 100 else 2
        number_widths = [max(len(name), *(len(c[i]) for c in counts)) for i, (name, _) in enumerate(COUNT_COLUMNS)]
        active_width = max(6, *(len(str(r["active"])) for r in rows))
        eta_width = max(3, *map(len, etas))
        available = width - sum(number_widths) - active_width - eta_width - gap * 8
        labels = [r["name"] + (" (no cb)" if not r["worker_progress"] else "") for r in rows]
        step_width = max(8, min(max(12, *map(get_cwidth, labels)), max(8, available // 2)))
        bar_width = max(8, available - step_width)
        widths = [step_width, *number_widths, active_width, eta_width, bar_width]
        headers = [[("bold", "Step")], *[[(f"bold fg:{COLORS[key]}" if key in COLORS else "bold", name)]
                   for name, key in COUNT_COLUMNS], [("bold", "Active")], [("bold", "ETA")], [("bold", "Progress")]]
        lines = [_line(headers, widths, right=(1, 2, 3, 4, 5, 6, 7), gap=gap)]
        for row, numbers, eta in zip(rows, counts, etas):
            suffix = " (no cb)" if not row["worker_progress"] else ""
            name = to_plain_text(_fit([("", row["name"])], max(0, step_width - len(suffix)))).rstrip()
            label = [("", name + suffix)]
            cells = [label, *[[("", n)] for n in numbers], [("", str(row["active"]))], [("", eta)],
                     _bar(row, bar_width, ascii_only=self.ascii_only)]
            lines.append(_line(cells, widths, right=(1, 2, 3, 4, 5, 6, 7), gap=gap))
        if hidden:
            lines.append([("", f"... {hidden} more steps")])
        return lines

    def operations(self, active, *, width=None, max_rows=None):
        if not active:
            return [[("", "No active operations")]]
        width = width or self.width or 116
        shown = active if max_rows is None else active[:max_rows]
        gap = 1 if width < 100 else 2
        statuses = [str(t["stats"].get("status") or ("starting" if t["stats"]["prefix"] == "starting" else "working")) for t in shown]
        status_width = min(16, max(6, *map(len, statuses)))
        elapsed_width = max(7, *(len(f"{t['stats']['elapsed']:.0f}s") for t in shown))
        eta_width = max(3, *(len(_duration(t["eta_seconds"])) for t in shown))
        remaining = max(3, width - status_width - elapsed_width - eta_width - 5 * gap)
        step_width = max(1, min(max(8, *(get_cwidth(t["step"]) for t in shown)), remaining // 4))
        id_width = max(1, min(max(12, *(get_cwidth(str(t["stats"].get("scene", t["scene"]))) for t in shown)), remaining // 2))
        widths = [step_width, id_width, status_width, elapsed_width, eta_width, max(1, remaining - step_width - id_width)]
        headers = [[("bold", label)] for label in ("Step", "ID", "Status", "Elapsed", "ETA", "Progress")]
        lines = [_line(headers, widths, right=(3, 4), gap=gap)]
        for task, status in zip(shown, statuses):
            stats, size = task["stats"], widths[-1]
            if stats["total"] is None:
                phase = int(monotonic() * 4) % 6
                colors = ["done" if (i - phase) % 6 < 3 else "run" for i in range(size)]
            else:
                done = round(size * min(1, stats["n"] / max(stats["total"], 1)))
                colors = ["done"] * done + ["run"] * (size - done)
            bar = [(f"fg:{COLORS[key]}", "#" if self.ascii_only else "▬") for key in colors]
            cells = [[("", task["step"])], [("", str(stats.get("scene", task["scene"])))],
                     [("", status)], [("", f"{stats['elapsed']:.0f}s")], [("", _duration(task["eta_seconds"]))], bar]
            lines.append(_line(cells, widths, right=(3, 4), gap=gap))
        if len(shown) < len(active):
            lines.append([("", f"+ {len(active) - len(shown)} other active operations")])
        return lines

    def render(self, *, width=None, include_messages=True):
        with self.lock:
            data = self.snapshot
        if data is None:
            return FormattedText([("", "Waiting for progress\n")])
        width = width or self.width or 120
        messages = [[("", m)] for m in data["messages"]] if include_messages else []
        message_lines = [*messages, []] if messages else []
        return FormattedText(
            [part for line in message_lines for part in _fit(line, width) + [("", "\n")]]
            + _panel(TITLE, [*self.table(data, width=width - 4), None,
                             *self.operations(data["active"], width=width - 4)],
                     width, ascii_only=self.ascii_only))

    def print_snapshot(self, *, include_messages=True):
        interactive = getattr(self.stream, "isatty", lambda: False)() and os.environ.get("TERM") != "dumb"
        output = create_output(stdout=self.stream) if interactive else None
        rendered = self.render(width=output.get_size().columns if output is not None and self.width is None else None,
                               include_messages=include_messages)
        encoding = getattr(self.stream, "encoding", None) or "utf-8"
        rendered = FormattedText([(style, text.encode(encoding, errors="replace").decode(encoding))
                                  for style, text in rendered])
        if interactive:
            print_formatted_text(rendered, output=output, end="")
        else:
            self.stream.write(to_plain_text(rendered))
            self.stream.flush()

    def _print_messages(self):
        """Write complete lines once, using prompt_toolkit's normal terminal output."""
        with self.lock:
            if not self.messages:
                return
            text = "\n".join(self.messages) + "\n"
            self.app.print_text(FormattedText([("", text)]))
            self.messages.clear()

    def _messages_printed(self, task):
        if not task.cancelled() and task.exception() is not None:
            self.error = task.exception()
        self.app.invalidate()

    def _after_render(self, app):
        if app.is_done:
            app.output.show_cursor()
            app.output.flush()

    def _prepare_screen(self, app):
        size = app.output.get_size()
        width = max(4, size.columns)
        with self.lock:
            if (self.messages and not self.stopping.is_set() and self.error is None
                    and (self._message_task is None or self._message_task.done())):
                # Run in the app's event-loop context. The library temporarily
                # erases the live box, writes normal terminal lines, and redraws.
                self._message_task = run_in_terminal(self._print_messages)
                self._message_task.add_done_callback(self._messages_printed)
            data = self.snapshot
            if data is None:
                self._pinned, self._pinned_height = FormattedText([("", "Waiting for progress")]), 1
                return
            active_limit = max(1, (size.rows - 10) // 2)
            operations = self.operations(data["active"], width=width - 4, max_rows=active_limit)
            row_limit = max(1, size.rows - len(operations) - 10)
            steps = self.table(data, width=width - 4, max_rows=row_limit)
            if size.rows < 13:
                # Keep the total and first operation visible on short terminals.
                steps = steps[:2]
                operations = operations[:2]
            self._pinned = FormattedText(_panel(TITLE, [*steps, None, *operations], width,
                                               ascii_only=self.ascii_only))
            # No trailing blank line inside the pinned window.
            self._pinned[-1] = (self._pinned[-1][0], self._pinned[-1][1].rstrip("\n"))
            self._pinned_height = len(steps) + len(operations) + 3

    def create_application(self, *, input=None, output=None):
        keys = KeyBindings()

        @keys.add("c-c")
        def interrupt(event):
            from _thread import interrupt_main
            interrupt_main()  # Preserve Ctrl-C processing cancellation on the calling main thread.

        body = Window(FormattedTextControl(lambda: self._pinned), height=lambda: self._pinned_height,
                      dont_extend_height=True, always_hide_cursor=True)
        self.app = Application(layout=Layout(body), key_bindings=keys,
                               full_screen=False, mouse_support=False, erase_when_done=True, refresh_interval=0.25,
                               before_render=self._prepare_screen, after_render=self._after_render,
                               input=input, output=output)
        return self.app

    def __enter__(self):
        if (getattr(self.stream, "isatty", lambda: False)() and sys.stdin.isatty()
                and os.environ.get("TERM") != "dumb"):
            self.create_application(input=create_input(stdin=sys.stdin), output=create_output(stdout=self.stream))

            def ready():
                self.ready.set()
                if self.stopping.is_set():
                    self.app.exit()

            def run():
                try:
                    self.app.run(pre_run=ready, handle_sigint=False, set_exception_handler=False)
                except EOFError:
                    pass  # Disconnected input falls back to the final static report.
                except Exception as exc:
                    self.error = exc
                finally:
                    self.app.output.show_cursor()
                    self.app.output.flush()
                    self.ready.set()
                    self.app.input.close()

            self.thread = Thread(target=run, name="vhr-dashboard", daemon=True)
            self.thread.start()
            self.ready.wait(5)
        return self

    def __exit__(self, *args):
        self.stopping.set()
        if self.thread is not None:
            def stop():
                if self.app.is_running and not self.app.is_done:
                    self.app.exit()
            try:
                if self.app.loop is not None:
                    self.app.loop.call_soon_threadsafe(stop)
            except RuntimeError:
                pass  # The input stream may already have closed the event loop.
            self.thread.join(5)
        if self.error is not None:
            warnings.warn(f"Interactive dashboard unavailable: {self.error}", RuntimeWarning, stacklevel=2)
        if self.app is not None:
            self._print_messages()  # Flush final lines that arrived just before shutdown.
        if self.snapshot is not None:
            self.print_snapshot(include_messages=self.app is None)
