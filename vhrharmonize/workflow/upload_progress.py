"""Grouped rsync progress data and a frontend sharing the workflow terminal UI."""

from copy import deepcopy
import os
import re
from time import monotonic

from .progress_terminal import COLORS, TerminalProgressDisplay, _bar, _duration, _line


RSYNC_FORMAT = "VHR:%i:%l:%n"
_ITEM = re.compile(r"^VHR:([^:]+):(\d+):(.*)$")
_PROGRESS = re.compile(r"^\s*([\d,]+)\s+(\d+|nan)%", re.IGNORECASE)


def _bytes(value):
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            return f"{value:.1f}{unit}" if unit != "B" else f"{value:.0f}B"
        value /= 1024


def _filename(value):
    # Rsync escapes control characters and, depending on locale, UTF-8 bytes.
    raw = re.sub(rb"\\#([0-7]{3})", lambda m: bytes([int(m[1], 8)]), os.fsencode(value))
    return os.fsdecode(raw).removeprefix("./")


class UploadProgress:
    """Track each destination once; snapshots contain group totals, never filenames.

    Bytes measure logical file data processed by rsync, not SSH wire traffic.
    Files already current are accounted separately after a successful batch.
    """

    def __init__(self, items, *, groups=(), path_mappings=None, callbacks=()):
        self.callbacks = callbacks
        self.started = monotonic()
        self.finished = None
        self.last_publish = 0
        self.status = "checking"
        self.rows, self.files = {}, {}
        self.batch, self.current = [], None
        labels = {os.path.abspath(filename): (group["step"], group["variable"])
                  for group in groups for filename in group["files"]}
        roots = sorted((path_mappings or {}).items(), key=lambda item: len(item[1]), reverse=True)
        for section, local, remote in items:
            variable = next((name for name, root in roots if remote == root.rstrip("/")
                             or remote.startswith(root.rstrip("/") + "/")), "paths")
            label = labels.get(os.path.abspath(local),
                               ("inputs" if section == "uploaded_input_paths" else "controls", variable))
            if os.path.isdir(local):
                entries = ((os.path.join(directory, filename),
                            remote.rstrip("/") + "/" + os.path.relpath(os.path.join(directory, filename), local).replace(os.sep, "/"))
                           for directory, _, filenames in os.walk(local) for filename in filenames)
            else:
                entries = [(local, remote)]
            for source, destination in entries:
                if destination in self.files:
                    if os.path.abspath(source) != self.files[destination]["source"]:
                        raise ValueError(f"Multiple local files map to one remote path: {destination}")
                    continue
                size = os.path.getsize(source)
                row = self.rows.setdefault(label, dict(step=label[0], variable=label[1], files=0,
                    done=0, current=0, total_bytes=0, transferred_bytes=0, current_bytes=0,
                    started=None, elapsed_seconds=0.0, status="queued"))
                row["files"] += 1
                row["total_bytes"] += size
                self.files[destination] = dict(source=os.path.abspath(source), size=size, row=label,
                                               n=0, done=False, changed=False)

    def begin(self, remote_root, items):
        self.remote_root = remote_root.rstrip("/")
        destinations = {item[2].rstrip("/") for item in items}
        directories = tuple(item[2].rstrip("/") + "/" for item in items if os.path.isdir(item[1]))
        self.batch = [name for name in self.files if name in destinations or name.startswith(directories)]
        self.current = None
        for name in self.batch:
            row = self.rows[self.files[name]["row"]]
            row["status"] = "checking"
            if row["started"] is None:
                row["started"] = monotonic()
        self.publish(force=True)

    def _advance(self, entry, n):
        n = min(entry["size"], max(entry["n"], n))
        self.rows[entry["row"]]["transferred_bytes"] += n - entry["n"]
        entry["n"] = n

    def _complete(self, entry):
        if entry["done"]:
            return
        row = self.rows[entry["row"]]
        if entry["changed"]:
            self._advance(entry, entry["size"])
        else:
            row["current"] += 1
            row["current_bytes"] += entry["size"]
        row["done"] += 1
        entry["done"] = True

    def consume(self, line):
        """Consume item markers and CR/LF progress records without printing them."""
        match = _ITEM.match(line)
        if match:
            flags, _, filename = match.groups()
            self.current = self.files.get(self.remote_root + "/" + _filename(filename))
            if self.current is not None and len(flags) > 1 and flags[1] == "f":
                self.current["changed"] = flags[0] in "<>ch"
                row = self.rows[self.current["row"]]
                row["status"] = "uploading" if self.current["changed"] else "checking"
                if self.current["changed"]:
                    self.status = "uploading"
                if not self.current["changed"]:
                    self._complete(self.current)
            else:
                self.current = None
        elif self.current is not None and (match := _PROGRESS.match(line)):
            self._advance(self.current, int(match[1].replace(",", "")))
            if match[2] == "100" or self.current["size"] == 0:
                self._complete(self.current)
        self.publish()

    def finish(self, *, error=None):
        for name in self.batch:
            entry = self.files[name]
            row = self.rows[entry["row"]]
            if error is None:
                self._complete(entry)
                row["status"] = "done" if row["done"] == row["files"] else "queued"
            else:
                row["status"] = "cancelled" if isinstance(error, KeyboardInterrupt) else "failed"
            row["elapsed_seconds"] = monotonic() - row["started"]
        self.status = ("cancelled" if isinstance(error, KeyboardInterrupt) else "failed") if error else (
            "done" if all(row["done"] == row["files"] for row in self.rows.values()) else "uploading")
        if self.status in {"done", "failed", "cancelled"}:
            self.finished = monotonic()
        self.publish(force=True)

    def snapshot(self):
        now = monotonic()
        rows = []
        for original in self.rows.values():
            row = deepcopy(original)
            elapsed = now - row["started"] if row["started"] is not None and row["status"] in {"checking", "uploading"} else row["elapsed_seconds"]
            row["elapsed_seconds"] = elapsed
            row["rate"] = row["transferred_bytes"] / elapsed if elapsed > 0 and row["transferred_bytes"] else None
            remaining = row["total_bytes"] - row["transferred_bytes"] - row["current_bytes"]
            row["eta_seconds"] = 0 if row["status"] == "done" else remaining / row["rate"] if row["rate"] else None
            row.pop("started")
            rows.append(row)
        total = {key: sum(row[key] for row in rows)
                 for key in ("files", "done", "current", "total_bytes", "transferred_bytes", "current_bytes")}
        elapsed = max(0, (self.finished if self.finished is not None else now) - self.started)
        rate = total["transferred_bytes"] / elapsed if elapsed and total["transferred_bytes"] else None
        total.update(step="total", variable="", status=self.status, elapsed_seconds=elapsed, rate=rate,
                     eta_seconds=0 if self.status == "done" else
                     (total["total_bytes"] - total["transferred_bytes"] - total["current_bytes"]) / rate if rate else None)
        return dict(version=1, kind="upload", status=self.status, total=total, rows=rows, active=[], messages=[])

    def publish(self, *, force=False):
        if not force and monotonic() - self.last_publish < 0.2:
            return
        self.last_publish = monotonic()
        snapshot = self.snapshot()
        for callback in self.callbacks:
            callback(deepcopy(snapshot))


class TerminalUploadDisplay(TerminalProgressDisplay):
    title = "VHRHarmonize HPC Upload"

    def __init__(self, **kwargs):
        super().__init__(read_input=False, **kwargs)

    def update(self, snapshot):
        with self.lock:
            self.snapshot = deepcopy(snapshot)
        if self.app is not None:
            self.app.invalidate()

    def table(self, data=None, *, width=None, max_rows=None):
        data = data or self.snapshot
        if data is None:
            return [[("", "Checking uploads")]]
        width = width or self.width or 116
        rows = data["rows"]
        hidden = 0
        if max_rows is not None and len(rows) > max_rows:
            ordered = sorted(rows, key=lambda row: row["status"] not in {"checking", "uploading"})
            rows = ordered[:max(0, max_rows - 1)]
            hidden = len(data["rows"]) - len(rows)
        rows = [data["total"], *rows]
        values = [[row["step"], row["variable"], f"{row['done']}/{row['files']}",
                   f"{_bytes(row['transferred_bytes'] + row['current_bytes'])}/{_bytes(row['total_bytes'])}",
                   _bytes(row["rate"]) + "/s" if row["rate"] else "TBD", row["status"],
                   "done" if row["status"] == "done" else _duration(row["eta_seconds"])] for row in rows]
        headers = ["Step", "From", "Files", "Bytes", "Speed", "Status", "ETA"]
        # Keep labels and the rightmost bar visible on narrower terminals.
        columns = [0, 1, 2, *([3] if width >= 100 else []), *([4] if width >= 85 else []), 5, *([6] if width >= 70 else [])]
        widths = [max(len(headers[i]), *(len(value[i]) for value in values)) for i in columns]
        widths[0] = min(widths[0], max(4, width // 7))
        widths[1] = min(widths[1], max(4, width // 5))
        bar_width = max(4, width - sum(widths) - len(columns))
        widths.append(bar_width)
        cells = [[("bold" + (f" fg:{COLORS['done']}" if i in {2, 3} else ""), headers[i])] for i in columns]
        lines = [_line([*cells, [("bold", "Progress")]], widths, gap=1)]
        for row, value in zip(rows, values):
            total = row["total_bytes"] or row["files"]
            reused = row["current_bytes"] if row["total_bytes"] else row["current"]
            done = row["transferred_bytes"] if row["total_bytes"] else row["done"] - row["current"]
            bar = _bar(dict(all=total, reused=reused, run=total - reused, done=done),
                       bar_width, ascii_only=self.ascii_only)
            lines.append(_line([*[[("", value[i])] for i in columns], bar], widths, gap=1))
        if hidden:
            lines.append([("", f"... {hidden} more steps")])
        return lines

    def operations(self, active, *, width=None, max_rows=None):
        total = self.snapshot["total"]
        return [[("", f"Elapsed {total['elapsed_seconds']:.0f}s  "),
                 (f"fg:{COLORS['done']}", f"Transferred {_bytes(total['transferred_bytes'])}  "),
                 (f"fg:{COLORS['reused']}", f"Already current {total['current']} files ({_bytes(total['current_bytes'])})")]]
