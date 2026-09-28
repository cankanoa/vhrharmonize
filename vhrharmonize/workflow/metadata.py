"""Append completed const/var contexts to explicitly configured JSON destinations."""

import json
import os
from pathlib import Path
from vhrharmonize.io.metadata import write_json


class FinalMetadataWriter:
    def __init__(self, *, delete_first=True):
        self.delete_first = delete_first
        self.written_paths = set()

    def append(self, filename, context):
        # Resolve aliases before remembering first writes or opening destinations.
        filename = os.path.realpath(os.path.abspath(os.path.expanduser(filename)))
        first = filename not in self.written_paths
        entries = []
        if os.path.exists(filename) and not (first and self.delete_first):
            entries = json.loads(Path(filename).read_text())
            if isinstance(entries, dict):
                entries = [entries]
            if not isinstance(entries, list):
                raise ValueError(f"Final metadata must contain a JSON object or array: {filename}")
        entries.append(context)
        # Atomic replacement also clears an old file on the first successful write.
        write_json(filename, entries)
        self.written_paths.add(filename)
