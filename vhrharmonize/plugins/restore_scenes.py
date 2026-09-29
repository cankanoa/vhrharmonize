"""Restore a materialized scene snapshot, used by generated HPC recipes."""

from copy import deepcopy
import os
from vhrharmonize.workflow.values import assign, lookup, path, remap_paths
from .base import FunctionPlugin
from vhrharmonize.io.progress import progress, reports_progress


@reports_progress(worker_progress=True)
def restore_scenes(records: list[dict], directory_locations: dict | None = None) -> dict:
    """Return saved scene variables and their declared file bookkeeping."""
    scenes = []
    for original in progress(records, desc="Restoring scenes", unit="scenes"):
        record = deepcopy(original)
        expanded_paths = {p: os.path.expanduser(p) for p in record.get("file_paths", [])}
        variables = remap_paths(record["context"]["var"], expanded_paths)
        for field in record.get("path_fields", []):
            assign(variables, field, path(lookup(variables, field), base_dir="."))
        record["id"] = expanded_paths.get(record["id"], record["id"])
        record["source_paths"] = [os.path.expanduser(p) for p in record["source_paths"]]
        scenes.append(
            {
                **variables,
                "restored_record": {key: record[key] for key in ("id", "source_paths")},
            }
        )
    return {"scenes": scenes}


class RestoreScenes(FunctionPlugin):
    target = "vhrharmonize.plugins.restore_scenes:restore_scenes"
    scene_records_return = "scenes"
    scene_id_return = "restored_record.id"
    source_file_protection_paths_return = "restored_record.source_paths"

    def run(self, *, params, shared):
        locations = params.get("directory_locations", {})
        self.temporary_directory_context_paths = tuple(locations.get("temp_dir", ()))
        self.output_directory_context_paths = tuple(locations.get("output_dir", ()))
        return super().run(params=params, shared=shared)
