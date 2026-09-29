"""Copy a file and explicitly configured companions to a workflow location."""

from pathlib import Path
import shutil
from .base import FunctionPlugin
from vhrharmonize.io.progress import progress, reports_progress


@reports_progress(worker_progress=True)
def copy_file(input_path: str, output_path: str, companions=None):
    """Copy input and companions to the explicitly supplied output path."""
    if not output_path:
        raise ValueError("output_path must be an explicit output path")
    output_path = str(output_path)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_path, output_path)
    for source in progress(companions or [], desc="Copying companion files", unit="files"):
        source = Path(source)
        # Companion extensions follow the copied raster stem.
        destination = Path(output_path).with_suffix(source.suffix)
        if destination != source:
            shutil.copy2(source, destination)
    return output_path


class FileSource(FunctionPlugin):
    input_dependency_paths = frozenset({"companions", "input_path"})
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset({"output_path"})
    output_dependency_paths = output_path_resolution_paths
    output_target_paths = output_path_resolution_paths
    output_parent_creation_paths = output_path_resolution_paths
    output_collision_check_paths = output_path_resolution_paths
    output_reuse_paths = output_path_resolution_paths
    output_validation_paths = output_path_resolution_paths
    output_invalid_removal_paths = output_path_resolution_paths
    output_overview_calculation_paths = frozenset({"output_path"})
    output_temporary_cleanup_paths = output_path_resolution_paths
    output_hpc_staging_paths = output_path_resolution_paths
    output_hpc_download_paths = output_path_resolution_paths
    output_context_checkpoint_paths = frozenset({"output_path"})

    target = "vhrharmonize.plugins.file_source:copy_file"
