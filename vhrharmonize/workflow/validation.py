"""Bounded file validation and a cache scoped to one planning pass."""

import os
from concurrent.futures import as_completed

from vhrharmonize.io.logging import _log
from vhrharmonize.io.validation import _existing_output_failures, _is_gdal_raster_path
from .concurrency import image_executor, _resolve_concurrent_processing


class ValidationPass:
    """Check each path/mode once; never retain results across execution."""

    def __init__(self, controls):
        self.controls = controls
        self.results = {}

    @staticmethod
    def key(filename, enabled):
        return os.path.normcase(os.path.abspath(filename)), bool(enabled)

    def prepare(self, requests):
        pending = {self.key(filename, enabled) for filename, enabled in requests} - self.results.keys()
        # Existence checks stay local and never open raster/vector contents.
        for filename, enabled in sorted(pending):
            exists = os.path.exists(filename)
            supported = _is_gdal_raster_path(filename) or os.path.splitext(filename)[1].lower() in {".json", ".gpkg"}
            if not enabled or not exists or not supported:
                self.results[filename, enabled] = None if exists else "missing"
        pending -= self.results.keys()
        if not pending:
            return
        options = dict(
            check_validity=True,
            validity_raster_grid=self.controls["validity_raster_grid"],
            validity_raster_check_stop=self.controls["validity_raster_check_stop"],
            tile_threads=self.controls["tile_threads"],
            log_to_console=False, step="core:validation",
        )
        total = len(pending)
        enabled = self.controls["log_to_console"]
        _log(f"Start | 0% (0/{total} files)", enabled=enabled, step="core:validation")
        done = 0
        last_percent = 0

        def completed(key, failures):
            nonlocal done, last_percent
            self.results[key] = failures.get(key[0])
            done += 1
            percent = 100 * done // total
            if percent != last_percent:
                _log(f"{percent}% ({done}/{total} files)", enabled=enabled, step="core:validation")
                last_percent = percent

        workers = _resolve_concurrent_processing(self.controls["concurrent_processing"])
        backend = self.controls["concurrent_processing_backend"]
        if workers == 1 and backend != "dask":
            for key in sorted(pending):
                completed(key, _existing_output_failures([key[0]], **options))
        else:
            with image_executor(workers, backend, self.controls.get("dask_scheduler")) as executor:
                futures = {executor.submit(_existing_output_failures, [key[0]], **options): key
                           for key in sorted(pending)}
                for future in as_completed(futures):
                    completed(futures[future], future.result())

    def failures(self, filenames, enabled):
        filenames = list(filenames)
        self.prepare((filename, enabled) for filename in filenames)
        return {filename: reason for filename in filenames
                if (reason := self.results[self.key(filename, enabled)]) is not None}
