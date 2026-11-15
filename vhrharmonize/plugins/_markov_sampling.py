"""Bounded raster sampling for SpectralMatch's Markov node workers.

Warped VRT blocks can expand to almost an entire VHR source image at coarse
GSDs. Use the VRT only to define the sampling grid; materialize each requested
window with the chunked GDAL warper instead of reading its virtual bands.
"""

from collections import OrderedDict
from contextlib import contextmanager
from importlib import import_module
from threading import RLock

import numpy as np
from osgeo import gdal


_native = import_module("spectralmatch.seamline.markov_triangles")
_native_solve_node = _native._solve_node
_dispatch_lock = RLock()
_sampling_lock = RLock()


class BoundedRasterSampler:
    """Sample the native map-aligned grids with a 32 MiB warp/cache budget each."""

    window_size = 256
    warp_memory_mb = 32
    cache_bytes = 32 * 1024 * 1024

    def __init__(self, crs):
        self.crs = crs
        self.sources = {}
        self.grids = {}
        self.tiles = OrderedDict()
        self.tile_bytes = 0

    def close(self):
        self.tiles.clear()
        self.tile_bytes = 0
        self.grids.clear()
        self.sources.clear()

    def _grid(self, path, gsd):
        if path not in self.sources:
            self.sources[path] = gdal.Open(path)
        key = (path, gsd)
        if key not in self.grids:
            self.grids[key] = self.sources[path] if gsd is None else gdal.Warp(
                "", self.sources[path], format="VRT", dstSRS=self.crs,
                xRes=gsd, yRes=gsd, targetAlignedPixels=True,
                resampleAlg="average", outputType=gdal.GDT_Float64,
                dstNodata=float("nan"),
            )
        return self.grids[key]

    def _tile(self, path, gsd, bands, tx, ty, grid):
        key = (path, gsd, bands, tx, ty)
        if key in self.tiles:
            self.tiles.move_to_end(key)
            return self.tiles[key]
        xoff, yoff = tx * self.window_size, ty * self.window_size
        width = min(self.window_size, grid.RasterXSize - xoff)
        height = min(self.window_size, grid.RasterYSize - yoff)
        dataset = grid
        if gsd is not None:
            transform = grid.GetGeoTransform()
            left, top = gdal.ApplyGeoTransform(transform, xoff, yoff)
            right, bottom = gdal.ApplyGeoTransform(transform, xoff + width, yoff + height)
            # Keep all source bands so GDAL's multiband NoData semantics remain
            # identical. MEM output is small; source reads are chunked by -wm.
            dataset = gdal.Warp(
                "", self.sources[path], format="MEM", dstSRS=self.crs,
                outputBounds=(left, bottom, right, top), width=width, height=height,
                resampleAlg="average", outputType=gdal.GDT_Float64,
                dstNodata=float("nan"), warpMemoryLimit=self.warp_memory_mb,
            )
            xoff = yoff = 0
        values = np.empty((height, width, len(bands)), dtype=np.float64)
        for channel, number in enumerate(bands):
            band = dataset.GetRasterBand(number)
            data = band.ReadAsArray(xoff, yoff, width, height)
            valid = band.GetMaskBand().ReadAsArray(xoff, yoff, width, height) != 0
            values[:, :, channel] = np.where(valid, data, np.nan)
        while self.tiles and self.tile_bytes + values.nbytes > self.cache_bytes:
            _, old = self.tiles.popitem(last=False)
            self.tile_bytes -= old.nbytes
        if values.nbytes <= self.cache_bytes:
            self.tiles[key] = values
            self.tile_bytes += values.nbytes
        return values

    def sample(self, path, coordinates, bands, gsd, native_fallback=False):
        """Preserve band order, masks, grid-edge tolerances and native fallback."""
        bands = tuple(bands)
        grid = self._grid(path, gsd)
        inverse = gdal.InvGeoTransform(grid.GetGeoTransform())
        x, y = coordinates.T
        col = inverse[0] + inverse[1] * x + inverse[2] * y
        row = inverse[3] + inverse[4] * x + inverse[5] * y
        for values, size in ((col, grid.RasterXSize), (row, grid.RasterYSize)):
            values[np.isclose(values, size, rtol=0, atol=1e-7)] = size - 1e-7
            values[np.isclose(values, 0, rtol=0, atol=1e-7)] = 0
        cols, rows = np.floor(col).astype(int), np.floor(row).astype(int)
        inside = np.flatnonzero(
            (cols >= 0) & (cols < grid.RasterXSize)
            & (rows >= 0) & (rows < grid.RasterYSize)
        )
        result = np.full((len(coordinates), len(bands)), np.nan)
        windows, groups = np.unique(
            np.column_stack((cols[inside] // self.window_size, rows[inside] // self.window_size)),
            axis=0, return_inverse=True,
        )
        order = np.argsort(groups, kind="stable")
        starts = np.r_[0, np.cumsum(np.bincount(groups, minlength=len(windows)))]
        for index, (tx, ty) in enumerate(windows):
            selected = inside[order[starts[index]:starts[index + 1]]]
            values = self._tile(path, gsd, bands, int(tx), int(ty), grid)
            result[selected] = values[
                rows[selected] - ty * self.window_size,
                cols[selected] - tx * self.window_size,
            ]
        missing = ~np.isfinite(result).all(axis=1)
        if native_fallback and missing.any():
            native = self.sample(path, coordinates[missing], bands, None)
            result[missing] = np.where(np.isfinite(result[missing]), result[missing], native)
        return result


def _solve_node(*args, **kwargs):
    """Install the sampler in the executing process, including spawn/Dask workers."""
    with _sampling_lock:
        previous = _native._RasterSampler
        _native._RasterSampler = BoundedRasterSampler
        try:
            return _native_solve_node(*args, **kwargs)
        finally:
            _native._RasterSampler = previous


@contextmanager
def bounded_markov_sampling():
    """Dispatch a picklable worker wrapper without modifying the dependency on disk."""
    with _dispatch_lock:
        previous = _native._solve_node
        _native._solve_node = _solve_node
        try:
            yield
        finally:
            _native._solve_node = previous
