"""Bounded raster sampling for SpectralMatch's Markov node workers.

Warped VRT blocks can expand to almost an entire VHR source image at coarse
GSDs. Use the VRT only to define the sampling grid; materialize tight windows
around nearby requested pixels with the chunked GDAL warper. Grow cached
windows by their missing strips, retaining the same grid and average values.
"""

from collections import OrderedDict
from contextlib import contextmanager
from importlib import import_module
from threading import RLock

import os

import numpy as np
from osgeo import gdal


_native = import_module("spectralmatch.seamline.markov_triangles")
_native_solve_node = _native._solve_node
_dispatch_lock = RLock()
_sampling_lock = RLock()
# Arrays only: never retain live GDAL handles across node processes/forks.
_worker_tiles = OrderedDict()
_worker_tile_bytes = 0
_worker_cache_limit = 32 * 1024 * 1024
_worker_cache_lock = RLock()


class BoundedRasterSampler:
    """Sample the native map-aligned grids with a 32 MiB warp/cache budget each."""

    window_size = 256
    coarse_window_size = 16
    warp_memory_mb = 32
    cache_bytes = 32 * 1024 * 1024

    def __init__(self, crs):
        self.crs = crs
        self.sources = {}
        self.grids = {}
        self.tiles = OrderedDict()
        self.tile_bytes = 0
        self.identities = {}

    def close(self):
        global _worker_tile_bytes
        with _worker_cache_lock:
            for key, value in self.tiles.items():
                identity = self.identities.get(key[0])
                if identity is None:
                    continue
                shared_key = (self.crs, identity, key)
                old = _worker_tiles.pop(shared_key, None)
                if old is not None:
                    _worker_tile_bytes -= old[2].nbytes
                size = value[2].nbytes
                while _worker_tiles and _worker_tile_bytes + size > _worker_cache_limit:
                    _, removed = _worker_tiles.popitem(last=False)
                    _worker_tile_bytes -= removed[2].nbytes
                if size <= _worker_cache_limit:
                    _worker_tiles[shared_key] = value
                    _worker_tile_bytes += size
        self.tiles.clear()
        self.tile_bytes = 0
        self.grids.clear()
        self.sources.clear()
        self.identities.clear()

    def _grid(self, path, gsd):
        if path not in self.sources:
            self.sources[path] = gdal.Open(path)
            try:
                # Include mask/metadata sidecars which can change validity or georeferencing.
                files = self.sources[path].GetFileList() or [path]
                self.identities[path] = tuple((name, st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns)
                                              for name in sorted(files) for st in [os.stat(name)])
            except OSError:
                self.identities[path] = None  # Non-filesystem sources are not shared.

        key = (path, gsd)
        if key not in self.grids:
            self.grids[key] = self.sources[path] if gsd is None else gdal.Warp(
                "", self.sources[path], format="VRT", dstSRS=self.crs,
                xRes=gsd, yRes=gsd, targetAlignedPixels=True,
                resampleAlg="average", outputType=gdal.GDT_Float64,
                dstNodata=float("nan"),
            )
        return self.grids[key]

    def _read_window(self, path, gsd, bands, window, grid):
        """Read exactly this grid-aligned window; retain all-band warp semantics."""
        xoff, yoff, width, height = window
        dataset = grid
        if gsd is not None:
            transform = grid.GetGeoTransform()
            left, top = gdal.ApplyGeoTransform(transform, xoff, yoff)
            right, bottom = gdal.ApplyGeoTransform(transform, xoff + width, yoff + height)
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
        return values

    def _tile(self, path, gsd, bands, tx, ty, grid, requested):
        key = (path, gsd, bands, tx, ty)
        old = self.tiles.get(key)
        if old is None and self.identities.get(path) is not None:
            with _worker_cache_lock:
                shared_key = (self.crs, self.identities[path], key)
                old = _worker_tiles.get(shared_key)
                if old is not None and old[2].nbytes > self.cache_bytes:
                    old = None
                if old is not None:
                    _worker_tiles.move_to_end(shared_key)
                    self._remember(key, old)
        xoff, yoff, width, height = requested
        if old is not None:
            ox, oy, values = old
            oh, ow = values.shape[:2]
            if ox <= xoff and oy <= yoff and xoff+width <= ox+ow and yoff+height <= oy+oh:
                self.tiles.move_to_end(key)
                return old
            # Grow a cached window using only its missing strips. Nearby
            # Laplacian requests do not re-average pixels already in the cache.
            right, bottom = max(xoff+width, ox+ow), max(yoff+height, oy+oh)
            xoff, yoff = min(xoff, ox), min(yoff, oy)
            width, height = right-xoff, bottom-yoff
            expanded = np.empty((height, width, len(bands)), dtype=np.float64)
            expanded[oy-yoff:oy-yoff+oh, ox-xoff:ox-xoff+ow] = values
            missing = [(xoff, yoff, width, oy-yoff),
                       (xoff, oy+oh, width, bottom-oy-oh),
                       (xoff, oy, ox-xoff, oh),
                       (ox+ow, oy, right-ox-ow, oh)]
            for x, y, w, h in missing:
                if w and h:
                    expanded[y-yoff:y-yoff+h, x-xoff:x-xoff+w] = self._read_window(
                        path, gsd, bands, (x, y, w, h), grid)
            self.tile_bytes -= values.nbytes
            del self.tiles[key]
            values = expanded
        else:
            values = self._read_window(path, gsd, bands, requested, grid)
        result = (xoff, yoff, values)
        self._remember(key, result)
        return result

    def _remember(self, key, value):
        old = self.tiles.pop(key, None)
        if old is not None:
            self.tile_bytes -= old[2].nbytes
        size = value[2].nbytes
        while self.tiles and self.tile_bytes + size > self.cache_bytes:
            _, removed = self.tiles.popitem(last=False)
            self.tile_bytes -= removed[2].nbytes
        if size <= self.cache_bytes:
            self.tiles[key] = value
            self.tile_bytes += size

    def _sampling_window_size(self, path, gsd):
        if gsd is None:
            return self.window_size
        transform = self.sources[path].GetGeoTransform()
        native = max(np.hypot(transform[1], transform[4]),
                     np.hypot(transform[2], transform[5]))
        # Near-native edge costs need many adjacent samples. Tiny bins multiply
        # warp setup there; coarse averages still use tightly bounded windows.
        return self.window_size if gsd <= 4 * native else self.coarse_window_size

    def _prefetch(self, path, gsd, bands, grid, requests):
        """Merge nearby *new* coarse requests, bounding unused pixels and memory."""
        if gsd is None:
            return
        pending = []
        for tx, ty, window in requests:
            key = (path, gsd, bands, tx, ty)
            shared_key = (self.crs, self.identities.get(path), key)
            with _worker_cache_lock:
                reused = shared_key in _worker_tiles
            if key not in self.tiles and not reused:
                pending.append((key, window))
        # Avoid evicting prefetched values before sample() can consume them.
        if sum(w*h*len(bands)*8 for _,(_,_,w,h) in pending) > self.cache_bytes:
            return
        groups = []
        max_pixels = self._sampling_window_size(path, gsd)**2
        for key, window in sorted(pending, key=lambda item: (item[1][1], item[1][0])):
            x,y,w,h = window
            if groups:
                box, area, members = groups[-1]
                ox,oy,ow,oh = box
                left,top,right,bottom = min(x,ox),min(y,oy),max(x+w,ox+ow),max(y+h,oy+oh)
                combined = (right-left)*(bottom-top)
                if combined <= max_pixels and combined <= 1.25*(area+w*h):
                    groups[-1] = ((left,top,right-left,bottom-top),area+w*h,members+[(key,window)])
                    continue
            groups.append((window,w*h,[(key,window)]))
        for (x,y,w,h), _, members in groups:
            if len(members) < 2:
                continue  # Existing on-demand reader handles individual requests.
            values = self._read_window(path,gsd,bands,(x,y,w,h),grid)
            for key,(ox,oy,ow,oh) in members:
                self._remember(key,(ox,oy,values[oy-y:oy-y+oh,ox-x:ox-x+ow].copy()))

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
        size = self._sampling_window_size(path, gsd)
        windows, groups = np.unique(
            np.column_stack((cols[inside] // size, rows[inside] // size)),
            axis=0, return_inverse=True,
        )
        order = np.argsort(groups, kind="stable")
        starts = np.r_[0, np.cumsum(np.bincount(groups, minlength=len(windows)))]
        requests = []
        for index, (tx, ty) in enumerate(windows):
            selected = inside[order[starts[index]:starts[index + 1]]]
            if gsd is None:
                xoff, yoff = int(tx)*size, int(ty)*size
                width = min(size, grid.RasterXSize-xoff)
                height = min(size, grid.RasterYSize-yoff)
            else:
                xoff, yoff = int(cols[selected].min()), int(rows[selected].min())
                width = int(cols[selected].max())-xoff+1
                height = int(rows[selected].max())-yoff+1
            requests.append((int(tx), int(ty), (xoff, yoff, width, height), selected))
        self._prefetch(path, gsd, bands, grid, [(tx,ty,win) for tx,ty,win,_ in requests])
        for tx, ty, window, selected in requests:
            xoff, yoff, values = self._tile(path, gsd, bands, tx, ty, grid, window)
            result[selected] = values[
                rows[selected] - yoff,
                cols[selected] - xoff,
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
