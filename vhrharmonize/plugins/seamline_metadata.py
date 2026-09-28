"""Build seamline metadata GeoPackages for downstream seamline ranking."""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from types import SimpleNamespace
from typing import Any, Dict, List, Mapping

import geopandas as gpd
from osgeo import gdal, ogr, osr
from shapely import make_valid
from shapely.affinity import affine_transform
from shapely.wkt import loads as wkt_loads

from vhrharmonize.workflow.concurrency import (
    _make_dask_client,
    _resolve_concurrent_processing,
    _resolve_concurrent_processing_backend,
)
from .base import FunctionPlugin
from vhrharmonize.io.logging import _log, _log_image_completed, _log_image_start, _log_step_start
from vhrharmonize.io.metadata import materialize_geometry


def _valid_data_polygon_from_image(path: str, *, eight_connected: bool = True) -> object:
    """Extract the largest valid-data polygon from a raster mask using GDAL."""
    dataset = gdal.Open(path, gdal.GA_ReadOnly)
    if dataset is None:
        raise RuntimeError(f"Cannot open {path}")
    geotransform = dataset.GetGeoTransform()
    band = dataset.GetRasterBand(1)
    if band is None:
        raise RuntimeError(f"Raster has no band 1: {path}")
    mask = band.GetMaskBand()
    if mask is None:
        raise RuntimeError(f"Raster has no valid mask band: {path}")

    datasource = ogr.GetDriverByName("MEM").CreateDataSource("mem")
    layer = datasource.CreateLayer("valid_data", geom_type=ogr.wkbPolygon)
    layer.CreateField(ogr.FieldDefn("val", ogr.OFTInteger))

    options = ["8CONNECTED=8"] if eight_connected else None
    if gdal.Polygonize(mask, None, layer, 0, options=options, callback=None) != gdal.CE_None:
        raise RuntimeError(f"Cannot polygonize valid-data mask: {path}")

    layer.ResetReading()
    best_area = -1.0
    best_geometry = None
    for feature in layer:
        if feature.GetField("val") != 255:
            continue
        geometry = feature.GetGeometryRef()
        if geometry is not None and geometry.GetArea() > best_area:
            best_area = geometry.GetArea()
            best_geometry = geometry.Clone()

    if best_geometry is None:
        raise ValueError(f"No valid-data polygon found: {path}")

    pixel_polygon = wkt_loads(best_geometry.ExportToWkt())
    dataset = None
    geometry = affine_transform(
        pixel_polygon,
        (
            geotransform[1],
            geotransform[2],
            geotransform[4],
            geotransform[5],
            geotransform[0],
            geotransform[3],
        ),
    )
    if not geometry.is_valid:
        geometry = make_valid(geometry, method="structure", keep_collapsed=False)
    if (
        geometry.is_empty
        or not geometry.is_valid
        or geometry.geom_type not in {"Polygon", "MultiPolygon"}
    ):
        raise ValueError(f"Cannot repair valid-data polygon: {path}")
    return geometry


def _calculate_seamline_metadata_geometry(
    image_path: str,
    source_metadata: Mapping[str, Any],
    footprint_source: str,
    calculate_bounds_eight_connected: bool,
    epsg: int,
    scene_basename: str,
    output_path: str,
    log_to_console: bool,
) -> tuple[str, object]:
    """Open inputs in the worker and return only the image basename and geometry."""
    _log_image_start(
        scene_basename,
        [image_path],
        [output_path],
        enabled=log_to_console,
        step="seamline_metadata",
    )
    if footprint_source == "calculate_bounds":
        geometry = _valid_data_polygon_from_image(
            image_path,
            eight_connected=calculate_bounds_eight_connected,
        )
        geometry = _project_image_geometry(image_path, geometry, epsg)
    elif footprint_source == "metadata":
        geometry = materialize_geometry(source_metadata["geometry"], epsg=epsg)
    else:
        raise ValueError(f"Unsupported seamline metadata footprint source: {footprint_source}")
    return os.path.basename(image_path), geometry


def _project_image_geometry(image_path, geometry, epsg):
    """The declared output CRS may differ from the imported raster's CRS."""
    from pyproj import CRS, Transformer
    from shapely.ops import transform

    dataset = gdal.Open(image_path, gdal.GA_ReadOnly)
    projection = dataset.GetProjectionRef() if dataset is not None else None
    dataset = None
    if not projection:
        raise ValueError(f"Raster has no CRS for its footprint: {image_path}")
    source_crs = CRS.from_wkt(projection)
    target_crs = CRS.from_epsg(epsg)
    if source_crs != target_crs:
        geometry = transform(
            Transformer.from_crs(source_crs, target_crs, always_xy=True).transform, geometry
        )
    return geometry


def _iter_seamline_metadata_results(
    tasks: list[tuple],
    *,
    worker_count: int,
    backend: str,
    dask_scheduler_file: str | None,
    dask_scheduler_address: str | None,
):
    """Yield footprints as workers finish, using the shared workflow backends."""
    if backend == "dask":
        from dask.distributed import as_completed as dask_as_completed

        client = _make_dask_client(
            SimpleNamespace(
                dask_scheduler_file=dask_scheduler_file,
                dask_scheduler_address=dask_scheduler_address,
            )
        )
        futures = []
        try:
            futures = [
                client.submit(_calculate_seamline_metadata_geometry, *task) for task in tasks
            ]
            for future in dask_as_completed(futures):
                yield future.result()
        except BaseException:
            if futures:
                client.cancel(futures)
            raise
        finally:
            client.close()
    elif worker_count <= 1 or len(tasks) <= 1:
        for task in tasks:
            yield _calculate_seamline_metadata_geometry(*task)
    else:
        with ProcessPoolExecutor(max_workers=min(worker_count, len(tasks))) as executor:
            futures = [
                executor.submit(_calculate_seamline_metadata_geometry, *task) for task in tasks
            ]
            try:
                for future in as_completed(futures):
                    yield future.result()
            finally:
                for future in futures:
                    future.cancel()


def _open_seamline_metadata_writer(output_path, layer, epsg, records, *, append):
    """Open a parent-only writer and establish fields before appending results."""
    if append:
        datasource = ogr.Open(output_path, update=1)
        output_layer = datasource.GetLayerByName(layer) if datasource is not None else None
    else:
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        if os.path.exists(output_path):
            os.remove(output_path)
        datasource = ogr.GetDriverByName("GPKG").CreateDataSource(output_path)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg)
        output_layer = (
            datasource.CreateLayer(layer, srs=srs, geom_type=ogr.wkbUnknown)
            if datasource is not None
            else None
        )
    if output_layer is None:
        raise RuntimeError(
            f"Cannot open seamline metadata layer for writing: {output_path}:{layer}"
        )

    # Inspect all metadata, including nulls, so completion order cannot change
    # the type chosen for a column (e.g. a null first result followed by a float).
    field_names = dict.fromkeys(key for record in records.values() for key in record)
    for name in field_names:
        if output_layer.GetLayerDefn().GetFieldIndex(name) >= 0:
            continue
        values = [record[name] for record in records.values() if record.get(name) is not None]
        if values and all(isinstance(value, (int, bool)) for value in values):
            field_type = ogr.OFTInteger64
        elif values and all(isinstance(value, (int, float)) for value in values):
            field_type = ogr.OFTReal
        else:
            field_type = ogr.OFTString
        if output_layer.CreateField(ogr.FieldDefn(name, field_type)) != ogr.OGRERR_NONE:
            raise RuntimeError(f"Cannot create seamline metadata field: {name}")
    return datasource, output_layer


def _write_seamline_metadata_record(datasource, output_layer, record, geometry):
    """Commit one returned footprint before counting it as processed."""
    if datasource.StartTransaction() != ogr.OGRERR_NONE:
        raise RuntimeError("Cannot start seamline metadata transaction.")
    try:
        feature = ogr.Feature(output_layer.GetLayerDefn())
        for name, value in record.items():
            if value is not None:
                feature.SetField(name, int(value) if isinstance(value, bool) else value)
        if feature.SetGeometry(ogr.CreateGeometryFromWkb(geometry.wkb)) != ogr.OGRERR_NONE:
            raise RuntimeError(f"Cannot set footprint geometry: {record['image_basename']}")
        if output_layer.CreateFeature(feature) != ogr.OGRERR_NONE:
            raise RuntimeError(f"Cannot write footprint: {record['image_basename']}")
        feature = None
        if datasource.CommitTransaction() != ogr.OGRERR_NONE:
            raise RuntimeError(f"Cannot commit footprint: {record['image_basename']}")
        datasource.FlushCache()
    except BaseException:
        datasource.RollbackTransaction()
        raise


def write_seamline_metadata_gpkg(
    image_paths: List[str],
    output_path: str,
    *,
    metadata_records: List[dict],
    layer: str = "seamline_metadata",
    image_field_name: str = "image_path",
    footprint_source: str = "calculate_bounds",
    calculate_bounds_eight_connected: bool = True,
    epsg: int = 4326,
    scene_total: int | None = None,
    reuse: bool = False,
    log_to_console: bool = False,
    concurrent_processing: int | str = 1,
    concurrent_processing_backend: str = "process_pool",
    dask_scheduler_file: str | None = None,
    dask_scheduler_address: str | None = None,
) -> str:
    """Calculate missing footprints in workers and commit each result in the parent.

    Matching image_basename values are already processed and never submitted.
    Completed records remain available for reuse if another worker fails.


    """

    _log_step_start("seamline_metadata", enabled=log_to_console)
    worker_count = _resolve_concurrent_processing(concurrent_processing)
    backend = _resolve_concurrent_processing_backend(concurrent_processing_backend)
    if backend == "dask" and worker_count != 1:
        raise ValueError(
            "concurrent_processing must be 1 when concurrent_processing_backend is 'dask'."
        )
    if footprint_source not in {"calculate_bounds", "metadata"}:
        raise ValueError(f"Unsupported seamline metadata footprint source: {footprint_source}")

    existing_image_basenames = set()
    append = reuse and os.path.exists(output_path)
    if append:
        existing_gdf = gpd.read_file(output_path, layer=layer)
        if "image_basename" not in existing_gdf.columns:
            raise ValueError("Existing seamline metadata is missing image_basename.")
        if image_field_name not in existing_gdf.columns:
            raise ValueError(f"Existing seamline metadata is missing {image_field_name}.")
        if existing_gdf.crs is None or existing_gdf.crs.to_epsg() != epsg:
            raise ValueError("Existing seamline metadata CRS does not match epsg.")
        existing_image_basenames.update(existing_gdf["image_basename"].dropna())
        del existing_gdf

    if len(image_paths) != len(metadata_records):
        raise ValueError("image_paths and metadata_records must have the same length")
    images = {
        os.path.basename(path): (path, metadata)
        for path, metadata in zip(image_paths, metadata_records)
    }
    if len(images) != len(set(image_paths)):
        raise ValueError("Seamline images need unique basenames")
    processed_image_basenames = existing_image_basenames.intersection(images)
    total = scene_total if scene_total is not None else len(images)
    records: Dict[str, Dict[str, Any]] = {}
    tasks = []
    for image_basename, (image_path, metadata) in images.items():
        if image_basename in processed_image_basenames:
            continue
        record = {
            key: value
            for key, value in metadata.items()
            if isinstance(value, (str, int, float, bool)) or value is None
        }
        record.update(
            {
                image_field_name: image_path,
                "image_basename": image_basename,
                "scene_basename": metadata.get("scene_id", image_basename),
            }
        )
        records[image_basename] = record
        tasks.append(
            (
                image_path,
                metadata,
                footprint_source,
                calculate_bounds_eight_connected,
                epsg,
                record["scene_basename"],
                output_path,
                log_to_console,
            )
        )

    if not tasks:
        if append:
            return output_path
        raise ValueError("No scene outputs were available for seamline metadata.")

    results = _iter_seamline_metadata_results(
        tasks,
        worker_count=worker_count,
        backend=backend,
        dask_scheduler_file=dask_scheduler_file,
        dask_scheduler_address=dask_scheduler_address,
    )
    datasource = output_layer = None
    try:
        for completed, (image_basename, geometry) in enumerate(results, start=1):
            if datasource is None:
                datasource, output_layer = _open_seamline_metadata_writer(
                    output_path,
                    layer,
                    epsg,
                    records,
                    append=append,
                )
            record = records[image_basename]
            _write_seamline_metadata_record(datasource, output_layer, record, geometry)
            processed_image_basenames.add(image_basename)
            _log_image_completed(
                record["scene_basename"],
                completed,
                total,
                processing_total=len(tasks),
                enabled=log_to_console,
                step="seamline_metadata",
            )
    finally:
        output_layer = None
        datasource = None
        results.close()
    return output_path


__all__ = [
    "SeamlineMetadata",
    "write_seamline_metadata_gpkg",
]


class SeamlineMetadata(FunctionPlugin):
    input_dependency_paths = frozenset({"image_paths"})
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset({"output_path"})
    output_dependency_paths = output_path_resolution_paths
    output_target_paths = output_path_resolution_paths
    output_parent_creation_paths = output_path_resolution_paths
    output_collision_check_paths = output_path_resolution_paths
    # The function owns record-level reuse and incomplete-output recovery.
    output_temporary_cleanup_paths = output_path_resolution_paths
    output_hpc_staging_paths = output_path_resolution_paths
    output_hpc_download_paths = output_path_resolution_paths
    output_context_checkpoint_paths = frozenset({"output_path"})

    scope = "aggregate"
    target = "vhrharmonize.plugins.seamline_metadata:write_seamline_metadata_gpkg"
