"""Build seamline metadata GeoPackages for downstream seamline ranking."""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from types import SimpleNamespace
from typing import Any, Dict, List, Mapping

import geopandas as gpd
from osgeo import gdal, ogr, osr
from shapely.affinity import affine_transform
from shapely.wkt import loads as wkt_loads

from vhrharmonize.preprocess.concurrency import (
    _make_dask_client, _resolve_concurrent_processing,
    _resolve_concurrent_processing_backend)
from vhrharmonize.preprocess.helpers import (_log, _log_image_completed,
                                             _log_image_start, _log_step_start)
from vhrharmonize.providers.standardized import materialize_scene_bounds



def _standardized_metadata_fields(metadata: object) -> Dict[str, Any]:
    """Return scalar standardized metadata fields using their existing names."""
    values = metadata.to_dict()
    return {
        key: value
        for key, value in values.items()
        if isinstance(value, (str, int, float, bool)) or value is None
    }


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
    return affine_transform(
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
    _log_image_start(scene_basename, [image_path], [output_path], enabled=log_to_console, step="seamline_metadata")
    if footprint_source == "calculate_bounds":
        geometry = _valid_data_polygon_from_image(
            image_path, eight_connected=calculate_bounds_eight_connected,
        )
    elif footprint_source == "package_bounds":
        geometry = materialize_scene_bounds(source_metadata, epsg=epsg)
    else:
        raise ValueError(f"Unsupported seamline metadata footprint source: {footprint_source}")
    return os.path.basename(image_path), geometry


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

        client = _make_dask_client(SimpleNamespace(
            dask_scheduler_file=dask_scheduler_file,
            dask_scheduler_address=dask_scheduler_address,
        ))
        futures = []
        try:
            futures = [client.submit(_calculate_seamline_metadata_geometry, *task) for task in tasks]
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
            futures = [executor.submit(_calculate_seamline_metadata_geometry, *task) for task in tasks]
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
        output_layer = datasource.CreateLayer(layer, srs=srs, geom_type=ogr.wkbUnknown) if datasource is not None else None
    if output_layer is None:
        raise RuntimeError(f"Cannot open seamline metadata layer for writing: {output_path}:{layer}")

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
    states: List[object],
    output_path: str,
    *,
    layer: str,
    image_field_name: str,
    footprint_source: str,
    calculate_bounds_eight_connected: bool,
    epsg: int,
    scene_total: int | None = None,
    run_from_existing_check_validity: bool = False,
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
        raise ValueError("concurrent_processing must be 1 when concurrent_processing_backend is 'dask'.")
    if footprint_source not in {"calculate_bounds", "package_bounds"}:
        raise ValueError(f"Unsupported seamline metadata footprint source: {footprint_source}")

    existing_image_basenames = set()
    append = run_from_existing_check_validity and os.path.exists(output_path)
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

    states_by_image_basename = {
        os.path.basename(state.current_files[0]): state
        for state in states if state.current_files
    }
    processed_image_basenames = existing_image_basenames.intersection(states_by_image_basename)
    total = scene_total if scene_total is not None else len(states_by_image_basename)
    if processed_image_basenames:
        _log(f"Already processed {len(processed_image_basenames)}/{total}", enabled=log_to_console, step="seamline_metadata")

    records: Dict[str, Dict[str, Any]] = {}
    tasks = []
    for image_basename, state in states_by_image_basename.items():
        if image_basename in processed_image_basenames:
            continue
        image_path = state.current_files[0]
        mul_image = state.scene.mul_image
        if mul_image is None or mul_image.standardized_metadata is None:
            raise ValueError(f"WorldView scene is missing standardized metadata: {state.scene.primary_basename}")
        record = {
            image_field_name: image_path,
            "image_basename": image_basename,
            "scene_basename": state.scene.primary_basename,
            "scene_id": state.scene.scene_id,
            "catalog_id": state.scene.catalog_id,
            "mul_basename": mul_image.basename,
            "mul_imd_file": mul_image.imd_file,
        }
        record.update(_standardized_metadata_fields(mul_image.standardized_metadata))
        records[image_basename] = record
        source_metadata = mul_image.standardized_metadata.source_metadata if footprint_source == "package_bounds" else {}
        tasks.append((
            image_path, source_metadata, footprint_source,
            calculate_bounds_eight_connected, epsg, state.scene.primary_basename,
            output_path, log_to_console,
        ))

    if not tasks:
        if append:
            return output_path
        raise ValueError("No scene outputs were available for seamline metadata.")

    results = _iter_seamline_metadata_results(
        tasks, worker_count=worker_count, backend=backend,
        dask_scheduler_file=dask_scheduler_file, dask_scheduler_address=dask_scheduler_address,
    )
    datasource = output_layer = None
    try:
        for image_basename, geometry in results:
            if datasource is None:
                datasource, output_layer = _open_seamline_metadata_writer(
                    output_path, layer, epsg, records, append=append,
                )
            record = records[image_basename]
            _write_seamline_metadata_record(datasource, output_layer, record, geometry)
            processed_image_basenames.add(image_basename)
            _log_image_completed(record["scene_basename"], len(processed_image_basenames), total, enabled=log_to_console, step="seamline_metadata")
    finally:
        output_layer = None
        datasource = None
        results.close()
    return output_path


__all__ = [
    "write_seamline_metadata_gpkg",
]
