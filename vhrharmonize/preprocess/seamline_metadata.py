"""Build seamline metadata GeoPackages for downstream seamline ranking."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import geopandas as gpd
from osgeo import gdal, ogr
from shapely.affinity import affine_transform
from shapely.wkt import loads as wkt_loads


def _standardized_metadata_fields(metadata: object) -> Dict[str, Any]:
    """Return scalar standardized metadata fields using their existing names."""
    values = metadata.to_dict()
    return {
        key: value
        for key, value in values.items()
        if isinstance(value, (str, int, float, bool)) or value is None
    }


def _scene_footprint_geometry(shp_file: str, *, epsg: int) -> object:
    """Read and merge a scene footprint geometry."""
    gdf = gpd.read_file(shp_file)
    if gdf.empty:
        raise ValueError(f"Empty scene footprint shapefile: {shp_file}")
    if gdf.crs is None:
        raise ValueError(f"Scene footprint shapefile has no CRS: {shp_file}")
    gdf = gdf.to_crs(epsg=epsg)
    geometries = gdf.geometry[gdf.geometry.notnull()]
    if geometries.empty:
        raise ValueError(f"Scene footprint shapefile has no geometry: {shp_file}")
    if hasattr(geometries, "union_all"):
        return geometries.union_all()
    return geometries.unary_union


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
    gdal.Polygonize(mask, None, layer, 0, options=options, callback=None)

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


def write_seamline_metadata_gpkg(
    states: List[object],
    output_path: str,
    *,
    layer: str,
    image_field_name: str,
    footprint_source: str,
    calculate_bounds_eight_connected: bool,
    epsg: int,
) -> str:
    """Write seamline footprint and IMD metadata polygons."""
    records: List[Dict[str, Any]] = []
    for state in states:
        if not state.current_files:
            continue
        mul_image = state.scene.mul_image
        if mul_image is None or mul_image.shp_file is None:
            raise ValueError(f"WorldView scene is missing multispectral shapefile: {state.scene.primary_basename}")
        if mul_image.standardized_metadata is None:
            raise ValueError(f"WorldView scene is missing standardized metadata: {state.scene.primary_basename}")

        image_path = state.current_files[0]
        if footprint_source == "calculate_bounds":
            geometry = _valid_data_polygon_from_image(
                image_path,
                eight_connected=calculate_bounds_eight_connected,
            )
        elif footprint_source == "package_bounds":
            geometry = _scene_footprint_geometry(mul_image.shp_file, epsg=epsg)
        else:
            raise ValueError(f"Unsupported seamline metadata footprint source: {footprint_source}")
        record: Dict[str, Any] = {
            image_field_name: image_path,
            "image_basename": os.path.basename(image_path),
            "scene_basename": state.scene.primary_basename,
            "scene_id": state.scene.scene_id,
            "catalog_id": state.scene.catalog_id,
            "mul_basename": mul_image.basename,
            "mul_imd_file": mul_image.imd_file,
            "mul_shp_file": mul_image.shp_file,
            "geometry": geometry,
        }
        record.update(_standardized_metadata_fields(mul_image.standardized_metadata))
        records.append(record)

    if not records:
        raise ValueError("No scene outputs were available for seamline metadata.")

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    if os.path.exists(output_path):
        os.remove(output_path)
    gdf = gpd.GeoDataFrame(records, geometry="geometry", crs=f"EPSG:{epsg}")
    gdf.to_file(output_path, layer=layer, driver="GPKG")
    return output_path


__all__ = ["write_seamline_metadata_gpkg"]
