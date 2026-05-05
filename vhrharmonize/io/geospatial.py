import os
import warnings
from typing import Iterable, Optional
import rasterio
import rasterio.mask
import geopandas as gpd
import numpy as np

from osgeo import ogr, osr
from shapely.geometry import mapping


def calculate_raster_overviews(
    input_image_path: str,
    overview_scales: Optional[Iterable[int]],
    *,
    resampling: str = "nearest",
    ) -> str:
    """Build internal raster overviews.
    Args:
        input_image_path: Input raster path.
        overview_scales: Requested overview decimation factors.
        resampling: Rasterio overview resampling method name.
    Returns:
        Input raster path after overview creation.
    """
    factors = []
    seen = set()
    for value in overview_scales or []:
        factor = int(value)
        if factor <= 1 or factor in seen:
            continue
        seen.add(factor)
        factors.append(factor)

    if not factors:
        return input_image_path

    with rasterio.open(input_image_path, "r+") as dataset:
        dataset.build_overviews(factors, rasterio.enums.Resampling[resampling])
        dataset.update_tags(ns="rio_overview", resampling=resampling)
    return input_image_path


def shp_to_gpkg(
    input_shp_path: str,
    output_gpkg_path: str,
    override_projection_epsg: Optional[int] = None
    ) -> None:
    """Convert a shapefile to a GeoPackage.
    Args:
        input_shp_path: Input shapefile path.
        output_gpkg_path: Output GeoPackage path.
        override_projection_epsg: Optional CRS EPSG override to assign without reprojection.
    Returns:
        None.
    """

    # Open the input shapefile (read-only)
    shp_driver = ogr.GetDriverByName("ESRI Shapefile")
    in_ds = shp_driver.Open(input_shp_path, 0)  # 0 = read-only
    if not in_ds:
        raise FileNotFoundError(f"Could not open Shapefile: {input_shp_path}")

    in_layer = in_ds.GetLayer()
    if not in_layer:
        raise RuntimeError("Could not get layer from Shapefile.")

    # Determine the layer's current spatial reference
    in_srs = in_layer.GetSpatialRef()

    # If user wants to override the CRS, just assign that EPSG
    if override_projection_epsg is not None:
        out_srs = osr.SpatialReference()
        out_srs.ImportFromEPSG(override_projection_epsg)
    else:
        out_srs = in_srs  # keep the original (which could be None)

    # Prepare output: if the GPKG already exists, remove it
    gpkg_driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(output_gpkg_path):
        os.remove(output_gpkg_path)

    if not os.path.exists(os.path.dirname(output_gpkg_path)):
        os.makedirs(os.path.dirname(output_gpkg_path))

    out_ds = gpkg_driver.CreateDataSource(output_gpkg_path)

    # Create the output layer with the (possibly overridden) CRS
    # Note: if out_srs is None, it will simply have "unknown" CRS in the GeoPackage
    out_layer = out_ds.CreateLayer(
        name="layer",
        srs=out_srs,
        geom_type=in_layer.GetGeomType()
    )

    # Copy fields from the input layer
    in_layer_defn = in_layer.GetLayerDefn()
    for i in range(in_layer_defn.GetFieldCount()):
        field_defn = in_layer_defn.GetFieldDefn(i)
        out_layer.CreateField(field_defn)

    # Copy features from input to output
    in_layer.ResetReading()
    for in_feature in in_layer:
        out_feature = ogr.Feature(out_layer.GetLayerDefn())

        # Copy geometry by cloning (no reprojection)
        geom = in_feature.GetGeometryRef()
        if geom is not None:
            out_feature.SetGeometry(geom.Clone())

        # Copy field attributes
        for i in range(in_layer_defn.GetFieldCount()):
            out_feature.SetField(
                in_layer_defn.GetFieldDefn(i).GetNameRef(),
                in_feature.GetField(i)
            )

        out_layer.CreateFeature(out_feature)
        out_feature = None

    # Cleanup
    out_ds = None
    in_ds = None

    if override_projection_epsg:
        print(f"Assigned EPSG:{override_projection_epsg} to output layer {output_gpkg_path}")
    else:
        print(f"No CRS override. Output saved to '{output_gpkg_path}'.")


def get_image_percentile_value(
    input_image_path: str,
    percentile: float = 50.0,
    mask: Optional[str] = None,
    override_mask_crs_epsg: Optional[int] = None,
    ) -> float:
    """Compute a raster percentile value.
    Args:
        input_image_path: Input raster path.
        percentile: Requested percentile in the inclusive range 0 to 100.
        mask: Optional vector mask path used to limit sampled pixels.
        override_mask_crs_epsg: Optional CRS EPSG override for the mask.
    Returns:
        Percentile value from valid raster pixels.
    """
    if percentile < 0 or percentile > 100:
        raise ValueError("percentile must be between 0 and 100.")

    with rasterio.open(input_image_path) as src:
        nodata = src.nodata
        collected = []

        if mask:
            mask_gdf = gpd.read_file(mask)

            if mask_gdf.crs is None:
                if override_mask_crs_epsg:
                    mask_gdf.set_crs(epsg=override_mask_crs_epsg, inplace=True)
                else:
                    raise ValueError(
                        f"Mask file '{mask}' has no CRS. You must specify 'override_mask_crs_epsg'."
                    )
            elif override_mask_crs_epsg:
                mask_gdf = mask_gdf.to_crs(epsg=override_mask_crs_epsg)

            geometries = [mapping(geom) for geom in mask_gdf.geometry]
        else:
            geometries = None

        for i in range(1, src.count + 1):
            if geometries:
                band_array, _ = rasterio.mask.mask(src, geometries, indexes=i, filled=False)
            else:
                band_array = src.read(i, masked=True)

            if nodata is not None:
                band_array = np.ma.masked_equal(band_array, nodata)

            if band_array.mask.all():
                continue

            valid = band_array.compressed() if np.ma.isMaskedArray(band_array) else band_array.ravel()
            if valid.size:
                collected.append(valid)

    if not collected:
        warnings.warn(
            f"Warning: Percentile value in file '{input_image_path}' could not be found. "
            "This could be because the mask is outside of the bounds of this image."
        )
        return float("nan")

    values = np.concatenate(collected)
    return float(np.percentile(values, percentile))
