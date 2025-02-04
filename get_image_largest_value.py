import os
from osgeo import gdal, ogr, osr
import numpy as np

def get_image_largest_value(input_image_path, mask=None, override_mask_crs_epsg=None):
    """
    Get the largest value in a raster image, optionally masked by a shapefile.

    Args:
        input_image_path (str): Path to the raster image.
        mask (str, optional): Path to a geographic mask (shapefile). Defaults to None.
        override_mask_crs_epsg (int, optional): EPSG code to override the mask's coordinate system. Defaults to None.

    Returns:
        float: The largest value in the raster (masked if mask is provided).
    """
    # Open the raster dataset
    dataset = gdal.Open(input_image_path, gdal.GA_ReadOnly)
    if not dataset:
        raise FileNotFoundError(f"Raster file not found: {input_image_path}")

    # Get raster geotransform and nodata value
    geo_transform = dataset.GetGeoTransform()
    projection = dataset.GetProjection()
    x_size = dataset.RasterXSize
    y_size = dataset.RasterYSize

    band = dataset.GetRasterBand(1)
    nodata_value = band.GetNoDataValue()

    # Read the raster data as an array
    raster_array = band.ReadAsArray()
    if nodata_value is not None:
        raster_mask = (raster_array != nodata_value)
    else:
        raster_mask = np.ones_like(raster_array, dtype=bool)

    # If a mask is provided, check for overlaps in valid data
    mask_array = None
    if mask:
        mask_ds = ogr.Open(mask)
        if not mask_ds:
            raise FileNotFoundError(f"Mask file not found: {mask}")

        mask_layer = mask_ds.GetLayer()

        # Override the mask's CRS if requested
        if override_mask_crs_epsg:
            target_srs = osr.SpatialReference()
            target_srs.ImportFromEPSG(override_mask_crs_epsg)

            # Create a memory layer with the corrected CRS
            driver = ogr.GetDriverByName("Memory")
            mem_ds = driver.CreateDataSource("in_memory")
            mem_layer = mem_ds.CreateLayer("mask", srs=target_srs, geom_type=mask_layer.GetGeomType())

            # Copy features from the original mask layer
            for feature in mask_layer:
                mem_layer.CreateFeature(feature.Clone())

            mask_layer = mem_layer  # Replace original mask layer with corrected one

        # Create an in-memory raster for the mask
        driver = gdal.GetDriverByName("MEM")
        mask_raster = driver.Create("", x_size, y_size, 1, gdal.GDT_Byte)
        mask_raster.SetGeoTransform(geo_transform)
        mask_raster.SetProjection(projection)

        # Rasterize the shapefile to the mask raster
        gdal.RasterizeLayer(mask_raster, [1], mask_layer, burn_values=[1])
        mask_array = mask_raster.ReadAsArray()

        # Check if the mask has true values where the raster has no valid data
        overlap_check = np.logical_and(mask_array == 1, ~raster_mask)
        if np.any(overlap_check):
            print("Warning: The provided mask contains areas where there is no valid data in the input image.")

        # Clean up
        mask_ds = None
        mask_raster = None
        if override_mask_crs_epsg:
            mem_ds = None

    # Loop through all bands and find the largest value
    largest_value = None
    for band_index in range(1, dataset.RasterCount + 1):
        band = dataset.GetRasterBand(band_index)
        band_array = band.ReadAsArray()

        if nodata_value is not None:
            band_array = np.ma.masked_where(band_array == nodata_value, band_array)

        if mask_array is not None:
            # Apply the mask to the band data
            band_array = np.ma.masked_where(mask_array == 0, band_array)

        # Find the maximum value in the current band
        band_max = np.max(band_array)
        largest_value = max(largest_value, band_max) if largest_value is not None else band_max

    # Clean up
    dataset = None

    return float(largest_value)

# Example Usage
# largest_value = get_image_largest_value("input_raster.tif", mask="mask.shp", override_mask_crs_epsg=4326)
# print("Largest value in the raster:", largest_value)
