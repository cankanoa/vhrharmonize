from osgeo import gdal, osr

def convert_vertical_datum(
        input_image_path,
        output_image_path,
        input_override_srs=None,
        output_override_srs=None
        ):
    """
    Convert between vertical datums (ellipsoidal and geometric height).

    Parameters:
    - input_image_path: Path to the input raster image.
    - output_image_path: Path to save the output raster.
    - input_override_srs: Override the input spatial reference system (WKT format).
    - output_override_srs: Override the output spatial reference system (WKT format).

    Throws:
    - RuntimeError: If no vertical datum can be determined for the input and no override is provided.
    """
    # Open the input image
    dataset = gdal.Open(input_image_path, gdal.GA_ReadOnly)
    if not dataset:
        raise FileNotFoundError(f"Unable to open input image: {input_image_path}")

    # Get the spatial reference of the input
    input_srs = None
    if input_override_srs:
        input_srs = osr.SpatialReference()
        input_srs.ImportFromWkt(input_override_srs)
    else:
        wkt = dataset.GetProjection()
        if wkt:
            input_srs = osr.SpatialReference()
            input_srs.ImportFromWkt(wkt)
        else:
            raise RuntimeError("Input raster does not have a spatial reference system, and no override was provided.")

    # Check if the input SRS has a vertical datum
    if not input_srs.IsVertical():
        raise RuntimeError("Input spatial reference system does not have a vertical datum.")

    # Handle output SRS
    if output_override_srs:
        output_srs = osr.SpatialReference()
        output_srs.ImportFromWkt(output_override_srs)
    else:
        output_srs = input_srs.Clone()
        print("No output SRS override provided; using the same as input SRS.")

    # Set up the transformation
    transform = osr.CoordinateTransformation(input_srs, output_srs)

    # Create the output raster
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.CreateCopy(output_image_path, dataset, strict=0)

    if not output_dataset:
        raise RuntimeError("Failed to create output dataset.")

    # Apply the transformation to each pixel
    band = dataset.GetRasterBand(1)
    output_band = output_dataset.GetRasterBand(1)
    transform_fn = lambda x, y, z: transform.TransformPoint(x, y, z)

    for y in range(band.YSize):
        scanline = band.ReadAsArray(0, y, band.XSize, 1)
        geotransform = dataset.GetGeoTransform()
        for x in range(band.XSize):
            px = geotransform[0] + x * geotransform[1] + y * geotransform[2]
            py = geotransform[3] + x * geotransform[4] + y * geotransform[5]
            pz = scanline[0, x]
            px_out, py_out, pz_out = transform_fn(px, py, pz)
            scanline[0, x] = pz_out
        output_band.WriteArray(scanline, 0, y)

    # Save and close the datasets
    dataset = None
    output_dataset = None
    print(f"Vertical datum conversion completed. Output saved to: {output_image_path}")


# if __name__ == "__main__":
#     # Paths to the input DEM and output converted DEM
#     dem_path = "/mnt/x/Imagery/Lidar/Big_Island/2018_PuuWaawaa/DEM/2018_2020_bigIsland_DEM_J970216_000_000.tif"
#     output_path = "/mnt/d/outputconvertdem.tif"
#
#     # Define input SRS: WGS84 with EGM96 geoid model
#     input_srs = osr.SpatialReference()
#     input_srs.ImportFromEPSG(4326)  # WGS84 horizontal
#     input_srs.SetAttrValue("VERT_CS", "EGM96")  # Add EGM96 as a vertical reference
#
#     # Define output SRS: WGS84 3D (ellipsoidal height)
#     output_srs = osr.SpatialReference()
#     output_srs.ImportFromEPSG(4979)  # WGS84 3D (latitude, longitude, and ellipsoidal height)
#
#     # Call the conversion function
#     convert_vertical_datum(
#         input_image_path=dem_path,
#         output_image_path=output_path,
#         input_override_srs=input_srs.ExportToWkt(),
#         output_override_srs=output_srs.ExportToWkt()
#     )