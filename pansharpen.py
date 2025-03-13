import orthority as oty
from osgeo import gdal

def pansharpen_image(
        input_low_resolution_path,
        input_high_resolution_path,
        output_image_path,
        change_nodata_value = None,
        ):
    pan_sharp = oty.PanSharpen(input_high_resolution_path, input_low_resolution_path)
    pan_sharp.process(output_image_path, write_mask=False)

    if change_nodata_value is not None:
        _change_nodata_value(output_image_path, change_nodata_value)

def _change_nodata_value(input_image_path, new_nodata_value):
    """
    Changes the NoData value of a raster in place and updates pixels with the original NoData value to the new one.

    :param input_image_path: Path to the input raster file.
    :param new_nodata_value: The new NoData value to be set.
    """
    # Open the input raster in update mode
    dataset = gdal.Open(input_image_path, gdal.GA_Update)
    if dataset is None:
        raise FileNotFoundError(f"Cannot open {input_image_path}")

    # Read the first band
    band = dataset.GetRasterBand(1)
    old_nodata_value = band.GetNoDataValue()

    if old_nodata_value is None:
        print("No existing NoData value found. Setting a new one.")

    # Read the raster data into a NumPy array
    raster_data = band.ReadAsArray()

    # Replace old NoData values with the new NoData value
    if old_nodata_value is not None:
        raster_data[raster_data == old_nodata_value] = new_nodata_value

    # Write modified data back to the same raster
    band.WriteArray(raster_data)
    band.SetNoDataValue(new_nodata_value)  # Update NoData metadata
    band.FlushCache()

    # Close dataset
    dataset = None

    print(f"NoData value updated in-place: {old_nodata_value} → {new_nodata_value}")