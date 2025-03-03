import os
from osgeo import gdal
from helper_functions import windows_to_wsl_path

def orthorectify_image_with_envi(
        input_image_path,
        output_image_path,
        dem_path,
        envi_engine,
        output_image_path_to_delete=None
        ):

    if output_image_path_to_delete:
        try:
            if os.path.exists(output_image_path_to_delete):
                os.remove(output_image_path_to_delete)
                print(f"Deleted existing output image: {output_image_path_to_delete}")
        except Exception as e:
            print(f"Error deleting output image {output_image_path_to_delete}: {e}")


    ortho_params = {
        'INPUT_RASTER': {
            'url': input_image_path,
            'factory': 'URLRaster'
        },
        'DEM_RASTER': {
            'url': dem_path,
            'factory': 'URLRaster'
        },
        'OUTPUT_RASTER_URI': output_image_path,
    }

    try:

        task = envi_engine.task('RPCOrthorectification')
        result = task.execute(ortho_params)

        print(f"Orthorectification complete for {input_image_path}")

    except Exception as e:
        print(f"Error processing {input_image_path}: {e}")

def pansharpen_image_with_envi(
        input_low_res,
        input_high_res,
        output_image_path,
        envi_engine,
        output_image_to_delete=None,
        output_image_masked_path=None
        ):
    # if output_image_to_delete and os.path.exists(output_image_to_delete):
    #     os.remove(output_image_to_delete)
    #     print(f"Deleted existing output: {output_image_to_delete}")

    low_ds, high_ds = gdal.Open(windows_to_wsl_path(input_low_res)), gdal.Open(windows_to_wsl_path(input_high_res))
    if not low_ds or not high_ds:
        print("Error opening input images.")
        return

    high_mask = get_nodata_mask(high_ds)  # Use only the high-resolution image as the mask

    pansharpen_params = {
        'INPUT_LOW_RESOLUTION_RASTER': {'url': input_low_res, 'factory': 'URLRaster'},
        'INPUT_HIGH_RESOLUTION_RASTER': {'url': input_high_res, 'factory': 'URLRaster'},
        'OUTPUT_RASTER_URI': output_image_path,
    }

    print(f"Starting Gram-Schmidt pansharpening for {input_low_res} & {input_high_res}")
    task = envi_engine.task('GramSchmidtPanSharpening')
    # task.execute(pansharpen_params)
    print(f"Pansharpening complete: {output_image_path}")

    if output_image_masked_path:
        ds = gdal.Open(windows_to_wsl_path(output_image_path), gdal.GA_ReadOnly)
        driver = gdal.GetDriverByName('GTiff')
        masked_ds = driver.Create(output_image_masked_path, ds.RasterXSize, ds.RasterYSize, ds.RasterCount, ds.GetRasterBand(1).DataType)
        masked_ds.SetGeoTransform(ds.GetGeoTransform())
        masked_ds.SetProjection(ds.GetProjection())

        for i in range(ds.RasterCount):
            band_data = ds.GetRasterBand(i + 1).ReadAsArray() * high_mask  # Apply only high-resolution mask
            masked_ds.GetRasterBand(i + 1).WriteArray(band_data)
            masked_ds.GetRasterBand(i + 1).SetNoDataValue(ds.GetRasterBand(i + 1).GetNoDataValue())

        masked_ds.FlushCache()
        ds, masked_ds = None, None  # Close datasets
        print(f"Masked output created and saved: {output_image_masked_path}")
    print("done")

def convert_dn_to_radiance_with_envi(
        pan_tif_file,
        pan_radiance_path,
        envi_engine,
        output_image_path_to_delete=None
        ):
    """
Convert DN to radiance using ENVI's RadiometricCalibration tool.
Copies RPC data from input image to output image using GDAL.
Deletes output_image_path_to_delete if provided and exists.
    """
    import os
    from osgeo import gdal

    if output_image_path_to_delete:
        try:
            if os.path.exists(output_image_path_to_delete):
                os.remove(output_image_path_to_delete)
                print(f"Deleted existing output image: {output_image_path_to_delete}")
        except Exception as e:
            print(f"Error deleting output image {output_image_path_to_delete}: {e}")

    # Define parameters for Radiometric Calibration
    radiometric_params = {
        'INPUT_RASTER': {
            'url': pan_tif_file,
            'factory': 'URLRaster'
        },
        'OUTPUT_RASTER_URI': pan_radiance_path
    }

    try:
        print(f"Starting Radiometric Calibration for: {pan_tif_file}")

        # Access the radiometric calibration task
        task = envi_engine.task("RadiometricCalibration")

        # Execute the calibration task with the defined parameters
        task.execute(radiometric_params)

        print(f"Radiometric calibration complete. Output saved to: {pan_radiance_path}")

        # Copy RPC metadata from input image to output image using GDAL
        src_ds = gdal.Open(pan_tif_file, gdal.GA_ReadOnly)
        dst_ds = gdal.Open(pan_radiance_path, gdal.GA_Update)
        if src_ds and dst_ds:
            rpc_metadata = src_ds.GetMetadata("RPC")
            if rpc_metadata:
                dst_ds.SetMetadata(rpc_metadata, "RPC")
                dst_ds.FlushCache()
                print("RPC metadata copied successfully.")
            dst_ds = None
        src_ds = None

    except Exception as e:
        print(f"Error during radiometric calibration: {e}")

def get_nodata_mask(
        dataset
        ):
    band = dataset.GetRasterBand(1)
    nodata_value = band.GetNoDataValue()
    array = band.ReadAsArray()
    return (array == nodata_value) if nodata_value is not None else np.ones_like(array, dtype=bool)