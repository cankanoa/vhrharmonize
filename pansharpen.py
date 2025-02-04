import os
from osgeo import gdal
import tempfile

# def MAIN_pansharpen_image(input_low_resolution_path, input_high_resolution_path, output_image_path):
#     """
#     Pansharpen an image using GDAL after aligning the extents of both images.
#
#     Args:
# input_low_resolution_path (str): Path to the low-resolution multispectral image.
# input_high_resolution_path (str): Path to the high-resolution panchromatic image.
# output_image_path (str): Path to save the pansharpened image.
#     """
#     # Ensure the output directory exists
#     output_dir = os.path.dirname(output_image_path)
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)
#
#     # Open input datasets
#     low_res_ds = gdal.Open(input_low_resolution_path, gdal.GA_ReadOnly)
#     high_res_ds = gdal.Open(input_high_resolution_path, gdal.GA_ReadOnly)
#
#     if not low_res_ds or not high_res_ds:
#         raise FileNotFoundError("One or both input files could not be opened.")
#
#     # Get extents of both datasets
#     low_res_gt = low_res_ds.GetGeoTransform()
#     high_res_gt = high_res_ds.GetGeoTransform()
#
#     # Compute the bounds of both datasets
#     low_res_bounds = (
#         low_res_gt[0],  # minX
#         low_res_gt[3] + low_res_gt[5] * low_res_ds.RasterYSize,  # minY
#         low_res_gt[0] + low_res_gt[1] * low_res_ds.RasterXSize,  # maxX
#         low_res_gt[3],  # maxY
#     )
#     high_res_bounds = (
#         high_res_gt[0],  # minX
#         high_res_gt[3] + high_res_gt[5] * high_res_ds.RasterYSize,  # minY
#         high_res_gt[0] + high_res_gt[1] * high_res_ds.RasterXSize,  # maxX
#         high_res_gt[3],  # maxY
#     )
#
#     # Calculate the smallest bounds
#     common_bounds = (
#         max(low_res_bounds[0], high_res_bounds[0]),  # minX
#         max(low_res_bounds[1], high_res_bounds[1]),  # minY
#         min(low_res_bounds[2], high_res_bounds[2]),  # maxX
#         min(low_res_bounds[3], high_res_bounds[3]),  # maxY
#     )
#
#     # Clip both datasets to the common extent
#     with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as low_res_clipped_file:
#         low_res_clipped_path = low_res_clipped_file.name
#         gdal.Warp(
#             destNameOrDestDS=low_res_clipped_path,
#             srcDSOrSrcDSTab=low_res_ds,
#             outputBounds=common_bounds,
#             format="GTiff"
#         )
#
#     with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as high_res_clipped_file:
#         high_res_clipped_path = high_res_clipped_file.name
#         gdal.Warp(
#             destNameOrDestDS=high_res_clipped_path,
#             srcDSOrSrcDSTab=high_res_ds,
#             outputBounds=common_bounds,
#             format="GTiff"
#         )
#
#     # Open clipped datasets
#     low_res_clipped_ds = gdal.Open(low_res_clipped_path, gdal.GA_ReadOnly)
#     high_res_clipped_ds = gdal.Open(high_res_clipped_path, gdal.GA_ReadOnly)
#
#     # Get the panchromatic (high-resolution) band
#     high_res_band = high_res_clipped_ds.GetRasterBand(1)
#     if not high_res_band:
#         raise RuntimeError("High-resolution image does not contain any bands.")
#
#     # Get all bands from the low-resolution dataset as a sequence
#     low_res_bands = [low_res_clipped_ds.GetRasterBand(i + 1) for i in range(low_res_clipped_ds.RasterCount)]
#     if not low_res_bands:
#         raise RuntimeError("Low-resolution image does not contain any bands.")
#
#     # Get nodata value from the first band of the high-resolution image
#     high_res_band_nodata = high_res_band.GetNoDataValue()
#     if high_res_band_nodata is None:
#         high_res_band_nodata = 0  # Default to 0 if no nodata value is set
#
#     # Perform pansharpening using GDAL
#     pansharpen = gdal.CreatePansharpenedVRT(
#         "",  # Empty string for in-memory VRT
#         high_res_band,
#         low_res_bands,
#     )
#
#     if pansharpen is None:
#         raise RuntimeError("Pansharpening failed.")
#
#     # Create output TIFF file
#     driver = gdal.GetDriverByName("GTiff")
#     output_ds = driver.CreateCopy(output_image_path, pansharpen, strict=0)
#
#     if output_ds is None:
#         raise RuntimeError("Failed to create the output file.")
#
#     # Set nodata value for all bands in the output
#     for i in range(output_ds.RasterCount):
#         output_ds.GetRasterBand(i + 1).SetNoDataValue(high_res_band_nodata)
#
#     # Carry over metadata from the low-resolution image
#     output_ds.SetMetadata(low_res_clipped_ds.GetMetadata())
#
#     # Close datasets
#     output_ds = None
#     low_res_clipped_ds = None
#     high_res_clipped_ds = None
#     low_res_ds = None
#     high_res_ds = None
#
#     # Clean up temporary files
#     os.remove(low_res_clipped_path)
#     os.remove(high_res_clipped_path)
#
#     print(f"Pansharpened image saved to: {output_image_path}")




import orthority as oty

# def MAIN_pansharpen_image(input_low_resolution_path, input_high_resolution_path, output_image_path):
#     pan_sharp = oty.PanSharpen(
#         input_high_resolution_path,  # panchromatic (high-res)
#         input_low_resolution_path    # multispectral (low-res)
#     )
#     pan_sharp.process(output_image_path)


from osgeo import gdal
from osgeo_utils import gdal_pansharpen

def MAIN_pansharpen_image(input_low_resolution_path, input_high_resolution_path, output_image_path):
    output_image_dir = os.path.dirname(output_image_path)
    if not os.path.exists(output_image_dir):
        os.makedirs(output_image_dir)

    # Determine how many bands the MS image has
    ms_ds = gdal.Open(input_low_resolution_path, gdal.GA_ReadOnly)
    ms_band_count = ms_ds.RasterCount
    ms_ds = None  # close after reading metadata

    # Generate band indices [1..ms_band_count]
    band_nums = list(range(1, ms_band_count + 1))

    # Use GDAL's gdal_pansharpen utility to perform pan-sharpening
    gdal_pansharpen.gdal_pansharpen(
        pan_name=input_high_resolution_path,
        spectral_names=[input_low_resolution_path],
        band_nums=band_nums,         # pansharpen all bands
        dst_filename=output_image_path,
        resampling="bilinear",
    )


import os
from envipyengine import Engine
from envipyengine.config import set as set_config

def MAIN_pansharpen_image_with_envi(input_low_resolution_path, input_high_resolution_path, output_image_path, envi_engine, output_image_path_to_delete=None):

    if output_image_path_to_delete:
        try:
            if os.path.exists(output_image_path_to_delete):
                os.remove(output_image_path_to_delete)
                print(f"Deleted existing output image: {output_image_path_to_delete}")
        except Exception as e:
            print(f"Error deleting output image {output_image_path_to_delete}: {e}")

    # Define parameters for Gram-Schmidt pansharpening
    pansharpen_params = {
        'INPUT_LOW_RESOLUTION_RASTER': {
            'url': input_low_resolution_path,
            'factory': 'URLRaster'
        },
        'INPUT_HIGH_RESOLUTION_RASTER': {
            'url': input_high_resolution_path,
            'factory': 'URLRaster'
        },
        # 'INPUT_HIGH_RESOLUTION_RASTER ':input_high_resolution_path,
        # 'INPUT_LOW_RESOLUTION_RASTER ': input_low_resolution_path,
        'OUTPUT_RASTER_URI': output_image_path,
    }

    try:
        print(f"Starting Gram-Schmidt pansharpening for:\n  MS Image: {input_low_resolution_path}\n  PAN Image: {input_high_resolution_path}")

        # Access the pansharpening task
        task = envi_engine.task('GramSchmidtPanSharpening')

        # Execute the pansharpening task with the defined parameters
        result = task.execute(pansharpen_params)

        print(f"Pansharpening complete. Output saved to: {output_image_path}")

    except Exception as e:
        print(f"Error during pansharpening: {e}")


import os
import tempfile
from osgeo import gdal

# def MAIN_pansharpen_image_with_envi(
#     input_low_resolution_path,
#     input_high_resolution_path,
#     output_image_path,
#     envi_engine
# ):
#     # Create a temporary filename but close & remove the file immediately
#     fd, temp_dat = tempfile.mkstemp(suffix=".dat")
#     os.close(fd)
#     if os.path.exists(temp_dat):
#         os.remove(temp_dat)
#
#     pansharpen_params = {
#         'INPUT_LOW_RESOLUTION_RASTER': {'url': input_low_resolution_path, 'factory': 'URLRaster'},
#         'INPUT_HIGH_RESOLUTION_RASTER': {'url': input_high_resolution_path, 'factory': 'URLRaster'},
#         'OUTPUT_RASTER_URI': temp_dat
#     }
#
#     try:
#         # Run ENVI Gram-Schmidt PanSharpening
#         task = envi_engine.task('GramSchmidtPanSharpening')
#         task.execute(pansharpen_params)
#
#         # Convert ENVI .dat output to GeoTIFF
#         gdal.Translate(output_image_path, temp_dat, format="GTiff")
#
#         # Remove ENVI output files
#         for ext in (".dat", ".hdr", ".aux.xml"):
#             f = temp_dat.replace(".dat", ext)
#             if os.path.exists(f):
#                 os.remove(f)
#
#         print(f"Pansharpened image saved to {output_image_path}")
#     except Exception as e:
#         print(f"Error during pansharpening: {e}")
