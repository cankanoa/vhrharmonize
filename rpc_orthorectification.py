import os
import io
import shutil
import subprocess
import tempfile
from osgeo import gdal, osr, ogr
import orthority as oty
from skimage.data import camera


# def MAIN_rpc_orthorectification(
#         input_image_path,
#         output_image_path,
#         dem_file,
#         output_epsg,
#         cutline_gpkg,
#         gcp_csv_text=None,
#         resample_alg='cubic'
# ):
#     # Open the input image
#     input_dataset = gdal.Open(input_image_path)
#     if not input_dataset:
#         raise FileNotFoundError(f"Input image not found: {input_image_path}")
#
#
#     # Verify that the GeoPackage (cutline) exists
#     if not os.path.exists(cutline_gpkg):
#         raise FileNotFoundError(f"Cutline GeoPackage not found: {cutline_gpkg}")
#
#     os.makedirs(os.path.dirname(output_image_path), exist_ok=True)
#     data_source = ogr.Open(cutline_gpkg)
#     first_layer = data_source.GetLayer(0)
#     first_layer_name = first_layer.GetName()
#
#
#     # Clip the DEM using the cutline and save to a temporary file
#     with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as temp_dem_file:
#         temp_dem_path = temp_dem_file.name
#         gdal.Warp(
#             temp_dem_path,                     # Output raster file
#             dem_file,                      # Input raster dataset
#             format="GTiff",                    # Output format
#             cutlineDSName=cutline_gpkg,                # Input vector file (GeoPackage)
#             cutlineLayer=first_layer_name,           # Layer name in the GeoPackage
#             cropToCutline=True,                # Crop the raster to the vector mask
#             dstNodata=-9999                    # Set the nodata value for the output raster
#         )
#
#     # Add GCPs if provided
#     if gcp_csv_text:
#         gcps = []
#         csv_reader = io.StringIO(gcp_csv_text)
#         header = next(csv_reader)  # Skip the header row
#         for line in csv_reader:
#             mapX, mapY, sourceX, sourceY, z = map(float, line.strip().split(','))
#             gcps.append(gdal.GCP(mapX, mapY, z, sourceX, sourceY))
#         input_dataset.SetGCPs(gcps, input_dataset.GetProjection())
#
#     # Set up the options for orthorectification
#     warp_options = gdal.WarpOptions(
#         format='GTiff',
#         # rpc=True,  # Enable RPC processing
#         dstSRS=f"EPSG:{output_epsg}", # 4326 is WGS84
#         resampleAlg=resample_alg,
#         transformerOptions=[f"RPC_DEM={temp_dem_path}"],  # Use DEM for orthorectification
#     )
#
#     # Perform the warp (orthorectification)
#     gdal.Warp(output_image_path, input_dataset, options=warp_options)
#
#     # Clean up temporary files
#     os.remove(temp_dem_path)
#     cutline_ds = None
#     input_dataset = None
#
#     print(f"Orthorectified image saved to: {output_image_path}")



def MAIN_rpc_orthorectification(input_image_path, output_image_path, dem_image_path, output_epsg):
    # Step 4: Perform orthorectification using GDAL
    warp_options = gdal.WarpOptions(
        format="GTiff",
        dstSRS=f"EPSG:{output_epsg}",
        rpc=True,
        transformerOptions=[f"RPC_DEM={dem_image_path}"],
        resampleAlg="cubic",
        copyMetadata=True,  # Copy metadata excluding original RPC
    )

    # Perform the warp
    gdal.Warp(
        destNameOrDestDS=output_image_path,
        srcDSOrSrcDSTab=input_image_path,
        options=warp_options,
    )



import os
from tempfile import NamedTemporaryFile
import os
import json
import re
from osgeo import gdal
from tempfile import NamedTemporaryFile

import os
from osgeo import gdal
from tempfile import NamedTemporaryFile

def MAIN_gcp_refined_rpc_orthorectification(
    input_image_path,
    output_image_path,
    gcp_geojson_file_path,
    dem_image_path,
    output_epsg,
    output_nodata_value=None,
    dtype=None,
    external_rpc_image_path=None
):
    import orthority as oty

    # Ensure the output directory exists
    output_image_dir = os.path.dirname(output_image_path)
    if not os.path.exists(output_image_dir):
        os.makedirs(output_image_dir)

    # Step 1: Refine the RPC model with Orthority
    if external_rpc_image_path:
        cameras = oty.RpcCameras.from_images([external_rpc_image_path])
    else:
        cameras = oty.RpcCameras.from_images([input_image_path])


    cameras.refine(gcp_geojson_file_path)
    camera = cameras.get(external_rpc_image_path)
    camera = camera._rpc  # Extract the refined RPC camera

    # Step 2: Create a temporary file for the modified input image
    with NamedTemporaryFile(delete=False, suffix=".tif") as temp_file:
        temp_image_path = temp_file.name

    # Build TranslateOptions to remove all metadata, optionally applying data type
    translate_kwargs = {
        "format": "GTiff",
        "metadataOptions": ["REMOVE_ALL"],  # Remove all metadata
    }

    if dtype is not None:
        gdal_dtype = gdal.GetDataTypeByName(dtype)
        if gdal_dtype is None:
            raise ValueError(f"Unsupported GDAL data type: {dtype}")
        translate_kwargs["outputType"] = gdal_dtype

    translate_options = gdal.TranslateOptions(**translate_kwargs)
    gdal.Translate(temp_image_path, input_image_path, options=translate_options)

    # Step 3: Add refined RPC metadata (and optional NoData) to the temporary file
    temp_dataset = gdal.Open(temp_image_path, gdal.GA_Update)

    # If we have a NoData value, apply it to each band
    if output_nodata_value is not None:
        for band_idx in range(1, temp_dataset.RasterCount + 1):
            band = temp_dataset.GetRasterBand(band_idx)
            band.SetNoDataValue(output_nodata_value)

    # Clear any existing RPC metadata
    temp_dataset.SetMetadata(None, "RPC")

    # Add the refined RPC metadata
    refined_rpc_metadata = {
        "LINE_OFF": str(camera.line_off),
        "LINE_SCALE": str(camera.line_scale),
        "SAMP_OFF": str(camera.samp_off),
        "SAMP_SCALE": str(camera.samp_scale),
        "LAT_OFF": str(camera.lat_off),
        "LAT_SCALE": str(camera.lat_scale),
        "LONG_OFF": str(camera.long_off),
        "LONG_SCALE": str(camera.long_scale),
        "HEIGHT_OFF": str(camera.height_off),
        "HEIGHT_SCALE": str(camera.height_scale),
        "LINE_NUM_COEFF": " ".join(map(str, camera.line_num_coeff)),
        "LINE_DEN_COEFF": " ".join(map(str, camera.line_den_coeff)),
        "SAMP_NUM_COEFF": " ".join(map(str, camera.samp_num_coeff)),
        "SAMP_DEN_COEFF": " ".join(map(str, camera.samp_den_coeff)),
    }
    temp_dataset.SetMetadata(refined_rpc_metadata, "RPC")

    # Close the temp dataset to ensure changes are written
    temp_dataset = None

    # Step 4: Perform orthorectification using GDAL
    warp_options = gdal.WarpOptions(
        format="GTiff",
        dstSRS=f"EPSG:{output_epsg}",
        rpc=True,
        transformerOptions=[f"RPC_DEM={dem_image_path}"],
        resampleAlg="bilinear",
        copyMetadata=True,  # Copies metadata excluding original RPC
    )

    # Perform the warp
    gdal.Warp(
        destNameOrDestDS=output_image_path,
        srcDSOrSrcDSTab=temp_image_path,
        options=warp_options,
    )

    # Clean up the temporary file
    os.remove(temp_image_path)

    print(f"Orthorectified image saved to {output_image_path}")












import pyproj

def qgis_gcps_to_csv(input_gcp_path, output_epsg=None):
    """
    Convert QGIS GCP format to GDAL-compatible CSV format, with optional map coordinate transformation.

    Parameters:
        input_gcp_path (str): Path to the QGIS GCP file.
        output_epsg (int or str, optional): EPSG code to transform map coordinates to.
                                            If None, map coordinates are not transformed.

    Returns:
        str: CSV text containing GDAL-compatible GCPs.
    """
    output_csv_text = "mapX,mapY,sourceX,sourceY,z\n"

    # Initialize transformer if output_epsg is provided
    if output_epsg is not None:
        try:
            # Assumption: Input map coordinates are in WGS84 (EPSG:4326)
            transformer = pyproj.Transformer.from_crs(
                "EPSG:4326",          # Input CRS
                f"EPSG:{output_epsg}",# Output CRS
                always_xy=True
            )
            print(f"Initialized transformer from EPSG:4326 to EPSG:{output_epsg}.")
        except Exception as e:
            print(f"Error initializing transformer: {e}")
            return output_csv_text
    else:
        transformer = None

    with open(input_gcp_path, 'r') as infile:
        for line_number, line in enumerate(infile, start=1):
            stripped_line = line.strip()

            # Skip metadata lines or empty lines
            if stripped_line.startswith('#') or stripped_line == '':
                continue

            first_char = stripped_line[0]
            if not (first_char.isdigit() or first_char == '-'):
                continue

            parts = stripped_line.split(',')
            if len(parts) >= 8:
                try:
                    mapX = float(parts[0])
                    mapY = float(parts[1])
                    sourceX = float(parts[2])
                    sourceY = float(parts[3])
                    # The remaining parts are enable, dX, dY, residual
                except ValueError:
                    print(f"Line {line_number}: Invalid numeric values. Skipping line.")
                    continue

                # Transform map coordinates if transformer is available
                if transformer:
                    try:
                        transformed_mapX, transformed_mapY = transformer.transform(mapX, mapY)
                        mapX_str = f"{transformed_mapX}"
                        mapY_str = f"{transformed_mapY}"
                    except Exception as e:
                        print(f"Line {line_number}: Error transforming coordinates ({mapX}, {mapY}): {e}")
                        # Optionally, skip this GCP or retain original coordinates
                        mapX_str = f"{mapX}"
                        mapY_str = f"{mapY}"
                else:
                    mapX_str = f"{mapX}"
                    mapY_str = f"{mapY}"

                # Format source coordinates without modification (pixel coordinates)
                sourceX_str = f"{sourceX}"
                sourceY_str = f"{sourceY}"

                # Append to CSV text
                output_csv_text += f"{mapX_str},{mapY_str},{sourceX_str},{sourceY_str},0\n"
            else:
                print(f"Line {line_number}: Insufficient columns. Expected at least 8, got {len(parts)}. Skipping line.")
                continue

    print("Converted QGIS GCP file to GDAL-compatible CSV text.")
    return output_csv_text


# import json
# from osgeo import gdal
#
# def qgis_gcps_to_geojson(input_image_path, qgis_gcp_file_path, file_name, dem_file_path, output_geojson_path, force_positive_pixel_values=False):
#     geojson = {
#         "type": "FeatureCollection",
#         "features": []
#     }
#
#     # Open the DEM file
#     dem_dataset = gdal.Open(dem_file_path)
#     if not dem_dataset:
#         raise FileNotFoundError(f"Could not open DEM file: {dem_file_path}")
#     dem_band = dem_dataset.GetRasterBand(1)
#     dem_transform = dem_dataset.GetGeoTransform()
#
#     def get_elevation(lon, lat):
#         # Transform geographic coordinates to pixel coordinates
#         x_pixel = int((lon - dem_transform[0]) / dem_transform[1])
#         y_pixel = int((lat - dem_transform[3]) / dem_transform[5])
#
#         try:
#             elevation = dem_band.ReadAsArray(x_pixel, y_pixel, 1, 1)[0, 0]
#             return float(elevation) if elevation != dem_band.GetNoDataValue() else 0.0
#         except:
#             return 0.0
#
#     with open(qgis_gcp_file_path, 'r') as file:
#         lines = file.readlines()
#
#     for line in lines:
#         line = line.strip()
#         if not line or line.startswith("#") or not (line[0].isdigit() or line[0] == '-'):
#             continue
#
#         values = line.split(",")
#         try:
#             map_x = float(values[0])  # Longitude
#             map_y = float(values[1])  # Latitude
#             source_x = float(values[2])  # Pixel X
#             source_y = float(values[3])  # Pixel Y
#
#             # Apply absolute value if force_positive_pixel_values is True
#             if force_positive_pixel_values:
#                 source_x = abs(source_x)
#                 source_y = abs(source_y)
#         except (IndexError, ValueError):
#             continue
#
#         # Get elevation from DEM file
#         elevation = get_elevation(map_x, map_y)
#
#         feature = {
#             "type": "Feature",
#             "properties": {
#                 "filename": file_name,
#                 "ji": [source_x, source_y],
#                 "id": f"gcp-point-{len(geojson['features']) + 1}",
#                 "info": f"gcp-point-{len(geojson['features']) + 1}"
#             },
#             "geometry": {
#                 "type": "Point",
#                 "coordinates": [map_x, map_y, elevation]
#             }
#         }
#
#         geojson["features"].append(feature)
#
#     with open(output_geojson_path, "w") as f:
#         json.dump(geojson, f, indent=4)
#     print(f"GeoJSON saved to {output_geojson_path}")
#
#     # Close the DEM dataset
#     dem_dataset = None



import json
from osgeo import gdal


import numpy as np
from rpcm import rpc_from_geotiff

# def geo_to_image_coords(lon, lat, alt, image_path):
#     """
#     Converts geographic coordinates (longitude, latitude, elevation) to pixel coordinates (column, row) using RPC metadata via rpcm.
#
#     Args:
#         lon (float or array-like): Longitude(s) of the geographic coordinate(s).
#         lat (float or array-like): Latitude(s) of the geographic coordinate(s).
#         alt (float or array-like): Elevation(s) of the geographic coordinate(s).
#         image_path (str): Path to the input image with RPC metadata.
#
#     Returns:
#         tuple: (pixel_x, pixel_y) corresponding to the image coordinates (column, row).
#     """
#     # Load RPC metadata from the GeoTIFF file
#     rpc = rpc_from_geotiff(image_path)
#
#     # Convert to numpy arrays for vectorized operations
#     lon = np.array(lon, dtype=float)
#     lat = np.array(lat, dtype=float)
#     alt = np.array(alt, dtype=float)
#
#     # Perform the projection to get pixel coordinates
#     x, y = rpc.projection(lon, lat, alt)
#
#     return x, y
#
#
# def qgis_gcps_to_geojson(input_image_path, qgis_gcp_file_path, file_name, dem_file_path, output_geojson_path, force_positive_pixel_values=False):
#     geojson = {
#         "type": "FeatureCollection",
#         "features": []
#     }
#
#     # Open the DEM file
#     dem_dataset = gdal.Open(dem_file_path)
#     if not dem_dataset:
#         raise FileNotFoundError(f"Could not open DEM file: {dem_file_path}")
#     dem_band = dem_dataset.GetRasterBand(1)
#     dem_transform = dem_dataset.GetGeoTransform()
#
#     def get_elevation(lon, lat):
#         x_pixel = int((lon - dem_transform[0]) / dem_transform[1])
#         y_pixel = int((lat - dem_transform[3]) / dem_transform[5])
#
#         try:
#             elevation = dem_band.ReadAsArray(x_pixel, y_pixel, 1, 1)[0, 0]
#             return float(elevation) if elevation != dem_band.GetNoDataValue() else 0.0
#         except:
#             return 0.0
#
#     # Open the input image and fetch RPC metadata
#     image_dataset = gdal.Open(input_image_path)
#     if not image_dataset:
#         raise FileNotFoundError(f"Could not open input image file: {input_image_path}")
#     rpc_model = image_dataset.GetMetadata("RPC")
#     if not rpc_model:
#         raise ValueError(f"Input image does not contain RPC metadata: {input_image_path}")
#
#     with open(qgis_gcp_file_path, 'r') as file:
#         lines = file.readlines()
#
#     for line in lines:
#         line = line.strip()
#         if not line or line.startswith("#") or not (line[0].isdigit() or line[0] == '-'):
#             continue
#
#         values = line.split(",")
#         try:
#             map_x = float(values[0])  # Longitude
#             map_y = float(values[1])  # Latitude
#
#             # Get elevation from DEM file
#             elevation = get_elevation(map_x, map_y)
#
#             # Convert geographic coordinates to pixel coordinates using RPCM
#             source_x, source_y = geo_to_image_coords(map_x, map_y, elevation, input_image_path)
#
#             # Apply absolute value if force_positive_pixel_values is True
#             if force_positive_pixel_values:
#                 source_x = abs(source_x)
#                 source_y = abs(source_y)
#
#         except (IndexError, ValueError) as e:
#             print(f"Skipping line due to error: {e}")
#             continue
#
#         feature = {
#             "type": "Feature",
#             "properties": {
#                 "filename": file_name,
#                 "ji": [source_x, source_y],
#                 "id": f"gcp-point-{len(geojson['features']) + 1}",
#                 "info": f"gcp-point-{len(geojson['features']) + 1}"
#             },
#             "geometry": {
#                 "type": "Point",
#                 "coordinates": [map_x, map_y, elevation]
#             }
#         }
#
#         geojson["features"].append(feature)
#
#     with open(output_geojson_path, "w") as f:
#         json.dump(geojson, f, indent=4)
#     print(f"GeoJSON saved to {output_geojson_path}")
#
#     # Close the datasets
#     dem_dataset = None
#     image_dataset = None

import numpy as np
from osgeo import gdal

def geo_to_image_coords(dataset, x, y):
    transform_options = ["METHOD=GCP_POLYNOMIAL"] # "METHOD=GCP_TPS" or "METHOD=GCP_POLYNOMIAL"

    transformer = gdal.Transformer(dataset, None, transform_options)
    if not transformer:
        raise RuntimeError("Failed to create GDAL Transformer with GCPs.")

    # TransformPoint( bDstToSrc, X_in, Y_in ).
    success, (px, py, pz) = transformer.TransformPoint(1, float(x), float(y))
    return success, (px, py, pz)

def qgis_gcps_to_geojson(input_image_path, qgis_gcp_file_path, file_name, dem_file_path, output_geojson_path, force_positive_pixel_values=False):
    import json

    geojson = {
        "type": "FeatureCollection",
        "features": []
    }

    # Open the DEM file
    dem_dataset = gdal.Open(dem_file_path)
    if not dem_dataset:
        raise FileNotFoundError(f"Could not open DEM file: {dem_file_path}")
    dem_band = dem_dataset.GetRasterBand(1)
    dem_transform = dem_dataset.GetGeoTransform()

    def get_elevation(lon, lat):
        # Convert lon/lat into pixel coords within the DEM
        x_pixel = int((lon - dem_transform[0]) / dem_transform[1])
        y_pixel = int((lat - dem_transform[3]) / dem_transform[5])
        try:
            elevation = dem_band.ReadAsArray(x_pixel, y_pixel, 1, 1)[0, 0]
            no_data = dem_band.GetNoDataValue()
            return float(elevation) if elevation != no_data else 0.0
        except:
            return 0.0

    # Open the input image
    image_dataset = gdal.Open(input_image_path)
    if not image_dataset:
        raise FileNotFoundError(f"Could not open input image file: {input_image_path}")

    # Parse QGIS GCP file
    with open(qgis_gcp_file_path, 'r') as file:
        lines = file.readlines()

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Skip if the line doesn't start with digit or '-'
        if not (line[0].isdigit() or line[0] == '-'):
            continue

        values = line.split(",")
        try:
            # map_x, map_y are geographic coords (long, lat)
            map_x = float(values[0])
            map_y = float(values[1])
            pixel_x = float(values[2])
            pixel_y = float(values[3])

            elevation = get_elevation(map_x, map_y)

            # Convert from geographic (lon, lat) to pixel coords
            success, (px, py, pz) = geo_to_image_coords(image_dataset, pixel_x, pixel_y)
            if not success:
                print(f"Transformation failed for {pixel_x}, {pixel_y}")
                continue

            if force_positive_pixel_values:
                px, py = abs(px), abs(py)

        except (IndexError, ValueError) as e:
            print(f"Skipping line due to error: {e}")
            continue

        feature = {
            "type": "Feature",
            "properties": {
                "filename": file_name,
                # 'ji' can hold your pixel coords
                "ji": [px, py],
                "id": f"gcp-point-{len(geojson['features']) + 1}",
                "info": f"gcp-point-{len(geojson['features']) + 1}"
            },
            "geometry": {
                "type": "Point",
                "coordinates": [map_x, map_y, elevation]
            }
        }

        geojson["features"].append(feature)

    with open(output_geojson_path, "w") as f:
        json.dump(geojson, f, indent=4)
    print(f"GeoJSON saved to {output_geojson_path}")

    # Cleanup
    dem_dataset = None
    image_dataset = None
