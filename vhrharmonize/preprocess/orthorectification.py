import json
import os
from tempfile import NamedTemporaryFile
from osgeo import gdal
from typing import Union, Tuple

import orthority as oty
import pyproj

from vhrharmonize.logging_utils import log


def resolve_output_resolution_for_crs(
    output_epsg: int,
    product_resolution: Union[float, None],
):
    """Return a safe output resolution for the requested target CRS."""
    if product_resolution is None:
        return None
    crs = pyproj.CRS.from_epsg(int(output_epsg))
    if crs.is_geographic:
        return None
    return product_resolution


def gcp_refined_rpc_orthorectification(
    input_image_path,
    output_image_path,
    dem_image_path,
    output_epsg,
    gcp_geojson_file_path=None,
    output_nodata_value=None,
    dtype=None,
    output_resolution: Union[float, Tuple[float, float]]=None,
    log_to_console: bool = False,
    ):
    """Orthorectify an image using RPC metadata, with optional GCP-based RPC refinement."""
    log("Running orthorectification", enabled=log_to_console, step="orthorectification")

    # Set resolution of output raster
    if isinstance(output_resolution, float):
        x_res = y_res = output_resolution
    elif isinstance(output_resolution, tuple):
        x_res, y_res = output_resolution
    else:
        x_res, y_res = None, None

    # Ensure the output directory exists
    output_image_dir = os.path.dirname(output_image_path)
    if not os.path.exists(output_image_dir):
        os.makedirs(output_image_dir)

    camera = None
    if gcp_geojson_file_path:
        # Step 1: Refine the RPC model with Orthority if GCPs are provided
        cameras = oty.RpcCameras.from_images([input_image_path])
        cameras.refine(gcp_geojson_file_path)
        camera = cameras.get(input_image_path)._rpc  # Extract the refined RPC camera

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

    # Step 3: Update the temporary file
    temp_dataset = gdal.Open(temp_image_path, gdal.GA_Update)
    input_dataset = gdal.Open(input_image_path, gdal.GA_ReadOnly)

    # Apply NoData value if specified (unconditionally)
    if output_nodata_value is not None:
        for band_idx in range(1, temp_dataset.RasterCount + 1):
            band = temp_dataset.GetRasterBand(band_idx)
            band.SetNoDataValue(output_nodata_value)

    # If RPC refinement was done, add refined RPC metadata
    if camera is not None:
        # Clear any existing RPC metadata and add the refined RPC metadata
        temp_dataset.SetMetadata(None, "RPC")
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
    else:
        # Preserve source RPC metadata for standard RPC orthorectification path.
        if input_dataset is not None:
            source_rpc_metadata = input_dataset.GetMetadata("RPC")
            if source_rpc_metadata:
                temp_dataset.SetMetadata(source_rpc_metadata, "RPC")

    # Close the dataset to write changes
    temp_dataset = None
    input_dataset = None

    # Step 4: Perform orthorectification using GDAL Warp
    warp_options = gdal.WarpOptions(
        format="GTiff",
        dstSRS=f"EPSG:{output_epsg}",
        rpc=True,
        xRes=x_res,
        yRes=y_res,
        transformerOptions=[f"RPC_DEM={dem_image_path}"],
        resampleAlg="bilinear",
        copyMetadata=True,  # Copies metadata excluding original RPC
    )
    warp_result = gdal.Warp(
        destNameOrDestDS=output_image_path,
        srcDSOrSrcDSTab=temp_image_path,
        options=warp_options,
    )
    if warp_result is None:
        raise RuntimeError(
            "GDAL Warp failed during orthorectification. "
            "Check that the atmospheric-corrected input still contains valid RPC metadata."
        )
    warp_result = None

    # Step 3: Update the temporary file
    dataset = gdal.Open(output_image_path, gdal.GA_Update)
    if dataset is None:
        raise RuntimeError(f"Orthorectification output was not created: {output_image_path}")

    # Apply NoData value if specified (unconditionally)
    if output_nodata_value is not None:
        for band_idx in range(1, dataset.RasterCount + 1):
            band = dataset.GetRasterBand(band_idx)
            band.SetNoDataValue(output_nodata_value)

    # Clean up the temporary file
    os.remove(temp_image_path)
    log("Wrote output", enabled=log_to_console, step="orthorectification")


def qgis_gcps_to_csv(
    input_gcp_path,
    output_epsg=None,
    log_to_console: bool = False,
    ):

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
        transformer = pyproj.Transformer.from_crs(
            "EPSG:4326",
            f"EPSG:{output_epsg}",
            always_xy=True
        )
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
                mapX = float(parts[0])
                mapY = float(parts[1])
                sourceX = float(parts[2])
                sourceY = float(parts[3])

                # Transform map coordinates if transformer is available
                if transformer:
                    transformed_mapX, transformed_mapY = transformer.transform(mapX, mapY)
                    mapX_str = f"{transformed_mapX}"
                    mapY_str = f"{transformed_mapY}"
                else:
                    mapX_str = f"{mapX}"
                    mapY_str = f"{mapY}"

                # Format source coordinates without modification (pixel coordinates)
                sourceX_str = f"{sourceX}"
                sourceY_str = f"{sourceY}"

                # Append to CSV text
                output_csv_text += f"{mapX_str},{mapY_str},{sourceX_str},{sourceY_str},0\n"
            else:
                continue

    log("Converted QGIS GCP text to CSV", enabled=log_to_console, step="orthorectification")
    return output_csv_text


def geo_to_image_coords(
    dataset,
    x,
    y
    ):
    """Transform map coordinates into image pixel coordinates using dataset GCPs."""

    transform_options = ["METHOD=GCP_POLYNOMIAL"] # "METHOD=GCP_TPS" or "METHOD=GCP_POLYNOMIAL"

    transformer = gdal.Transformer(dataset, None, transform_options)
    if not transformer:
        raise RuntimeError("Failed to create GDAL Transformer with GCPs.")

    # TransformPoint( bDstToSrc, X_in, Y_in ).
    success, (px, py, pz) = transformer.TransformPoint(1, float(x), float(y))
    return success, (px, py, pz)


def qgis_gcps_to_geojson(
    input_image_path,
    qgis_gcp_file_path,
    file_name,
    dem_file_path,
    output_geojson_path,
    force_positive_pixel_values=False,
    log_to_console: bool = False,
    ):
    """Convert a QGIS GCP text file into Orthority-compatible GeoJSON control points."""

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
        if (
            x_pixel < 0
            or y_pixel < 0
            or x_pixel >= dem_dataset.RasterXSize
            or y_pixel >= dem_dataset.RasterYSize
        ):
            return 0.0
        elevation = dem_band.ReadAsArray(x_pixel, y_pixel, 1, 1)[0, 0]
        no_data = dem_band.GetNoDataValue()
        return float(elevation) if elevation != no_data else 0.0

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
        if len(values) < 4:
            continue
        map_x = float(values[0])
        map_y = float(values[1])
        pixel_x = float(values[2])
        pixel_y = float(values[3])

        elevation = get_elevation(map_x, map_y)

        success, (px, py, pz) = geo_to_image_coords(image_dataset, pixel_x, pixel_y)
        if not success:
            continue

        if force_positive_pixel_values:
            px, py = abs(px), abs(py)

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
    log("Wrote GCP GeoJSON", enabled=log_to_console, step="orthorectification")

    # Cleanup
    dem_dataset = None
    image_dataset = None
