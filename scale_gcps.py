import json
from osgeo import gdal

def scale_gcps_geojson(current_gcp_image, desired_scale_image, input_geojson_path, output_geojson_path, replace_filename=None):
    """
    Scale the 'ji' attribute in a GeoJSON based on the resolution difference between two images.

    Parameters:
    current_gcp_image (str): Path to the image with the current GCPs.
    desired_scale_image (str): Path to the image to which GCPs should be scaled.
    input_geojson_path (str): Path to the input GeoJSON file.
    output_geojson_path (str): Path to save the output GeoJSON file with scaled GCPs.
    replace_filename (str, optional): New filename to replace the 'filename' attribute in the GeoJSON.

    Returns:
    None
    """

    # Open the current and desired images
    current_dataset = gdal.Open(current_gcp_image)
    desired_dataset = gdal.Open(desired_scale_image)

    if not current_dataset or not desired_dataset:
        raise FileNotFoundError("One or both of the image paths are invalid.")

    # Get the dimensions of the images
    current_width = current_dataset.RasterXSize
    current_height = current_dataset.RasterYSize
    desired_width = desired_dataset.RasterXSize
    desired_height = desired_dataset.RasterYSize

    # Calculate scaling factors
    x_scale = desired_width / current_width
    y_scale = desired_height / current_height

    # Load the input GeoJSON
    with open(input_geojson_path, 'r') as f:
        geojson_data = json.load(f)

    # Iterate through the features and scale the 'ji' attribute
    for feature in geojson_data['features']:
        # Scale the 'ji' attribute
        if 'ji' in feature['properties']:
            ji = feature['properties']['ji']
            scaled_ji = [ji[0] * x_scale, ji[1] * y_scale]
            feature['properties']['ji'] = scaled_ji

        # Replace the 'filename' attribute if replace_filename is provided
        if replace_filename is not None:
            feature['properties']['filename'] = replace_filename

    # Save the updated GeoJSON to the output path
    with open(output_geojson_path, 'w') as f:
        json.dump(geojson_data, f, indent=4)

    print(f"Scaled GeoJSON saved to {output_geojson_path}")
