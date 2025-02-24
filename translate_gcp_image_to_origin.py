from osgeo import gdal
import os

def translate_gcp_image_to_origin(input_image_path, output_image_path):
    """
Replace all existing GCPs in the input dataset with four GCPs that
place the top-left corner at (0,0), top-right at (width,0), bottom-left
at (0,height), and bottom-right at (width,height) in a local coordinate system.

The output is saved as a new file at 'output_image_path'.
"""

    output_dir = os.path.dirname(output_image_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 1. Open the input dataset
    src_ds = gdal.Open(input_image_path, gdal.GA_ReadOnly)
    if not src_ds:
        raise IOError(f"Cannot open file: {input_image_path}")

    # 2. Get the dimensions of the raster
    width = src_ds.RasterXSize
    height = src_ds.RasterYSize

    # 3. Define the four corner GCPs in pixel/line (col/row) vs. map coordinates
    #    GCP args: (map_x, map_y, elevation, pixel_col, pixel_row)
    top_left = gdal.GCP(0, 0, 0, 0, 0)
    top_right = gdal.GCP(width, 0, 0, width, 0)
    bottom_left = gdal.GCP(0, -height, 0, 0, height)
    bottom_right = gdal.GCP(width, -height, 0, width, height)

    new_gcps = [top_left, top_right, bottom_left, bottom_right]


    # 5. Use gdal.Translate to create a new GeoTIFF with these new GCPs
    translate_options = gdal.TranslateOptions(
        GCPs=new_gcps,
        format='GTiff'
    )

    gdal.Translate(destName=output_image_path, srcDS=src_ds, options=translate_options)

    # 6. Close the dataset
    src_ds = None

    print(f"Translated to bounds: ((0,0), ({width},0), (0,{height}), ({width,height})) for image{output_image_path}")
