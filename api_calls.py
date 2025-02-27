import requests

def download_modis_water_vapor(date, zoom_level, tile_row, tile_col, tile_matrix_set="250m", output_format="png", output_file="modis_water_vapor.png"):
    """
    Downloads a MOD05_L2 (Water Vapor) tile from NASA GIBS using the WMTS REST API.

    Parameters:
    - date (str): Date in "YYYY-MM-DD" format.
    - zoom_level (int): Zoom level (TileMatrix) for the requested tile.
    - tile_row (int): Tile row index.
    - tile_col (int): Tile column index.
    - tile_matrix_set (str): Tile matrix set (default is "250m").
    - output_format (str): Image format, e.g., "png", "jpg".
    - output_file (str): Output filename to save the image.
    
    Returns:
    - Saves the image file to the specified output location.
    """

    base_url = "https://gibs.earthdata.nasa.gov/wmts/epsg4326/best"
    layer = "MOD05_L2"

    url = f"{base_url}/{layer}/default/{date}/{tile_matrix_set}/{zoom_level}/{tile_row}/{tile_col}.{output_format}"

    response = requests.get(url, stream=True)

    if response.status_code == 200:
        with open(output_file, "wb") as file:
            for chunk in response.iter_content(chunk_size=1024):
                file.write(chunk)
        print(f"Downloaded: {output_file}")
    else:
        print(f"Failed to download tile. HTTP Status Code: {response.status_code}")

# Example usage
download_modis_water_vapor(date="2024-02-25", zoom_level=6, tile_row=13, tile_col=36)