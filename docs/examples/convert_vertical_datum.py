# Paths to the input DEM and output converted DEM
dem_path = "/mnt/x/Imagery/Lidar/Big_Island/2018_PuuWaawaa/DEM/2018_2020_bigIsland_DEM_J970216_000_000.tif"
output_path = "/mnt/d/outputconvertdem.tif"

# Define input SRS: WGS84 with EGM96 geoid model
input_srs = osr.SpatialReference()
input_srs.ImportFromEPSG(4326)  # WGS84 horizontal
input_srs.SetAttrValue("VERT_CS", "EGM96")  # Add EGM96 as a vertical reference

# Define output SRS: WGS84 3D (ellipsoidal height)
output_srs = osr.SpatialReference()
output_srs.ImportFromEPSG(4979)  # WGS84 3D (latitude, longitude, and ellipsoidal height)

# Call the conversion function
convert_vertical_datum(
    input_image_path=dem_path,
    output_image_path=output_path,
    input_override_srs=input_srs.ExportToWkt(),
    output_override_srs=output_srs.ExportToWkt()
)