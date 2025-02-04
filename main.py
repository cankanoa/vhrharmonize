# ------------------- Imports
import os
import json
from tqdm import tqdm
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from osgeo import gdal

from SatelliteProcess.replace_band_continuous_value import replace_band_continuous_values_in_largest_segment
from rpc_orthorectification import qgis_gcps_to_csv, qgis_gcps_to_geojson, MAIN_gcp_refined_rpc_orthorectification
from test_flaash_params import create_test_flaash_params
from translate_gcp_image_to_origin import MAIN_translate_gcp_image_to_origin
from rpc_orthorectification import MAIN_rpc_orthorectification
from pansharpen import MAIN_pansharpen_image
from flaash import create_flaash_params, MAIN_run_flaash
from pyproj import Transformer
from envipyengine import Engine
from flaash import convert_dn_to_radiance_with_envi
from global_match import process_global_histogram_matching
from local_match import process_local_histogram_matching
from scale_gcps import scale_gcps_geojson
from pansharpen import MAIN_pansharpen_image_with_envi
from get_image_largest_value import get_image_largest_value

import envipyengine.config
from plot_pixel_distribution import plot_pixel_distribution
envipyengine.config.set('engine', "/mnt/c/Program Files/Harris/ENVI57/IDL89/bin/bin.x86_64/taskengine.exe")
envi_engine = Engine('ENVI')
envi_engine.tasks()






# Troubleshooting
# Original files are found by searching through all files in input_folders_array for any ".IMD" endings, as such, they must be present and other IMD files should not be added





# ------------------- Define variables
# Finds IMD files
input_folders_array = [
    "/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445319010",
    "/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445318010",
    # '/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20170705_35cm_WV03_BAB_016445286010'
]

# Only process these iamges if set, if not set, all images will be processed
filter_basenames = [
    '17DEC08211758-M1BS-016445319010_01_P003',
    '17DEC08211800-M1BS-016445319010_01_P004',
    '17DEC08211801-M1BS-016445319010_01_P005',
    '17DEC08211840-M1BS-016445318010_01_P015',
    '17DEC08211841-M1BS-016445318010_01_P016',
    ]

epsg = 6635
nodata_value = -9999
dtype = 'int16'

# Required to create a Root_WV03.txt file at each satellite image root folder; optional add override to params per scene and photo in the following format:
# {     "ParamsOverridesPerScene": {
#         "WATER_VAPOR_PRESET": 0.89
#     },
#     "ParamsOverridesPerPhoto": {
#         "17DEC08211758-M1BS-016445319010_01_P003": {
#             "Key": "value"}}}






# ------------------- Start processing
def main_process_imagery(input_folders_array):
    for input_folder in (input_folders_array):
        root_paths = find_roots(input_folder)
        for root_folder_path, root_file_path in root_paths:
            found_default_files = find_files(root_folder_path, root_file_path)
            for count, (photo_basename, found_default_file) in enumerate(found_default_files.items(), start=0):
                print(f"Processing {found_default_file['mul_photo_basename']}")

                # -------------------- Define variables per image
                print('Fetching variables from files')
                root_folder_path = found_default_file['root_folder_path']
                root_file_path = found_default_file['root_file_path']
                # photo_basename = found_default_file['photo_basename']

                mul_photo_basename = found_default_file['mul_photo_basename']
                mul_scene_basename = found_default_file['mul_scene_basename']
                mul_til_file = found_default_file['mul_til_file']
                mul_tif_file = found_default_file['mul_tif_file']
                mul_imd_file = found_default_file['mul_imd_file']
                mul_photo_info_file = found_default_file['mul_photo_info_file']

                pan_scene_basename = found_default_file['pan_scene_basename']
                pan_til_file = found_default_file['pan_til_file']
                pan_tif_file = found_default_file['pan_tif_file']
                pan_imd_file = found_default_file['pan_imd_file']
                pan_photo_basename = found_default_file['pan_photo_basename']
                pan_photo_info_file = found_default_file['pan_photo_info_file']

                params_overrides_scene, mul_params_overrides_photo, mul_imd_data = get_metadata_from_files(root_file_path, mul_imd_file, mul_photo_basename)
                params_overrides_scene, pan_params_overrides_photo, pan_imd_data = get_metadata_from_files(root_file_path, pan_imd_file, pan_photo_basename)

                mul_shp_path = found_default_file['mul_shp_file']
                pan_shp_path = found_default_file['pan_shp_file']


                # -------------------- Set any custom resources per image (array must match legth of found_default_files)
                print('Setting custom variables')
                if count == 0:
                    ortho_from_oab_reference_imagery = [
                        # '/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445319010/FLAASH/17DEC08211758-M1BS-016445319010_01_P003_FLAASH.dat',
                        # '/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445319010/FLAASH/17DEC08211800-M1BS-016445319010_01_P004_FLAASH.dat',
                        # '/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445319010/FLAASH/17DEC08211801-M1BS-016445319010_01_P005_FLAASH.dat',
                        # '/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445318010/20171208_36cm_WV03_BAB_016445318010_FLAASH/17DEC08211840-M1BS-016445318010_01_P015_FLAASH.dat',
                        # '/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445318010/20171208_36cm_WV03_BAB_016445318010_FLAASH/17DEC08211841-M1BS-016445318010_01_P016_FLAASH.dat',
                    ]
                    reference_image_path = [
                        '/mnt/x/PROJECTS_2/Big_Island/LandCover/Input/satellite/20171019_36cm_WV03_OAB_Georef/20171019_36cm_WV03_OAB_Pansharp_Quac_Georef.tif',
                        '/mnt/x/PROJECTS_2/Big_Island/LandCover/Input/satellite/20171019_36cm_WV03_OAB_Georef/20171019_36cm_WV03_OAB_Pansharp_Quac_Georef.tif',
                        '/mnt/x/PROJECTS_2/Big_Island/LandCover/Input/satellite/20171019_36cm_WV03_OAB_Georef/20171019_36cm_WV03_OAB_Pansharp_Quac_Georef.tif',
                        '/mnt/x/PROJECTS_2/Big_Island/LandCover/Input/satellite/20171208_36cm_WV03_OAB_Georef/20171208_36cm_WV03_OAB_Pansharp_Quac_Georef.tif',
                        '/mnt/x/PROJECTS_2/Big_Island/LandCover/Input/satellite/20171208_36cm_WV03_OAB_Georef/20171208_36cm_WV03_OAB_Pansharp_Quac_Georef.tif'
                    ]

                    # DEM file needs to be in WGS84 elipsoidal height
                    # dem_file_path = '/mnt/x/Imagery/Lidar/Big_Island/2018_PuuWaawaa/DEM/2018_2020_bigIsland_DEM_J970216_000_000.tif'#'/mnt/x/Imagery/Lidar/Big_Island/2018_PuuWaawaa/DEM/2018_2020_bigIsland_DEM_J970216_000_000.tif'
                    dem_file_path = '/mnt/d/demlast.tif'
                    # dem_file_path = '/mnt/d/dem_WGS84_Elipsoid.tif'
                    # dem_file_path = '/mnt/d/srtm_05_09/srtm_05_09.tif'
                    # dem_file_path = '/mnt/x/Imagery/Elevation/rasters_SRTMGL1Ellip/output_SRTMGL1Ellip.tif'



                # -------------------- Start processing per image
                print('----------Starting init')# -------------------- Initialize files
                # These shp files from maxar dont have a projection, this adds the correct projectionfor shp_path in [found_default_file['mul_shp_file'], found_default_file['pan_shp_file']]:
                mul_gpkg_path = os.path.join(input_folder, 'Mul_Footprint', f"{mul_photo_basename}.gpkg")
                shp_gpkg_path = os.path.join(input_folder, 'Pan_Footprint', f"{pan_photo_basename}.gpkg")
                mul_original_path = os.path.join(input_folder, 'Mul_Original', f"{mul_photo_basename}.tif")
                pan_original_path = os.path.join(input_folder, 'Pan_Original', f"{pan_photo_basename}.tif")

                # shp_to_gpkg(mul_shp_path, mul_gpkg_path, 4326)
                # shp_to_gpkg(pan_shp_path, shp_gpkg_path, 4326)

                # copy_tif_file(mul_tif_file, mul_original_path, nodata_value, dtype)
                # copy_tif_file(pan_tif_file, pan_original_path, nodata_value, dtype)

                # MAIN_translate_gcp_image_to_origin(mul_til_file, mul_origin_path)
                # MAIN_translate_gcp_image_to_origin(pan_til_file, pan_origin_path)


                print('----------Starting flaash') # -------------------- FLAASH Multispectral image
                mul_flaash_image_path = os.path.join(root_folder_path,'Mul_FLAASH', f"{mul_photo_basename}_FLAASH.dat")
                mul_flaash_params_path = os.path.join(root_folder_path,'Mul_FLAASH', f"{mul_photo_basename}_FLAASH_Params.txt")

                # Create flaash params
                # flaash_params, output_params_path = create_flaash_params(wsl_to_windows_path(mul_tif_file), wsl_to_windows_path(mul_flaash_image_path), mul_flaash_params_path, params_overrides_scene, mul_params_overrides_photo, mul_imd_data, mul_shp_path, dem_file_path)
                # print(flaash_params, output_params_path)

                # Normal run flaash
                # MAIN_run_flaash(flaash_params, output_params_path, envi_engine, mul_flaash_image_path) # -------------------- RUN
                # convert_dat_to_tif(input_image_path=mul_flaash_dat_path, output_image_path=mul_flaash_path, mask=mul_gpkg_path, nodata_value=-9999, delete_dat_file=False)

                # Test flaash
                # if True:
                #     all_output_paths= []
                #     # Create test flaash params (optional)
                #     test_flaash_params_array = create_test_flaash_params(flaash_params, output_params_path)
                #     # print(test_flaash_params_array)
                #
                #     # Run test flaash (optional)
                #     for (test_params, test_output_params_path) in tqdm(test_flaash_params_array, desc="FLAASH"):
                #         print(test_params, test_output_params_path)
                #         MAIN_run_flaash(test_params, test_output_params_path, envi_engine)
                #         all_output_paths.append(test_params['OUTPUT_RASTER_URI'])
                #     # Or
                #     all_output_paths = parallel_flaash(test_flaash_params_array, envi_engine, max_workers=6)
                #
                #     for output_path in all_output_paths:
                #         plot_pixel_distribution(windows_to_wsl_path(output_path), save_plot_path=windows_to_wsl_path(os.path.splitext(output_path)[0] + '_plot.jpeg'))





                print('----------Starting panchomatic convertion to radiance')
                pan_radiance_path = os.path.join(root_folder_path, f'Pan_Radiance', f"{pan_photo_basename}_Radiance.tif")
                # These values come from Maxar and are updated regularly, update these values with the ones from their website
                radiance_gain = 0.130354235325
                radiance_offset = 5.505

                convert_dn_to_radiance(pan_tif_file, pan_radiance_path, radiance_gain, radiance_offset)
                # convert_dn_to_radiance_with_envi(wsl_to_windows_path(pan_tif_file), wsl_to_windows_path(pan_radiance_path), envi_engine, pan_radiance_path)





                print('----------Starting create gcp') # Step 2 -------------------- Create GCPs for image
                gcp_folder_name = 'OrthoFromUSGSLidar'
                if not os.path.exists(os.path.join(root_folder_path, gcp_folder_name)): os.makedirs(os.path.join(root_folder_path, gcp_folder_name))
                mul_origin_path = os.path.join(root_folder_path,f'Mul_Origin', f"{mul_photo_basename}_Origin.tif")
                pan_origin_path = os.path.join(root_folder_path,f'Pan_Origin', f"{pan_photo_basename}_Origin.tif")
                qgis_gcp_file_path = os.path.join(root_folder_path, gcp_folder_name, f'{mul_photo_basename}.TIF.points')
                mul_gcp_geojson_file_path = os.path.join(root_folder_path, gcp_folder_name, f'{mul_photo_basename}_points.geojson')
                pan_gcp_geojson_file_path = os.path.join(root_folder_path, gcp_folder_name, f'{pan_photo_basename}_points.geojson')

                # *Manual create GCPs in QGIS*; if a reference image is available they can be automatically created; Create GCPs in QGIS and them in this folder with the name f'{mul_origin_path}.points';

                # Converts QGIS '.points' file to Orthority GCP geojon format for refinement for mul and pan
                qgis_gcps_to_geojson(mul_tif_file, qgis_gcp_file_path, os.path.basename(mul_flaash_image_path), dem_file_path, mul_gcp_geojson_file_path, force_positive_pixel_values=True)
                qgis_gcps_to_geojson(pan_tif_file, qgis_gcp_file_path, os.path.basename(pan_radiance_path), dem_file_path, pan_gcp_geojson_file_path, force_positive_pixel_values=True)






                print('----------Starting orthorectify multispectral image') # Step 3 -------------------- Orthorectify Mulispectral Images
                # mul_ortho_default_rpc_model_path = os.path.join(root_folder_path, f'Mul_OrthoFromDefaultRPC', f"{mul_photo_basename}_OrthoFromDefaultRPC.tif")
                mul_flaash_ortho_path = os.path.join(root_folder_path, f'Mul_FLAASH_{gcp_folder_name}', f"{mul_photo_basename}_FLAASH_{gcp_folder_name}.tif")


                # Helpful command to put inpud dem into WGS84 vertical datum: gdalwarp -s_srs "+proj=longlat +datum=WGS84 +no_defs +geoidgrids=/mnt/c/Users/admin/Downloads/usa_geoid2012b/usa_geoid2012/g2012a_hawaii.gtx" -t_srs "+proj=longlat +datum=WGS84 +no_def" /mnt/d/demwgs84.tif /mnt/d/demwgs84_VWGS84.tif
                # Orthorectification info: https://up42.com/blog/how-to-perform-orthorectification-a-practical-guide
                # Ensure that the dem is in elipsoidal height, may need to convert it with a datum. proj geiods are here: https://download.osgeo.org/proj/vdatum/
                MAIN_gcp_refined_rpc_orthorectification(mul_flaash_image_path, mul_flaash_ortho_path, mul_gcp_geojson_file_path, dem_file_path, epsg)#, output_nodata_value=nodata_value, dtype='int16')






                print('----------Starting orthorectify panchromatic image') # Step 4 -------------------- Orthorectify Panchromatic Images
                pan_radiance_ortho_path = os.path.join(root_folder_path, f'Pan_Radiance_{gcp_folder_name}', f"{pan_photo_basename}_Radiance_{gcp_folder_name}.tif")
                pan_ortho_default_rpc_model_path = os.path.join(root_folder_path, f'Pan_OrthoFromDefaultRPC', f"{pan_photo_basename}_OrthoFromDefaultRPC.tif")

                MAIN_gcp_refined_rpc_orthorectification(pan_radiance_path, pan_radiance_ortho_path, pan_gcp_geojson_file_path, dem_file_path, epsg)#, output_nodata_value=nodata_value, dtype='int16')





                print('----------Starting pansharpen multispectral image') # Step 5 -------------------- Pansharp Mulispectral Images with Panchromatic Images
                # From: https://doi.org/10.5194/isprsarchives-XL-1-W1-239-2013
                mul_pansharp_path = os.path.join(root_folder_path, f'Mul_FLAASH_{gcp_folder_name}_Pansharp', f"{mul_photo_basename}_FLAASH_{gcp_folder_name}_Pansharps.tif")

                MAIN_pansharpen_image(mul_flaash_ortho_path, pan_radiance_ortho_path, mul_pansharp_path) # -------------------- RUN
                # MAIN_pansharpen_image_with_envi(wsl_to_windows_path(mul_flaash_ortho_path), wsl_to_windows_path(pan_radiance_ortho_path), wsl_to_windows_path(mul_pansharp_path), envi_engine, output_image_path_to_delete=mul_pansharp_path)






                # ------------------- Clean up dataset
                mul_cleaned_path = os.path.join(root_file_path, 'Cleaned', f'{mul_photo_basename}_Cleaned.tif')

                # replace_band_continuous_values_in_largest_segment(mul_pansharp_path, )

    print('Done with main process iamgery')

def main_match_imagery():
    # -------------------- Global Histogram Match Mulispectral Images
    input_folder = '/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445319010/Mul_FLAASH_OrthoFromUSGSLidar_Pansharp'
    global_folder = '/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages//GlobalMatch'
    output_global_basename = "_GlobalMatch"
    custom_mean_factor = 3 # Defualt 1; 3 works well sometimes
    custom_std_factor = 1 # Defualt 1

    input_image_paths_array = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.lower().endswith('.tif')]
    process_global_histogram_matching(input_image_paths_array, global_folder, output_global_basename, custom_mean_factor, custom_std_factor)






    # -------------------- Local Histogram Match Mulispectral Images
    input_image_paths_array = [os.path.join(f'{global_folder}/images', f) for f in os.listdir(f'{global_folder}/images') if f.lower().endswith('.tif')]
    local_folder = '/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445319010/LocalMatch'
    output_local_basename = "_LocalMatch"

    process_local_histogram_matching(
        input_image_paths_array,
        local_folder,
        output_local_basename,
        target_blocks_per_image = 100,
        global_nodata_value=-9999,
    )





    # if False: # -------------------- Mosaic Mulispectral Images



    print('Done with main match imagery')


def run_flaash_wrapper(args):
    """
Wrapper function (needed because executors only picklable callables).
Returns the 'OUTPUT_RASTER_URI' so we can collect it later.
    """
    test_params, test_output_params_path, envi_engine = args
    # Call your original function
    MAIN_run_flaash(test_params, test_output_params_path, envi_engine)
    return test_params["OUTPUT_RASTER_URI"]

def parallel_flaash(test_flaash_params_array, envi_engine, max_workers=4):
    """
Executes run_flaash in parallel on the items in test_flaash_params_array.
    """
    all_output_paths = []

    # Prepare the arguments for each parallel call
    tasks = [
        (test_params, test_output_params_path, envi_engine)
        for test_params, test_output_params_path in test_flaash_params_array
    ]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_flaash_wrapper, task) for task in tasks]

        # Use tqdm to display progress
        for future in tqdm(as_completed(futures), total=len(futures), desc="FLAASH"):
            # Get the result (the OUTPUT_RASTER_URI)
            output_uri = future.result()
            all_output_paths.append(output_uri)

    return all_output_paths


def change_last_folder_and_add_suffix(filepath, suffix):
    # Normalize the path to remove any trailing slashes
    filepath = os.path.normpath(filepath)

    # Split the path into directory and file components
    dir_path, file_name = os.path.split(filepath)

    if not dir_path:
        raise ValueError("The provided path does not contain a directory.")

    # Split the directory path into head and last folder
    head, last_folder = os.path.split(dir_path)

    if not last_folder:
        raise ValueError("The provided path does not contain a folder to replace.")

    # Replace the last folder with the suffix
    new_dir = os.path.join(head, suffix)

    # Create the new directory if it doesn't exist
    os.makedirs(new_dir, exist_ok=True)

    # Split the file name into name and extension
    name, ext = os.path.splitext(file_name)

    # Append the suffix to the file name
    new_file_name = f"{name}_{suffix}{ext}"

    # Join the new directory and new file name to form the new path
    new_path = os.path.join(new_dir, new_file_name)

    return new_path

def wsl_to_windows_path(path):
    wsl_pattern = r"^/mnt/([a-zA-Z])/(.*)$"
    windows_path = re.sub(wsl_pattern, r"\1:\\\2", path).replace("/", "\\")
    return windows_path

def windows_to_wsl_path(path):
    windows_pattern = r"^([a-zA-Z]):\\(.*)$"
    wsl_path = re.sub(windows_pattern, r"/mnt/\1/\2", path).replace("\\", "/")
    return wsl_path

def find_roots(input_folder):
    root_paths = []

    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.startswith("Root"):
                root_file_path = os.path.join(root, file)  # Full path to the "Root" file
                root_folder_path = os.path.dirname(root_file_path)  # Folder containing the "Root" file

                # Add the pair [root_folder_path, root_file_path] if not already in the list
                if [root_folder_path, root_file_path] not in root_paths:
                    root_paths.append([root_folder_path, root_file_path])

    print("Input array item scene count(Root files):", len(root_paths))
    return root_paths


def find_files(root_folder_path, root_file_path):
    """
    Find files with specific extensions and add their details to a JSON file.

    Parameters:
        root_folder_path (str): The root folder to search for files.
        root_file_path (str): Path to a JSON file where found files will be added.

    Returns:
        list: A list of found files with their details.
    """
    pattern = r'^(\d{2}[A-Z]{3}\d{2})(\d{6})-[A-Z0-9]+-(\d{12}).*$'
    replacement = r'\1_\2_\3'


    found_files_by_basename = {}

    # --------------- Find files in the default root folder
    for root, dirs, files in os.walk(root_folder_path):
        for file in files:
            if file.endswith(".IMD") and "M1BS" in file:
                mul_photo_basename = os.path.splitext(file)[0]
                photo_basename = re.sub(pattern, replacement, mul_photo_basename)
                if photo_basename not in found_files_by_basename:
                    found_files_by_basename[photo_basename] = {}
                found_files_by_basename[photo_basename].update({
                    'mul_imd_file': os.path.join(root, file),
                    'mul_tif_file': os.path.join(root, mul_photo_basename + ".TIF") if os.path.exists(os.path.join(root, mul_photo_basename + ".TIF")) else None,
                    'mul_photo_basename': mul_photo_basename,
                    'mul_scene_basename': os.path.basename(root_folder_path),
                    'mul_til_file': os.path.join(root, mul_photo_basename + ".TIL") if os.path.exists(os.path.join(root, mul_photo_basename + ".TIL")) else None,
                    'mul_photo_info_file': os.path.join(root, "PhotoInfo.txt") if os.path.exists(os.path.join(root, "PhotoInfo.txt")) else None,
                    'root_folder_path': root_folder_path,
                    'root_file_path': os.path.join(root_folder_path, "Root_WV03.txt")
                })


            if file.endswith(".IMD") and "P1BS" in file:
                pan_photo_basename = os.path.splitext(file)[0]
                photo_basename = re.sub(pattern, replacement, pan_photo_basename)
                if photo_basename not in found_files_by_basename:
                    found_files_by_basename[photo_basename] = {}
                found_files_by_basename[photo_basename].update({
                    'pan_imd_file': os.path.join(root, file),
                    'pan_tif_file': os.path.join(root, pan_photo_basename + ".TIF") if os.path.exists(os.path.join(root, pan_photo_basename + ".TIF")) else None,
                    'pan_photo_basename': pan_photo_basename,
                    'pan_scene_basename': os.path.basename(root_folder_path),
                    'pan_til_file': os.path.join(root, pan_photo_basename + ".TIL") if os.path.exists(os.path.join(root, pan_photo_basename + ".TIL")) else None,
                    'pan_photo_info_file': os.path.join(root, "PhotoInfo.txt") if os.path.exists(os.path.join(root, "PhotoInfo.txt")) else None,
                    'root_folder_path': root_folder_path,
                    'root_file_path': os.path.join(root_folder_path, "Root_WV03.txt")
                })

            if file.endswith(".shp") and "M1BS" in file:
                mul_gis_photo_basename = os.path.splitext(file)[0]
                photo_basename = re.sub(pattern, replacement, mul_gis_photo_basename)
                if photo_basename not in found_files_by_basename:
                    found_files_by_basename[photo_basename] = {}
                found_files_by_basename[photo_basename].update({
                    'mul_shp_file': os.path.join(root, file),
                })

            if file.endswith(".shp") and "P1BS" in file:
                pan_gis_photo_basename = os.path.splitext(file)[0]
                photo_basename = re.sub(pattern, replacement, pan_gis_photo_basename)
                if photo_basename not in found_files_by_basename:
                    found_files_by_basename[photo_basename] = {}
                found_files_by_basename[photo_basename].update({
                    'pan_shp_file': os.path.join(root, file),
                })

    # Filter by filter_basenames
    if filter_basenames:
        found_files_by_basename = {
            k: v for k, v in found_files_by_basename.items()
            if v.get('mul_photo_basename') in filter_basenames or v.get('pan_photo_basename') in filter_basenames
        }


    return found_files_by_basename




def find_subfolder_files(root_folder_path, subfolders_array, filter_basenames=None):
    """
    Finds and processes subfolders directly under the root folder that match the given names.

    Args:
        root_folder_path (str): The root directory to search for subfolders.
        subfolders_array (list of str): List of subfolder names to look for.
        filter_basenames (list of str, optional): List of basenames to filter files by.

    Returns:
        dict: A dictionary where keys are matching subfolder names,
              and values are lists of '.dat' files with their full paths found in those subfolders.
    """
    found_subfolders = {}

    # Iterate over all subfolders in the root folder
    for subfolder_name in os.listdir(root_folder_path):
        subfolder_path = os.path.join(root_folder_path, subfolder_name)

        # Check if it's a directory and matches the subfolders array
        if os.path.isdir(subfolder_path) and subfolder_name in subfolders_array:
            # List all .dat files in the subfolder with full paths
            subfolder_files = [
                os.path.join(subfolder_path, f) for f in os.listdir(subfolder_path)
                if os.path.isfile(os.path.join(subfolder_path, f)) and f.endswith('.dat')
            ]

            # Filter files by basenames if filter_basenames is specified
            if filter_basenames:
                subfolder_files = [
                    f for f in subfolder_files
                    if any(basename in os.path.basename(f) for basename in filter_basenames)
                ]

            # Add to the found subfolders dictionary
            found_subfolders[subfolder_name] = subfolder_files

    return found_subfolders




import re
import json

def get_metadata_from_files(root_file_path, imd_file, photo_basename):
    """
Extract metadata from the provided files.

Args:
root_file_path (str): Path to the root file containing scene-level overrides.
photo_info_file (str): Path to the photo info file containing photo-level overrides.
imd_file (str): Path to the IMD file containing metadata.

Returns:
dict: ParamsOverridesPerScene from root file.
dict: ParamsOverridesPerPhoto from photo info file.
dict: Metadata extracted from the IMD file.
    """

    # 1. Extract scene-level overrides
    params_overrides_scene = {}
    # 2. Extract photo-level overrides
    params_overrides_photo = {}

    if root_file_path:
        with open(root_file_path, 'r') as file:
            data = json.load(file)
            params_overrides_scene = data.get("ParamsOverridesPerScene", {})

            # Locate the dictionary of all photo overrides
            all_photos_overrides = data.get("ParamsOverridesPerPhoto", {})

            # Extract overrides for the specific photo_basename
            params_overrides_photo = all_photos_overrides.get(photo_basename, {})

    # Extract metadata from IMD file
    imd_data = {}
    key_mapping = {
        "SOLAR_AZIMUTH": "meanSunAz",
        "SOLAR_ZENITH": "meanSunEl",
        "LOS_AZIMUTH": "meanSatAz",
        "LOS_ZENITH": "meanOffNadirViewAngle",
        'x_res': 'meanProductRowGSD',
        'y_res': 'meanProductColGSD',
    }
    regex_pattern = r"(?P<key>{})\s*=\s*(?P<value>-?\d+(\.\d+)?)".format("|".join(key_mapping.values()))

    if imd_file:
        with open(imd_file, "r") as file:
            imd_content = file.read()

            matches = re.finditer(regex_pattern, imd_content)
            for match in matches:
                imd_key = match.group("key")
                value = float(match.group("value"))

                for python_key, imd_key_match in key_mapping.items():
                    if imd_key == imd_key_match:
                        if python_key == "SOLAR_ZENITH":
                            imd_data[python_key] = 90 - value
                        elif python_key == "LOS_ZENITH":
                            imd_data[python_key] = 180 - value
                        else:
                            imd_data[python_key] = value
                        break

    return params_overrides_scene, params_overrides_photo, imd_data


from osgeo import ogr, osr

def shp_to_gpkg(input_shp_path, output_gpkg_path, override_projection_epsg=None):
    """
    Reads a shapefile (which might not have a .prj),
    optionally forces a new EPSG code on it,
    then saves it as a GeoPackage (no in-memory layers used).

    Parameters
    ----------
    input_shp_path : str
    Path to the input shapefile (possibly missing a .prj).
    output_gpkg_path : str
    Path to the output GeoPackage file.
    override_projection_epsg : int, optional
    If provided, the layer's CRS is forced to this EPSG without reprojecting geometries.
    """

    # Open the input shapefile (read-only)
    shp_driver = ogr.GetDriverByName("ESRI Shapefile")
    in_ds = shp_driver.Open(input_shp_path, 0)  # 0 = read-only
    if not in_ds:
        raise FileNotFoundError(f"Could not open Shapefile: {input_shp_path}")

    in_layer = in_ds.GetLayer()
    if not in_layer:
        raise RuntimeError("Could not get layer from Shapefile.")

    # Determine the layer's current spatial reference
    in_srs = in_layer.GetSpatialRef()

    # If user wants to override the CRS, just assign that EPSG
    if override_projection_epsg is not None:
        out_srs = osr.SpatialReference()
        out_srs.ImportFromEPSG(override_projection_epsg)
    else:
        out_srs = in_srs  # keep the original (which could be None)

    # Prepare output: if the GPKG already exists, remove it
    gpkg_driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(output_gpkg_path):
        os.remove(output_gpkg_path)

    if not os.path.exists(os.path.dirname(output_gpkg_path)):
        os.makedirs(os.path.dirname(output_gpkg_path))

    out_ds = gpkg_driver.CreateDataSource(output_gpkg_path)

    # Create the output layer with the (possibly overridden) CRS
    # Note: if out_srs is None, it will simply have "unknown" CRS in the GeoPackage
    out_layer = out_ds.CreateLayer(
        name="layer",
        srs=out_srs,
        geom_type=in_layer.GetGeomType()
    )

    # Copy fields from the input layer
    in_layer_defn = in_layer.GetLayerDefn()
    for i in range(in_layer_defn.GetFieldCount()):
        field_defn = in_layer_defn.GetFieldDefn(i)
        out_layer.CreateField(field_defn)

    # Copy features from input to output
    in_layer.ResetReading()
    for in_feature in in_layer:
        out_feature = ogr.Feature(out_layer.GetLayerDefn())

        # Copy geometry by cloning (no reprojection)
        geom = in_feature.GetGeometryRef()
        if geom is not None:
            out_feature.SetGeometry(geom.Clone())

        # Copy field attributes
        for i in range(in_layer_defn.GetFieldCount()):
            out_feature.SetField(
                in_layer_defn.GetFieldDefn(i).GetNameRef(),
                in_feature.GetField(i)
            )

        out_layer.CreateFeature(out_feature)
        out_feature = None

    # Cleanup
    out_ds = None
    in_ds = None

    if override_projection_epsg:
        print(f"Assigned EPSG:{override_projection_epsg} to output layer {output_gpkg_path}")
    else:
        print(f"No CRS override. Output saved to '{output_gpkg_path}'.")

def convert_dat_to_tif(input_image_path, output_image_path, mask, nodata_value, delete_dat_file=False):
    """Convert a .dat raster to .tif, apply mask, set NoData value, and clean up files."""
    gdal.Warp(
        output_image_path, input_image_path, format="GTiff", dstNodata=nodata_value, srcNodata=nodata_value
    )
    # Delete the original file and associated .hdr and .enp files if output was created successfully
    if os.path.exists(output_image_path) and delete_dat_file:
        os.remove(input_image_path)
        base_name, _ = os.path.splitext(input_image_path)
        for ext in [".hdr", ".enp"]:
            extra_file = base_name + ext
            if os.path.exists(extra_file):
                os.remove(extra_file)

def copy_tif_file(input_image_path, output_image_path, output_nodata_value=None, output_dtype=None):
    """
    Copies a GeoTIFF file to a new location without altering its properties using gdal.Translate.
    Optionally sets a NoData value and allows specifying a GDAL dtype.

    Parameters
    ----------
    input_image_path : str
    Path to the input GeoTIFF file.
    output_image_path : str
    Path to the output GeoTIFF file.
    output_nodata_value : float, optional
    If set, applies the NoData value to all bands.
    output_dtype : str, optional
    GDAL data type (e.g., 'Int16', 'Float32'). If None, keeps the input dtype.
    """
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_image_path), exist_ok=True)

    # Open the input dataset
    src_ds = gdal.Open(input_image_path, gdal.GA_ReadOnly)
    if not src_ds:
        raise ValueError(f"Unable to open input image: {input_image_path}")

    # Define options for gdal.Translate
    translate_options = gdal.TranslateOptions(
        outputType=gdal.GetDataTypeByName(output_dtype) if output_dtype else None,  # Convert output_dtype if specified
        noData=output_nodata_value if output_nodata_value is not None else None,  # Apply NoData if set
        format="GTiff",  # Ensure output format is GeoTIFF
        metadataOptions=["COPY_SRC_METADATA"],  # Preserve metadata
    )

    # Translate (copy) the image
    gdal.Translate(output_image_path, src_ds, options=translate_options)

    # Cleanup
    src_ds = None

    print(f"Image copied to: {output_image_path} with dtype {output_dtype or 'unchanged'}")

from osgeo import gdal, gdal_array
import numpy as np

def spectral_scale_image(input_image_path, output_image_path, input_image_scale, output_image_scale):
    """
Scales an image from one range to another while preserving NoData values and copying all metadata.

Parameters:
input_image_path (str): Path to the input image file.
output_image_path (str): Path to the output image file.
input_image_scale (tuple): Tuple (min, max) defining the scale of the input image.
output_image_scale (tuple): Tuple (min, max) defining the scale of the output image.
    """
    # Open input image
    dataset = gdal.Open(input_image_path, gdal.GA_ReadOnly)
    if dataset is None:
        raise ValueError("Could not open input image")

    # Get image metadata
    metadata = dataset.GetMetadata()
    nodata_value = dataset.GetRasterBand(1).GetNoDataValue()

    # Read the raster bands
    bands = [dataset.GetRasterBand(i+1).ReadAsArray().astype(np.float32) for i in range(dataset.RasterCount)]

    # Scale function
    def scale_band(band):
        mask = (band == nodata_value)
        scaled_band = (band - input_image_scale[0]) / (input_image_scale[1] - input_image_scale[0])
        scaled_band = scaled_band * (output_image_scale[1] - output_image_scale[0]) + output_image_scale[0]
        scaled_band[mask] = nodata_value  # Restore NoData values
        return scaled_band

    # Apply scaling to each band
    scaled_bands = [scale_band(band) for band in bands]

    # Create output image
    driver = gdal.GetDriverByName("GTiff")
    out_dataset = driver.Create(output_image_path, dataset.RasterXSize, dataset.RasterYSize, dataset.RasterCount, gdal.GDT_Float32)

    # Copy geotransform and projection
    out_dataset.SetGeoTransform(dataset.GetGeoTransform())
    out_dataset.SetProjection(dataset.GetProjection())
    out_dataset.SetMetadata(metadata)

    # Write scaled bands
    for i, scaled_band in enumerate(scaled_bands):
        out_band = out_dataset.GetRasterBand(i+1)
        out_band.WriteArray(scaled_band)
        out_band.SetNoDataValue(nodata_value)

    # Close datasets
    out_dataset = None
    dataset = None


from osgeo import gdal, gdal_array

from osgeo import gdal
import os

def convert_dn_to_radiance(input_image_path, output_image_path, gain, offset):
    ds = gdal.Open(input_image_path, gdal.GA_ReadOnly)
    if not ds:
        raise RuntimeError(f"Failed to open {input_image_path}")

    os.makedirs(os.path.dirname(output_image_path), exist_ok=True)

    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(output_image_path, ds.RasterXSize, ds.RasterYSize, ds.RasterCount, gdal.GDT_Float32)
    out_ds.SetGeoTransform(ds.GetGeoTransform())
    out_ds.SetProjection(ds.GetProjection())
    out_ds.SetMetadata(ds.GetMetadata())

    # Copy RPC metadata
    rpc_metadata = ds.GetMetadata("RPC")
    if rpc_metadata:
        out_ds.SetMetadata(rpc_metadata, "RPC")

    for band in range(1, ds.RasterCount + 1):
        src_band = ds.GetRasterBand(band)
        data = src_band.ReadAsArray()
        radiance = data * gain + offset

        out_band = out_ds.GetRasterBand(band)
        out_band.WriteArray(radiance)

        # Copy NoData value
        nodata = src_band.GetNoDataValue()
        if nodata is not None:
            out_band.SetNoDataValue(nodata)

    out_ds.FlushCache()
    del ds, out_ds






main_process_imagery(input_folders_array)
# main_match_imagery()




