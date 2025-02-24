# Apply FLAASH to multiple photos in an input folder and output multiple atmosphericlly corrected images

# Update these parameters in the code as needed:
    # parameters
        # Check weather here: https://www.wunderground.com/history; visibility only ever gets to 10miles here but is presumably higher.
    # input_folder
        # Input path in WINDOWS format (but add extra "\" to parse them correctly or preface string with 'r')
    # output_folder
        # Input path in WINDOWS format (but add extra "\" to parse them correctly or preface string with 'r')
        # Either manually create output folder or it will be created from the output folder path; it should be in the following format: 20171208_36cm_WV03_BAB_FLAASH
    # key_mapping
    # math(math in CreateFLAASHParameters function)

# --------------------Define inputs and outputs manually
# input_folders_array = [
#     "/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445319010",
#     "/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445318010"
# ]

# output_folders_array = [
#     "/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445319010/20171208_36cm_WV03_BAB_016445319010_FLAASH",
#     "/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445318010/20171208_36cm_WV03_BAB_016445318010_FLAASH"
# ]

# --------------------Start script
import os
from envipyengine import Engine
from pprint import pprint
import envipyengine.config
from concurrent.futures import ProcessPoolExecutor, as_completed
import re
from tqdm import tqdm
import json
from helper_functions import get_image_largest_value
# Start the ENVI engine
envipyengine.config.set('engine', "/mnt/c/Program Files/Harris/ENVI57/IDL89/bin/bin.x86_64/taskengine.exe")
envi_engine = Engine('ENVI')
envi_engine.tasks()

# Resources to find currect parameters:
    # Reference ENVI GUI for better understanding of settings
    # https://www.nv5geospatialsoftware.com/docs/Flaash.html
    # https://www.nv5geospatialsoftware.com/Support/Self-Help-Tools/Help-Articles/Help-Articles-Detail/ArtMID/10220/ArticleID/24080/Atmospheric-correction-FLAASH-vs-QUAC-which-one-should-I-use
    # http://essay.utwente.nl/83442/1/davaadorj.pdf
    # https://www.nv5geospatialsoftware.com/docs/Flaash.html
    # https://www.nv5geospatialsoftware.com/portals/0/pdfs/envi/Flaash_Module.pdf
    # https://www.sciencedirect.com/science/article/pii/S0924271621000095
    # https://dg-cms-uploads-production.s3.amazonaws.com/uploads/document/file/106/ISD_External.pdf
    # https://www.mdpi.com/1424-8220/16/10/1624
    # worldview.earthdata.nasa.gov (product: 'MOD05_L2' (Water Vapor))
    # https://github.com/dawhite/MCTK (MODIS Conversion Toolkit (MCTK))
    # https://www.wunderground.com/history


def create_flaash_params(input_image_path, output_image_path, output_params_path, params_overrides_scene, params_overrides_photo, imd_data, mask=None, dem_file=None):
    """
    Create FLAASH parameters based on the provided inputs.

    Args:
        input_image_path (str): Path to the input image.
        output_image_path (str): Path to the output image.
        output_params_path (str): Path to the output parameters.
        params_overrides_scene (dict): Scene-specific overrides.
        params_overrides_photo (dict): Photo-specific overrides.
        imd_data (dict): Metadata extracted from the IMD file.

    Returns:
        dict: Final FLAASH parameters.
        str: Path to the output parameters file.
    """

    largest_elevation = get_image_largest_value(dem_file, mask, override_mask_crs_epsg=4326)/1000
    flaash_params = {
        'INPUT_RASTER': {
            'url': input_image_path,
            'factory': 'URLRaster'
        },
        'SENSOR_TYPE': None,
        'INPUT_SCALE': None,
        'OUTPUT_SCALE': None,
        'CALIBRATION_FILE': None,
        'CALIBRATION_FORMAT': None,
        'CALIBRATION_UNITS': None,
        'LAT_LONG': None,
        'SENSOR_ALTITUDE': None,
        'DATE_TIME': None,
        'USE_ADJACENCY': None,
        'DEFAULT_VISIBILITY': None,
        'USE_POLISHING': None,
        'POLISHING_RESOLUTION': None,
        'SENSOR_AUTOCALIBRATION': None,
        'SENSOR_CAL_PRECISION': None,
        'SENSOR_CAL_FEATURE_LIST': None,
        'GROUND_ELEVATION': largest_elevation,
        'SOLAR_AZIMUTH': imd_data.get("SOLAR_AZIMUTH"),
        'SOLAR_ZENITH': imd_data.get("SOLAR_ZENITH"),
        'LOS_AZIMUTH': imd_data.get("LOS_AZIMUTH"),
        'LOS_ZENITH': imd_data.get("LOS_ZENITH"),
        'IFOV': None,
        'MODTRAN_ATM': 'Mid-Latitude Summer',
        'MODTRAN_AER': 'Maritime',
        'MODTRAN_RES': 5.0,
        'MODTRAN_MSCAT': "DISORT",
        'CO2_MIXING': None,
        'WATER_ABS_CHOICE': None,
        'WATER_MULT': None,
        'WATER_VAPOR_PRESET': None,
        'USE_AEROSOL': 'Disabled',
        'AEROSOL_SCALE_HT': None,
        'AER_BAND_RATIO': 0.5,
        'AER_BAND_WAVL': None,
        'AER_REFERENCE_VALUE': None,
        'AER_REFERENCE_PIXEL': None,
        'AER_BANDLOW_WAVL': 425,
        'AER_BANDLOW_MAXREFL': None,
        'AER_BANDHIGH_WAVL': 660,
        'AER_BANDHIGH_MAXREFL': 0.2,
        'CLOUD_RASTER_URI': None,
        'WATER_RASTER_URI': None,
        'OUTPUT_RASTER_URI': output_image_path,
        'CLOUD_RASTER': None,
        'WATER_RASTER': None,
        'OUTPUT_RASTER': None,
    }

    for key, value in params_overrides_scene.items():
        flaash_params[key] = value

    for key, value in params_overrides_photo.items():
        flaash_params[key] = value

    final_flaash_params = {key: value for key, value in flaash_params.items() if value is not None}
    return final_flaash_params, output_params_path





def run_flaash(flaash_params, output_params_path, envi_engine, output_image_path_to_delete=None):
    print('Processing photo with params:', flaash_params)
    if output_image_path_to_delete:
        try:
            if os.path.exists(output_image_path_to_delete):
                os.remove(output_image_path_to_delete)
                print(f"Deleted existing output image: {output_image_path_to_delete}")
        except Exception as e:
            print(f"Error deleting output image {output_image_path_to_delete}: {e}")

    try:
        task = envi_engine.task('FLAASH')
        result = task.execute(flaash_params)
        with open(output_params_path, 'w') as file: file.write(str(flaash_params))
        print(f"Photo complete")

    except Exception as e:
        print(f"Error processing: {flaash_params}: {e}")

def run_flaash_wrapper(args):
    """
    Wrapper function (needed because executors only picklable callables).
    Returns the 'OUTPUT_RASTER_URI' so we can collect it later.
    """
    test_params, test_output_params_path, envi_engine = args
    # Call your original function
    run_flaash(test_params, test_output_params_path, envi_engine)
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
