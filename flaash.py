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

import os
import itertools

def create_test_flaash_params(flaash_params, output_params_path):
    test_params={ #'OUTPUT_SCALE': {None,2, 3, 4}, 'AER_BAND_RATIO': {None, 0.3, 4}})
                 'SENSOR_TYPE': None,
                 'INPUT_SCALE': None,
                 'OUTPUT_SCALE': None,
                 'CALIBRATION_FILE': None,
                 'CALIBRATION_FORMAT': None,
                 'CALIBRATION_UNITS': None,
                 'LAT_LONG': None,
                 'SENSOR_ALTITUDE': None,
                 'DATE_TIME': None,
                 'USE_ADJACENCY': None, # {'Disabled','Legacy exponential scattering kernel','Wavelength-dependent scattering kernel'},
                 'DEFAULT_VISIBILITY': None, # {10, 40, 100},
                 'USE_POLISHING': None, #'Disabled', 'Polish using reference materials', 'Polish using statistical detection of spectral artifacts'
                 'POLISHING_RESOLUTION': None,
                 'SENSOR_AUTOCALIBRATION': None, #{True, False},
                 'SENSOR_CAL_PRECISION': None,
                 'SENSOR_CAL_FEATURE_LIST': None,
                 'GROUND_ELEVATION': {None, 0, 1.6, 2, 3, 4, 5},
                 'SOLAR_AZIMUTH': None,
                 'SOLAR_ZENITH': None,
                 'LOS_AZIMUTH': None,
                 'LOS_ZENITH': None,
                 'IFOV': None,
                 'MODTRAN_ATM': None, # {'Tropical Atmosphere','Mid-Latitude Summer','Mid-Latitude Winter','Sub-Arctic Summer','Sub-Arctic Winter','1976 US Standard Atmosphere'},
                 'MODTRAN_AER': None, # {'No Aerosol', 'High-Visibility Rural', 'Low-Visibility Rural', 'Maritime', 'Urban', 'Tropospheric'},
                 'MODTRAN_RES': {15},
                 'MODTRAN_MSCAT': None,
                 'CO2_MIXING': None,
                 'WATER_ABS_CHOICE': None,
                 'WATER_MULT': None,
                 'WATER_VAPOR_PRESET': None, #{None, 0.89, 0.2, 1, 1.3, 3, 5}, #Should be defined per scene or even better, per photo; we want about 2.6 for PuuWaawaa on 20241208 based on MODIS data, based on Mid latitude Summer model (water vapor = 2.92 g/cm^2), 2.6=2.92*x, x=0.89
                 'USE_AEROSOL': {'Disabled'}, #DONE # Retrieval: 'Disabled', 'Automatic Selection', 'Vegetation Based Retrieval', 'Water Based Retrieval', 'Wavelength Dependent Water Based Retrieval', 'Linear Regression Retrieval'
                 'AEROSOL_SCALE_HT': {2, 1, 2.5, 0.5},
                 'AER_BAND_RATIO': {0.5, 0.2, 0.9},
                 'AER_BAND_WAVL': None,
                 'AER_REFERENCE_VALUE': None,
                 'AER_REFERENCE_PIXEL': None,
                 'AER_BANDLOW_WAVL': None, #{425}, # 425 (for single def)
                 'AER_BANDLOW_MAXREFL': None,
                 'AER_BANDHIGH_WAVL': None, # {660}, # 660 (for single def)
                 'AER_BANDHIGH_MAXREFL': None, #{0.1, 0.2, 0.5}, #0.02,
                 }

    # ----------------------- 1) Exclude OUTPUT_RASTER_URI from testable params -----------------------
    testable_params = {
        k: v
        for k, v in test_params.items()
        if k != "OUTPUT_RASTER_URI" and v is not None and len(v) > 0
    }

    # If there are no testable params, just return one combination with the original flaash_params
    if not testable_params:
        return [(flaash_params, output_params_path)]

    # ----------------------- 2) Prepare for permutations of testable params -------------------------
    param_items = list(testable_params.items())
    # Example: [("AER_BAND_RATIO", [0.5,0.4,0.3]), ("USE_AEROSOL", ["Disabled","Auto"])]

    param_names = [item[0] for item in param_items]         # e.g. ["AER_BAND_RATIO", "USE_AEROSOL"]
    param_values_lists = [item[1] for item in param_items]  # e.g. [[0.5,0.4,0.3], ["Disabled","Auto"]]

    # Cartesian product of all values => e.g. 3 x 2 = 6 combos
    all_combinations = list(itertools.product(*param_values_lists))

    # ----------------------- 3) Base output paths for Raster and Params -----------------------------
    # Raster base
    base_output_raster_uri = flaash_params.get("OUTPUT_RASTER_URI", "")
    base_raster_dir = os.path.dirname(base_output_raster_uri)
    base_raster_filename = os.path.basename(base_output_raster_uri)
    raster_name, raster_ext = os.path.splitext(base_raster_filename)

    # Params base
    base_params_dir = os.path.dirname(output_params_path)
    base_params_filename = os.path.basename(output_params_path)
    params_name, params_ext = os.path.splitext(base_params_filename)

    # ----------------------- 4) Build updated param dictionaries + param file paths -----------------
    result_list = []

    for combo in all_combinations:
        # Start with a copy of the base dictionary
        combo_params = dict(flaash_params)

        # Build a short suffix for each param-value in this combination
        short_suffix_parts = []
        for name, val in zip(param_names, combo):
            # name[:4], str(val)[:4]
            short_suffix_parts.append(f"{name[:3]}-{str(val)[:3]}")

            # Overwrite that parameter in the dictionary
            combo_params[name] = val

        # Combine them into one short suffix
        # e.g. "AER_-0.5_USE_-Disa"
        short_suffix = "_".join(short_suffix_parts)

        # Create the new OUTPUT_RASTER_URI
        new_raster_filename = f"{raster_name}_{short_suffix}{raster_ext}"
        new_raster_path = os.path.join(base_raster_dir, new_raster_filename)
        combo_params["OUTPUT_RASTER_URI"] = new_raster_path

        # Create the new output_params_path
        new_params_filename = f"{params_name}_{short_suffix}{params_ext}"
        new_params_path = os.path.join(base_params_dir, new_params_filename)

        # Append to results
        result_list.append((combo_params, new_params_path))

    return result_list
    