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
