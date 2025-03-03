import os
import re
import json

def find_roots(
        input_folder
        ):
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

def find_files(
        root_folder_path,
        root_file_path,
        filter_basenames
        ):
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
                    'root_folder_path': root_folder_path,
                    'root_file_path': os.path.join(root_folder_path, "Root_WV.txt")
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
                    'root_folder_path': root_folder_path,
                    'root_file_path': os.path.join(root_folder_path, "Root_WV.txt")
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

def find_subfolder_files(
        root_folder_path,
        subfolders_array,
        filter_basenames=None
        ):
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

def get_metadata_from_files(
        root_file_path,
        imd_file,
        photo_basename
        ):
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