import os
import re
import json
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple

_WV_BAND_ORDER = [
    "BAND_C",
    "BAND_B",
    "BAND_G",
    "BAND_Y",
    "BAND_R",
    "BAND_RE",
    "BAND_N",
    "BAND_N2",
]

_WV_RADIOMETRIC_GAIN_OFFSET_2018V0: Dict[str, Dict[str, Tuple[float, float]]] = {
    # Source: Maxar radiometric calibration update tables (2018v0) for WV2/WV3 VNIR bands.
    # Formula: radiance = gain * (DN * absCalFactor/effectiveBandwidth) + offset
    "WV02": {
        "BAND_C": (0.977, -5.552),
        "BAND_B": (0.987, -4.443),
        "BAND_G": (0.987, -2.523),
        "BAND_Y": (0.987, -4.929),
        "BAND_R": (0.987, -3.914),
        "BAND_RE": (0.987, -5.675),
        "BAND_N": (0.987, -1.546),
        "BAND_N2": (0.987, -3.564),
    },
    "WV03": {
        "BAND_C": (0.905, -8.604),
        "BAND_B": (0.940, -5.809),
        "BAND_G": (0.938, -4.996),
        "BAND_Y": (0.962, -3.649),
        "BAND_R": (0.964, -3.021),
        "BAND_RE": (1.000, -4.521),
        "BAND_N": (0.961, -5.522),
        "BAND_N2": (0.978, -2.992),
    },
}


def find_roots(
    input_folder
    ):
    """Discover scene root folders and optional Root* override files.

    Discovery behavior:
    1. If Root* files are present, use those folder/file pairs.
    2. Otherwise, auto-discover scene folders from Maxar bundle structure
       (M1BS/P1BS IMD files) and return root file as ``None``.
    """

    root_paths = []

    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.startswith("Root"):
                root_file_path = os.path.join(root, file)  # Full path to the "Root" file
                root_folder_path = os.path.dirname(root_file_path)  # Folder containing the "Root" file

                # Add the pair [root_folder_path, root_file_path] if not already in the list
                if [root_folder_path, root_file_path] not in root_paths:
                    root_paths.append([root_folder_path, root_file_path])

    if root_paths:
        print("Input array item scene count(Root files):", len(root_paths))
        return root_paths

    # Fallback discovery for raw Maxar data without Root_WV.txt.
    scene_flags: Dict[str, Dict[str, bool]] = {}
    for root, _, files in os.walk(input_folder):
        has_mul = any(f.endswith(".IMD") and "M1BS" in f for f in files)
        has_pan = any(f.endswith(".IMD") and "P1BS" in f for f in files)
        if not (has_mul or has_pan):
            continue

        scene_root = root
        base = os.path.basename(root).upper()
        if base.endswith("_MUL") or base.endswith("_PAN"):
            scene_root = os.path.dirname(root)

        flags = scene_flags.setdefault(scene_root, {"mul": False, "pan": False})
        flags["mul"] = flags["mul"] or has_mul
        flags["pan"] = flags["pan"] or has_pan

    for scene_root, flags in sorted(scene_flags.items()):
        if flags["mul"] and flags["pan"]:
            root_paths.append([scene_root, None])

    print("Input array item scene count(auto-discovered roots):", len(root_paths))
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
                    'root_file_path': root_file_path,
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
                    'root_file_path': root_file_path,
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

    if root_file_path and os.path.exists(root_file_path):
        with open(root_file_path, 'r') as file:
            content = file.read().strip()  # Read and strip whitespace to check for emptiness
            if not content:
                params_overrides_scene = None
                params_overrides_photo = None
            else:
                try:
                    data = json.loads(content)
                except json.JSONDecodeError:
                    data = {}

                try:
                    params_overrides_scene = data.get("ParamsOverridesPerScene", None)
                except (AttributeError, TypeError):
                    params_overrides_scene = None

                try:
                    all_photos_overrides = data.get("ParamsOverridesPerPhoto", {})
                    params_overrides_photo = all_photos_overrides.get(photo_basename, None)
                except (AttributeError, TypeError):
                    params_overrides_photo = None

    # Extract metadata from IMD file
    imd_data = {}
    key_mapping = {
        "SOLAR_AZIMUTH": "meanSunAz",
        "SOLAR_ZENITH": "meanSunEl",
        "LOS_AZIMUTH": "meanSatAz",
        "LOS_ZENITH": "meanOffNadirViewAngle",
        'x_res': 'meanProductRowGSD',
        'y_res': 'meanProductColGSD',
        'product_res': 'meanProductGSD',
    }
    regex_pattern = r"(?P<key>{})\s*=\s*(?P<value>-?\d+(\.\d+)?)".format("|".join(key_mapping.values()))

    if imd_file:
        with open(imd_file, "r") as file:
            imd_content = file.read()

            # Prefer IMD acquisition timestamp for any time-dependent processing
            # (for example NASA atmosphere queries), because basename date tokens
            # can be inconsistent across products.
            acq_dt = None
            for key in ("firstLineTime", "TLCTime", "earliestAcqTime", "latestAcqTime"):
                m = re.search(
                    rf"{key}\s*=\s*(?:\"([^\"]+)\"|([0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}}T[0-9:.+-]+Z?))\s*;",
                    imd_content,
                )
                if not m:
                    continue
                raw = (m.group(1) or m.group(2) or "").strip()
                try:
                    acq_dt = datetime.fromisoformat(raw.replace("Z", "+00:00"))
                    break
                except ValueError:
                    continue
            if acq_dt is not None:
                if acq_dt.tzinfo is None:
                    acq_dt = acq_dt.replace(tzinfo=timezone.utc)
                imd_data["ACQUISITION_DATETIME_UTC"] = acq_dt.astimezone(timezone.utc).isoformat()

            matches = re.finditer(regex_pattern, imd_content)
            for match in matches:
                imd_key = match.group("key")
                value = float(match.group("value"))

                for python_key, imd_key_match in key_mapping.items():
                    if imd_key == imd_key_match:
                        if python_key == "SOLAR_ZENITH":
                            imd_data[python_key] = 90 - value
                        elif python_key == "LOS_ZENITH":
                            # FLAASH expects line-of-sight zenith in the [90, 270] convention.
                            # IMD meanOffNadirViewAngle is from nadir (0-90), so convert accordingly.
                            imd_data[python_key] = 180 - value
                        else:
                            imd_data[python_key] = value
                        break
            offnadir_match = re.search(
                r"meanOffNadirViewAngle\s*=\s*(-?\d+(\.\d+)?)",
                imd_content,
            )
            if offnadir_match:
                # Py6S Geometry.User().view_z expects 0..90 degrees from nadir.
                imd_data["VIEW_ZENITH"] = float(offnadir_match.group(1))

            satid_match = re.search(r'satId\s*=\s*"([^"]+)"', imd_content)
            sensor_id = None
            if satid_match:
                sensor_id = satid_match.group(1)
                imd_data["sensor_id"] = sensor_id

            band_pattern = re.compile(
                r"BEGIN_GROUP\s*=\s*(BAND_[A-Z0-9]+)(.*?)END_GROUP\s*=\s*\1",
                re.DOTALL,
            )
            abs_pattern = re.compile(r"absCalFactor\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")
            bw_pattern = re.compile(r"effectiveBandwidth\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")

            band_factors: Dict[str, float] = {}
            discovered_band_order: List[str] = []
            for match in band_pattern.finditer(imd_content):
                band_name = match.group(1)
                body = match.group(2)
                abs_match = abs_pattern.search(body)
                bw_match = bw_pattern.search(body)
                if not abs_match or not bw_match:
                    continue
                abs_cal = float(abs_match.group(1))
                eff_bw = float(bw_match.group(1))
                if eff_bw == 0.0:
                    continue
                discovered_band_order.append(band_name)
                band_factors[band_name] = abs_cal / eff_bw

            if band_factors and sensor_id and sensor_id.upper().startswith("WV"):
                ordered_names = [b for b in _WV_BAND_ORDER if b in band_factors]
                if not ordered_names:
                    ordered_names = [b for b in discovered_band_order if b in band_factors]
                imd_data["worldview_band_order"] = ordered_names
                imd_data["worldview_dn_to_radiance_factors"] = [band_factors[b] for b in ordered_names]

                sensor_key = sensor_id.upper()
                gain_offsets = _WV_RADIOMETRIC_GAIN_OFFSET_2018V0.get(sensor_key)
                if gain_offsets and all(b in gain_offsets for b in ordered_names):
                    gains = [gain_offsets[b][0] for b in ordered_names]
                    offsets = [gain_offsets[b][1] for b in ordered_names]
                    imd_data["worldview_radiometric_adjustment_version"] = "2018v0"
                    imd_data["worldview_dn_to_radiance_gains"] = gains
                    imd_data["worldview_dn_to_radiance_offsets"] = offsets

    return params_overrides_scene, params_overrides_photo, imd_data
