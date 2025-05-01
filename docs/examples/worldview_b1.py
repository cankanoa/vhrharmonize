# ------------------- Imports
import os
from vhrharmonize import(
    find_files,
    find_roots,
    get_metadata_from_files,
    shp_to_gpkg,
    get_image_largest_value,
    wsl_to_windows_path,
    qgis_gcps_to_geojson,
    gcp_refined_rpc_orthorectification,
    pansharpen_image,
    run_flaash
    )

from envipyengine import Engine
import envipyengine.config
envipyengine.config.set('engine', "/path/to/envi/taskengine.exe")
envi_engine = Engine('ENVI')
envi_engine.tasks()

# ------------------- Define variables
# Code will search these folders for all Root_WV.txt files to use as roots for each scene
input_folders_array = [
    "/folder/to/search/for/root_wv.txt files within",
    "/second/path/to/search",
]

# Only process these iamges if set, if not set, all images will be processed (use either mul or pan photo name)
filter_basenames = [
    'only_process_these_images_by_basename',
]

# Some other requires params
epsg = 4326 # or a more specific one, for example 6635 for the Big Island
nodata_value = -9999
dtype = 'int16'
calculation_dtype_precision = 'float32'

# These values come from Maxar and are updated regularly, update these values with the ones from their website if using the conversion from dn to radiance
pan_radiance_gain = 0.130354235325
pan_radiance_offset = 5.505

# Get WATER_VAPOR_PRESET from: www.worldview.earthdata.nasa.gov (product: 'MOD05_L2' (Water Vapor))
# Required to create a Root_WV.txt file at each satellite image root folder; optional: add override to params per scene and photo in the following format:
# {  "ParamsOverridesPerScene": {
#       "WATER_VAPOR_PRESET": 0.89
#   },
#   "ParamsOverridesPerPhoto": {
#       "17DEC08211758-M1BS-016445319010_01_P003": {
#           "Key": "value"
#       }
#   }
# }

# DEM file needs to be in WGS84 elipsoidal height and cover all images that will be processed
# https://portal.opentopography.org/raster?opentopoID=OTSRTM.082016.4326.1
dem_file_path = '/path/to/dem.tif'


# ------------------- Start processing
def run_automated_image_preprocessing():
    for input_folder in (input_folders_array):
        root_paths = find_roots(input_folder)
        for root_folder_path, root_file_path in root_paths:
            found_default_files = find_files(root_folder_path, root_file_path, filter_basenames)
            for count, (photo_basename, found_default_file) in enumerate(found_default_files.items(), start=0):
                print(f"Processing {found_default_file['mul_photo_basename']}")





                print('----------Fetching variables from files') # -------------------- Define variables per image
                root_folder_path = found_default_file['root_folder_path']
                root_file_path = found_default_file['root_file_path']
                # photo_basename = found_default_file['photo_basename']

                mul_photo_basename = found_default_file['mul_photo_basename']
                mul_scene_basename = found_default_file['mul_scene_basename']
                mul_til_file = found_default_file['mul_til_file']
                mul_tif_file = found_default_file['mul_tif_file']
                mul_imd_file = found_default_file['mul_imd_file']

                pan_scene_basename = found_default_file['pan_scene_basename']
                pan_til_file = found_default_file['pan_til_file']
                pan_tif_file = found_default_file['pan_tif_file']
                pan_imd_file = found_default_file['pan_imd_file']
                pan_photo_basename = found_default_file['pan_photo_basename']

                params_overrides_scene, mul_params_overrides_photo, mul_imd_data = get_metadata_from_files(root_file_path, mul_imd_file, mul_photo_basename)
                _, pan_params_overrides_photo, pan_imd_data = get_metadata_from_files(root_file_path, pan_imd_file, pan_photo_basename)

                mul_shp_path = found_default_file['mul_shp_file']
                pan_shp_path = found_default_file['pan_shp_file']





                print('----------Setting custom variables') # -------------------- Set any custom resources per image (array must match legth of found_default_files)
                # if count == 0:






                print('----------Starting init') # -------------------- Initialize files
                # These shp files from maxar dont have a projection, this adds the correct projectionfor shp_path in [found_default_file['mul_shp_file'], found_default_file['pan_shp_file']]:
                mul_gpkg_path = os.path.join(input_folder, 'Mul_Footprint', f"{mul_photo_basename}.gpkg")
                shp_gpkg_path = os.path.join(input_folder, 'Pan_Footprint', f"{pan_photo_basename}.gpkg")
                mul_original_path = os.path.join(input_folder, 'Mul_Original', f"{mul_photo_basename}.tif")
                pan_original_path = os.path.join(input_folder, 'Pan_Original', f"{pan_photo_basename}.tif")

                shp_to_gpkg(mul_shp_path, mul_gpkg_path, 4326)
                shp_to_gpkg(pan_shp_path, shp_gpkg_path, 4326)

                # copy_tif_file(mul_tif_file, mul_original_path, nodata_value, dtype)
                # copy_tif_file(pan_tif_file, pan_original_path, nodata_value, dtype)

                # translate_gcp_image_to_origin(mul_til_file, mul_origin_path)
                # translate_gcp_image_to_origin(pan_til_file, pan_origin_path)





                print('----------Starting flaash') # -------------------- FLAASH Multispectral image
                mul_flaash_image_path = os.path.join(root_folder_path,'Mul_FLAASH', f"{mul_photo_basename}_FLAASH.dat")
                mul_flaash_params_path = os.path.join(root_folder_path,'Mul_FLAASH', f"{mul_photo_basename}_FLAASH_Params.txt")

                # Create flaash params
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
                # https://github.com/dawhite/MCTK (MODIS Conversion Toolkit (MCTK))
                # https://www.wunderground.com/history

                flaash_params = {
                    'INPUT_RASTER': {
                        'url': (mul_tif_file),
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
                    'GROUND_ELEVATION': get_image_largest_value(dem_file_path, mul_gpkg_path)/1000,
                    'SOLAR_AZIMUTH': mul_imd_data.get("SOLAR_AZIMUTH"),
                    'SOLAR_ZENITH': mul_imd_data.get("SOLAR_ZENITH"),
                    'LOS_AZIMUTH': mul_imd_data.get("LOS_AZIMUTH"),
                    'LOS_ZENITH': mul_imd_data.get("LOS_ZENITH"),
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
                    'OUTPUT_RASTER_URI': (mul_flaash_image_path),
                    'CLOUD_RASTER': None,
                    'WATER_RASTER': None,
                    'OUTPUT_RASTER': None,
                }

                if params_overrides_scene is not None: flaash_params.update(params_overrides_scene)
                if mul_params_overrides_photo is not None: flaash_params.update(mul_params_overrides_photo)
                flaash_params = {k: v for k, v in flaash_params.items() if v is not None}

                run_flaash(flaash_params, mul_flaash_params_path, envi_engine, mul_flaash_image_path) # -------------------- RUN






                print('----------Starting create gcp') # Step 2 -------------------- Create GCPs for image
                gcp_folder_name = "OrthoFromUSGSLidar" #'OrthoFromUSGSLidar', 'OrthoFromDefaultRPC'
                if not os.path.exists(os.path.join(root_folder_path, gcp_folder_name)): os.makedirs(os.path.join(root_folder_path, gcp_folder_name))
                mul_origin_path = os.path.join(root_folder_path,f'Mul_Origin', f"{mul_photo_basename}_Origin.tif")
                pan_origin_path = os.path.join(root_folder_path,f'Pan_Origin', f"{pan_photo_basename}_Origin.tif")
                qgis_gcp_file_path = os.path.join(root_folder_path, gcp_folder_name, f'{mul_photo_basename}.TIF.points')
                mul_gcp_geojson_file_path = os.path.join(root_folder_path, gcp_folder_name, f'{mul_photo_basename}_points.geojson')
                pan_gcp_geojson_file_path = os.path.join(root_folder_path, gcp_folder_name, f'{pan_photo_basename}_points.geojson')

                # *Manual create GCPs in QGIS*; if a reference image is available they can be automatically created; Create GCPs in QGIS and them in this folder with the name f'{mul_origin_path}.points';

                # Converts QGIS '.points' file to Orthority GCP geojon format for refinement for mul and pan
                qgis_gcps_to_geojson(mul_tif_file, qgis_gcp_file_path, os.path.basename(mul_flaash_image_path), dem_file_path, mul_gcp_geojson_file_path, force_positive_pixel_values=True)
                qgis_gcps_to_geojson(pan_tif_file, qgis_gcp_file_path, os.path.basename(pan_tif_file), dem_file_path, pan_gcp_geojson_file_path, force_positive_pixel_values=True)






                print('----------Starting orthorectify multispectral image') # Step 3 -------------------- Orthorectify Mulispectral Images
                # mul_ortho_default_rpc_model_path = os.path.join(root_folder_path, f'Mul_OrthoFromDefaultRPC', f"{mul_photo_basename}_OrthoFromDefaultRPC.tif")
                mul_flaash_ortho_path = os.path.join(root_folder_path, f'Mul_FLAASH_{gcp_folder_name}', f"{mul_photo_basename}_FLAASH_{gcp_folder_name}.tif")
                envi_mul_flaash_ortho_path = os.path.join(root_folder_path, f'Mul_FLAASH_{gcp_folder_name}', f"ENVI_{mul_photo_basename}_FLAASH_{gcp_folder_name}.dat")

                # Helpful command to put input dem into WGS84 vertical datum: gdalwarp -s_srs "+proj=longlat +datum=WGS84 +no_defs +geoidgrids=/mnt/c/Users/admin/Downloads/usa_geoid2012b/usa_geoid2012/g2012a_hawaii.gtx" -t_srs "+proj=longlat +datum=WGS84 +no_def" /mnt/d/demwgs84.tif /mnt/d/demwgs84_VWGS84.tif
                # Orthorectification info: https://up42.com/blog/how-to-perform-orthorectification-a-practical-guide
                # Ensure that the dem is in elipsoidal height, may need to convert it with a datum. proj geiods are here: https://download.osgeo.org/proj/vdatum/
                # If GCPs are provided this function will use them for RPC refinement, if not, it will use the default RPC model
                gcp_refined_rpc_orthorectification(mul_flaash_image_path, mul_flaash_ortho_path, dem_file_path, epsg, output_nodata_value=nodata_value, dtype='int16', gcp_geojson_file_path=mul_gcp_geojson_file_path, output_resolution=pan_imd_data.get("product_res") )






                print('----------Starting orthorectify panchromatic image') # Step 4 -------------------- Orthorectify Panchromatic Images
                pan_ortho_path = os.path.join(root_folder_path, f'Pan_{gcp_folder_name}', f"{pan_photo_basename}_{gcp_folder_name}.tif")

                gcp_refined_rpc_orthorectification(pan_tif_file, pan_ortho_path, dem_file_path, epsg, output_nodata_value=nodata_value, dtype='int16', gcp_geojson_file_path=pan_gcp_geojson_file_path, output_resolution=pan_imd_data.get("product_res"))






                print('----------Starting pansharpen multispectral image') # Step 5 -------------------- Pansharp Mulispectral Images with Panchromatic Images
                # From: https://doi.org/10.5194/isprsarchives-XL-1-W1-239-2013
                mul_pansharp_path = os.path.join(root_folder_path, f'Mul_FLAASH_{gcp_folder_name}_Pansharp', f"{mul_photo_basename}_FLAASH_{gcp_folder_name}_Pansharps.tif")

                pansharpen_image(mul_flaash_ortho_path, pan_ortho_path, mul_pansharp_path, change_nodata_value=nodata_value)






                # Optionally create another ortho from the multispectral imagery that is not FLAASHed that may work better with cloud masking algorythms like omnicloudmask
                mul_ortho_path = os.path.join(root_folder_path, f'Mul_{gcp_folder_name}', f"{mul_photo_basename}_{gcp_folder_name}.tif")

                # gcp_refined_rpc_orthorectification(mul_tif_file, mul_ortho_path, dem_file_path, epsg, output_nodata_value=nodata_value, dtype='UInt16', gcp_geojson_file_path=mul_gcp_geojson_file_path, output_resolution=mul_imd_data.get("product_res"))






                print('Done with image')

    print('Done with main process imagery')


run_automated_image_preprocessing()
