import shutil

import orthority as oty
from osgeo import gdal

# image_path = '/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445319010/016445319010_01_003/016445319010_01_003/016445319010_01/016445319010_01_P005_MUL/17DEC08211801-M1BS-016445319010_01_P005.TIF'
image_path = '/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445319010/016445319010_01_003/016445319010_01_003/016445319010_01/016445319010_01_P003_MUL/17DEC08211758-M1BS-016445319010_01_P003.TIF'
# image_path = '/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445319010/Mul_Origin/17DEC08211758-M1BS-016445319010_01_P003_Origin.tif'

# image_path = "/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445319010/Mul_Origin/17DEC08211758-M1BS-016445319010_01_P003_Origin.tif"
# image_path = "/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445319010/Mul_Origin/17DEC08211758-M1BS-016445319010_01_P003_Origin.tif"
# image_path = '/mnt/d/temp/Image3Strip.tif'

# dem_path = "/mnt/x/Imagery/Lidar/Big_Island/2018_PuuWaawaa/DEM/2018_2020_bigIsland_DEM_J970216_000_000.tif"
# dem_path = '/mnt/d/srtm_05_09/srtm_05_09.tif'
# dem_path = '/mnt/d/demPA11_VWGS84.tif'
dem_path = '/mnt/d/dem_WGS84_Elipsoid.tif'
# dem_path = '/mnt/d/demlast.tif'
# dem_path = '/mnt/x/Imagery/Elevation/rasters_SRTMGL1Ellip/output_SRTMGL1Ellip.tif'


output_path = "/mnt/d/temp/15.tif"

gcp_path = '/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445319010/OrthoFromUSGSLidar/17DEC08211758-M1BS-016445319010_01_P003_points.geojson'
# gcp_path = "/mnt/s/Satellite_Imagery/Big_Island/Unprocessed/PuuWaawaaImages/20171208_36cm_WV03_BAB_016445319010/OrthoFromUSGSLidar/17DEC08211758-M1BS-016445319010_01_P003_Origin.tif.points"

# gcps_csv = '''mapX,mapY,sourceX,sourceY,z
# -155.88268224285898,19.852630821961572,759.0807081888352,-3383.7924066441033,0
# -155.8847164117608,19.851857923380926,612.3164895994481,-3445.6344544785657,0
# -155.88424580962285,19.851408016744028,645.3431118132443,-3484.334241067552,0
# -155.89264667817832,19.871142396195562,74.5582228689888,-1822.8339914874734,0
# -155.89192524057069,19.87141193080013,126.71813766975728,-1802.1212716628,0
# -155.86371306049654,19.796590964861206,2002.0527569211017,-8098.154190448766,0
# -155.8654402242259,19.797426179745557,1881.4800684838146,-8027.235129235269,0
# -155.8665691543928,19.798575185131803,1802.829581643047,-7926.780057108802,0
# -155.86645190774502,19.800174314781415,1814.2696283835649,-7794.448423772881,0
# -155.8859297425226,19.810347060368198,443.0224527009493,-6914.610454955455,0
# -155.8885322863178,19.80921409036554,254.67225788913666,-7005.838229995421,0
# -155.88587758656365,19.79265976169955,407.52897912771243,-8394.139017651907,0
# -155.86995921376038,19.80752885915965,1576.7294252884017,-7174.977258973007,0
# -155.87066073632187,19.80762122680466,1528.4906442201388,-7165.351665131555,0
# -155.8713687811009,19.81626478524517,1493.6312211524116,-6441.811981418693,0
# -155.87867207055103,19.82725761938007,997.1167965905919,-5509.316355804454,0
# -155.87760815116212,19.82660945223859,1072.04931326194,-5568.074141271933,0
# -155.88419124952588,19.839332509462896,626.6967672890803,-4493.658962612439,0
# -155.88178495179787,19.838693180563062,795.9805084530228,-4550.919266965376,0
# -155.88254614585824,19.846992300765777,760.8124952687676,-3855.4908043717296,0
# -155.88031049433772,19.84774197940573,921.7038341908274,-3795.905795373736,0
# -155.8422153186368,19.807396755820804,3551.3998175384004,-7227.862812543452,0
# -155.84288649673385,19.805629349726082,3502.2321898453306,-7374.935000385826,0
# -155.8426177286293,19.805537999118677,3522.8790825740693,-7383.353525099627,0
# -155.84328669370052,19.806009120589852,3474.3957541128357,-7342.981475924685,0
# -155.84544256384643,19.80768250690549,3322.1638067824315,-7199.707210491068,0
# -155.85687852817674,19.813682938226744,2517.6384450996675,-6680.7512243767505,0
# -155.85372145414027,19.811542869792305,2735.9233099942667,-6865.8114437705835,0
# -155.8464080723875,19.810780073975703,3253.5842853289846,-6939.5549460512675,0
# -155.84723721423035,19.811355841342667,3194.1710715838663,-6890.783258006431,0
# -155.88729767031674,19.823244717996836,375.3103845110496,-5832.031357817753,0
# -155.88571268024464,19.825144635045795,479.1283692998263,-5681.258494456099,0
# -155.88609877947687,19.832417550006692,476.92595330447307,-5068.070920645699,0
# -155.8737406557961,19.800700922183893,1293.551369081056,-7740.211083066013,0
# -155.83229941936023,19.793690458806562,4239.549295453114,-8389.847891473153,0
# -155.83337857463846,19.79424260713023,4163.971694776855,-8341.033596520669,0
# -155.83701509645456,19.79424596113735,3902.0008395674768,-8335.43891219004,0
# -155.8412164572948,19.79244705121494,2872.5154863042053,-8384.13016845262,0
# -155.85142157620797,19.799080957468234,2884.408688876767,-7908.839240669221,0
# -155.85685104976582,19.802221813055624,2502.3477394351194,-7638.385766881705,0
# -155.85651312282127,19.801050240660636,2523.3042718927363,-7737.576870329943,0
# -155.85532516050347,19.80944220156393,2618.564360829268,-7037.4628254403,0
# -155.85882624259142,19.827015867560014,2421.11330687231,-5561.502399028914,0
# -155.8596794062029,19.831197799421535,2365.983738646612,-5208.702917958509,0
# -155.86150484690933,19.83516372225365,2239.797469016504,-4877.012995534097,0
# -155.86764396614223,19.84105906343972,1812.9187689179507,-4373.438469924695,0
# -155.87153692266006,19.83353624173874,1519.0237472329957,-4996.815704639361,0
# -155.8772697981889,19.83891651521337,1122.3729557316863,-4538.971240539065,0
# -155.8899547201828,19.846070040440342,229.62885274962662,-3922.1689255349875,0
# -155.89260739400567,19.855646924103297,52.09990800695956,-3117.765787464898,0
# -155.884059063127,19.85963366883499,648.6974263707488,-2803.036670526702,0
# -155.88909435831272,19.862208099294268,313.5290097480124,-2574.0674986417007,0
# -155.88855269065587,19.800819485463474,235.0723019302906,-7707.572038052415,0
# -155.87337245653197,19.853641675437768,1431.1293360581562,-3312.9627144414185,0'''

# gdalwarp -s_srs "+proj=utm +zone=5 +datum=NAD83 +ellps=GRS80 +no_defs +geoidgrids=/mnt/c/Users/admin/Downloads/usa_geoid2012b/usa_geoid2012/g2012a_hawaii.gtx"          -t_srs "+proj=longlat +datum=WGS84 +no_defs"          /mnt/x/Imagery/Lidar/Big_Island/2018_PuuWaawaa/DEM/2018_2020_bigIsland_DEM_J970216_000_000.tif          /mnt/d/demlast.tif

def orthorectify_with_gcps(input_image_path, output_image_path, gcp_geojson_file_path, dem_image_path, output_epsg):
    # Create cameras from image RPC tags
    cameras = oty.RpcCameras.from_images([input_image_path])

    # Refine the RPC model with the provided GCPs
    cameras.refine(gcp_geojson_file_path)

    # Get the refined camera model
    camera = cameras.get(input_image_path)

    # Perform RPC rectification and save the output
    ortho = oty.Ortho(
        src_file=input_image_path,
        dem_file=dem_image_path,
        camera=camera,
        crs=f"EPSG:{output_epsg}"
    )

    ortho.process(
        ortho_file=output_image_path,
        resolution=None,
        interp="cubic",
        dem_interp="cubic",
        per_band=False,
        overwrite=True,
        progress=True
    )

    print(f"Orthorectified image saved to {output_image_path}")

orthorectify_with_gcps(image_path, output_path, gcp_path, dem_path, 4326)

# from osgeo import gdal
# import os
# from tempfile import NamedTemporaryFile
# def rpc_orthorectification(input_image_path, output_image_path, gcp_geojson_file_path, dem_image_path, output_epsg):
#     import orthority as oty
#
#     output_image_dir = os.path.dirname(output_image_path)
#     if not os.path.exists(output_image_dir):
#         os.makedirs(output_image_dir)
#
#     cameras = oty.RpcCameras.from_images([input_image_path])
#     cameras.refine(gcp_geojson_file_path)
#     camera = cameras.get(input_image_path)
#     camera = camera._rpc
#
#     with NamedTemporaryFile(delete=False, suffix=".tif") as temp_file:
#         temp_image_path = temp_file.name
#
#     # Copy the input image to the temporary file and explicitly remove all metadata
#     translate_options = gdal.TranslateOptions(
#         format="GTiff",
#         metadataOptions=["REMOVE_ALL"],  # Remove all metadata
#     )
#     gdal.Translate(temp_image_path, input_image_path, options=translate_options)
#
#     temp_dataset = gdal.Open(temp_image_path, gdal.GA_Update)
#
#     temp_dataset.SetMetadata(None)  # Remove all metadata
#     temp_dataset.SetMetadata(None, "RPC")  # Remove any existing RPC metadata
#
#     # final_output_path = '/mnt/d/temp/tempfile.tif'
#     # translate_options_no_metadata = gdal.TranslateOptions(
#     #     format="GTiff",
#     #     metadataOptions=["REMOVE_ALL"],  # Ensure all metadata is stripped
#     #     creationOptions=["COMPRESS=LZW"]  # Optional: Compress for efficiency
#     # )
#     # gdal.Translate(final_output_path, temp_image_path, options=translate_options_no_metadata)
#
#     # shutil.copy(temp_image_path, '/mnt/d/temp/tempfile.tif')
#
#     # Add the refined RPC metadata
#     refined_rpc_metadata = {
#         "LINE_OFF": str(camera.line_off),
#         "LINE_SCALE": str(camera.line_scale),
#         "SAMP_OFF": str(camera.samp_off),
#         "SAMP_SCALE": str(camera.samp_scale),
#         "LAT_OFF": str(camera.lat_off),
#         "LAT_SCALE": str(camera.lat_scale),
#         "LONG_OFF": str(camera.long_off),
#         "LONG_SCALE": str(camera.long_scale),
#         "HEIGHT_OFF": str(camera.height_off),
#         "HEIGHT_SCALE": str(camera.height_scale),
#         "LINE_NUM_COEFF": " ".join(map(str, camera.line_num_coeff)),
#         "LINE_DEN_COEFF": " ".join(map(str, camera.line_den_coeff)),
#         "SAMP_NUM_COEFF": " ".join(map(str, camera.samp_num_coeff)),
#         "SAMP_DEN_COEFF": " ".join(map(str, camera.samp_den_coeff)),
#     }
#     temp_dataset.SetMetadata(refined_rpc_metadata, "RPC")
#
#
#     temp_dataset = None  # Save and close the dataset
#
#     # Step 4: Perform orthorectification using GDAL
#     warp_options = gdal.WarpOptions(
#         format="GTiff",
#         dstSRS=f"EPSG:{output_epsg}",
#         rpc=True,
#         transformerOptions=[f"RPC_DEM={dem_image_path}"],
#         resampleAlg="cubic",
#         copyMetadata=True,  # Copy metadata excluding original RPC
#     )
#
#     # Perform the warp
#     gdal.Warp(
#         destNameOrDestDS=output_image_path,
#         srcDSOrSrcDSTab=temp_image_path,
#         options=warp_options,
#     )
#
#     # Clean up the temporary file
#     os.remove(temp_image_path)
#
#     print(f"Orthorectified image saved to {output_image_path}")
#
# rpc_orthorectification(image_path, output_path, gcp_path, dem_path, 4326)