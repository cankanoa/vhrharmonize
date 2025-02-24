import re
import os
from osgeo import gdal, gdal_array
import numpy as np
import os
import numpy as np
import rasterio
from rasterio import open as rio_open
from scipy.ndimage import gaussian_filter
from osgeo import ogr, osr
from osgeo import gdal, ogr, osr
import numpy as np

def convert_dn_to_radiance(input_image_path, output_image_path, gain, offset, dtype=None):
    ds = gdal.Open(input_image_path, gdal.GA_ReadOnly)
    if not ds:
        raise RuntimeError(f"Failed to open {input_image_path}")

    os.makedirs(os.path.dirname(output_image_path), exist_ok=True)

    if dtype is None:
        dtype = ds.GetRasterBand(1).DataType  # Get the datatype from the first band
    else:
        gdal_dtype = gdal.GetDataTypeByName(dtype)
        if gdal_dtype is None:
            raise ValueError(f"Unsupported GDAL data type: {dtype}")
        dtype = gdal_dtype  # Convert dtype to GDAL data type

    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(output_image_path, ds.RasterXSize, ds.RasterYSize, ds.RasterCount, dtype)
    out_ds.SetMetadata(ds.GetMetadata())

    # Copy RPC metadata
    rpc_metadata = ds.GetMetadata("RPC")
    if rpc_metadata:
        out_ds.SetMetadata(rpc_metadata, "RPC")

    gcps = ds.GetGCPs()
    if gcps:
        out_ds.SetGCPs(gcps, ds.GetGCPProjection())

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


def spectral_scale_image(input_image_path, output_image_path, input_image_scale, output_image_scale):
    """
    Scales an image from one range to another while preserving NoData values,
    copying the default metadata, and also copying RPC metadata.
    """

    os.makedirs(os.path.dirname(output_image_path), exist_ok=True)

    # Open input image
    dataset = gdal.Open(input_image_path, gdal.GA_ReadOnly)
    if dataset is None:
        raise ValueError("Could not open input image")

    # Get image metadata
    metadata = dataset.GetMetadata()
    nodata_value = dataset.GetRasterBand(1).GetNoDataValue()

    # Read the raster bands
    bands = [dataset.GetRasterBand(i+1).ReadAsArray().astype(np.float32)
             for i in range(dataset.RasterCount)]

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
    out_dataset = driver.Create(output_image_path,
                                dataset.RasterXSize,
                                dataset.RasterYSize,
                                dataset.RasterCount,
                                gdal.GDT_Float32)

    # Copy geotransform, projection, and default metadata
    out_dataset.SetGeoTransform(dataset.GetGeoTransform())
    out_dataset.SetProjection(dataset.GetProjection())
    out_dataset.SetMetadata(metadata)

    # Copy RPC metadata if present
    rpc_metadata = dataset.GetMetadata("RPC")
    if rpc_metadata:
        out_dataset.SetMetadata(rpc_metadata, "RPC")

    # Write scaled bands
    for i, scaled_band in enumerate(scaled_bands):
        out_band = out_dataset.GetRasterBand(i+1)
        out_band.WriteArray(scaled_band)
        if nodata_value is not None:
            out_band.SetNoDataValue(float(nodata_value))

    # Close datasets
    out_dataset = None
    dataset = None

def stretch_spectral_values(
        input_image_paths_array,
        output_image_folder,
        output_basename,
        stretch_dictionary,
        smoothing=0,
        dtype_override=None,
        offset=0
):
    """
    For each image in 'input_image_paths_array':
    1. Print "Processing image: [path]" once per image.
    2. Retrieve metadata, including nodata and dtype.
    3. For each band:
    a) Convert band data to float for processing.
    b) Mask out nodata so it doesn't affect stretching.
    c) Print original band stats (mean, std, min, max).
    d) For each (key,value) in stretch_dictionary:
    - If ends with '%', interpret as input percentile (include negative).
    - If ends with '@', interpret as input percentile (exclude negative).
    - Otherwise, interpret as a float.
    e) Use np.interp to stretch valid pixels from input breakpoints
    to output breakpoints (WITHOUT sorting).
    f) Optionally apply smoothing if smoothing > 0.
    g) Add the fixed offset (offset param) to valid pixels.
    h) Print stretched band stats (mean, std, min, max).
    i) Restore nodata pixels.
    j) Cast the final stretched band back to original dtype or override dtype if provided.
    4. Write the adjusted data to a new image in 'output_image_folder'.
    5. Print "Saved to: [output_path]".

    NOTES:
    - By removing sorting, you must ensure the dictionary keys (after interpretation)
    are in ascending order for np.interp to work correctly.
    - The new '@' suffix means "percentile ignoring negatives."
    - 'dtype_override' allows forcing an output dtype (e.g. 'float32').
    - 'offset' is added to the valid pixels after interpolation (and smoothing).
    """

    def parse_breakpoint_value(band_arr, valid_mask, val):
        """
        Convert the key or value 'val' into a numeric.

        Cases:
        1. If val is int or float, return as-is.
        2. If val ends with '@', interpret as percentile ignoring negative pixels.
        3. If val ends with '%', interpret as percentile with all valid pixels (including negatives).
        4. Else, parse it as a literal float.
        """
        if isinstance(val, (int, float)):
            # Already numeric
            return val

        val_str = str(val).strip()

        # CASE A: Suffix "@" => percentile ignoring negative values
        if val_str.endswith('@'):
            # everything before '@' is the percentile
            percentile = float(val_str[:-1])  # remove '@'
            # copy valid data and mask out negatives
            data_no_neg = band_arr[valid_mask].copy()
            data_no_neg[data_no_neg < 0] = np.nan
            return np.nanpercentile(data_no_neg, percentile)

        # CASE B: Suffix "%" => percentile including negative values
        elif val_str.endswith('%'):
            percentile = float(val_str[:-1])  # remove '%'
            return np.nanpercentile(band_arr[valid_mask], percentile)

        # CASE C: Otherwise, parse as float
        else:
            return float(val_str)

    # Ensure output folder exists
    os.makedirs(output_image_folder, exist_ok=True)

    # Process each image
    for input_image_path in input_image_paths_array:
        print(f"Processing image: {input_image_path}")  # <--- PRINT ONCE PER IMAGE

        with rasterio.open(input_image_path) as src:
            profile = src.profile.copy()

            # Original dtype from the source
            original_dtype = profile['dtype']
            # Decide final output dtype
            final_dtype = dtype_override if dtype_override else original_dtype

            nodata_value = profile.get('nodata', None)

            out_bands = []

            for bidx in range(1, src.count + 1):
                band_data = src.read(bidx).astype(np.float32, copy=False)

                # Create mask for nodata
                if nodata_value is not None:
                    mask = (band_data == nodata_value)
                else:
                    # If no explicit nodata, you could treat NaNs as masked or none at all
                    mask = np.isnan(band_data)

                valid_mask = ~mask

                # Compute original stats (only for valid pixels)
                orig_valid = band_data[valid_mask]
                original_mean = float(np.nanmean(orig_valid)) if len(orig_valid) > 0 else np.nan
                original_std  = float(np.nanstd(orig_valid))  if len(orig_valid) > 0 else np.nan
                original_min  = float(np.nanmin(orig_valid))  if len(orig_valid) > 0 else np.nan
                original_max  = float(np.nanmax(orig_valid))  if len(orig_valid) > 0 else np.nan

                # Build input_vals/output_vals in dictionary order (NO sorting!)
                input_vals = []
                output_vals = []

                for k, v in stretch_dictionary.items():
                    in_val = parse_breakpoint_value(band_data, valid_mask, k)
                    out_val = parse_breakpoint_value(band_data, valid_mask, v)
                    input_vals.append(in_val)
                    output_vals.append(out_val)

                # Interpolate only valid pixels
                band_stretched = np.full_like(band_data, np.nan, dtype=np.float32)
                band_stretched[valid_mask] = np.interp(
                    band_data[valid_mask],
                    input_vals,
                    output_vals
                )

                # Optional smoothing
                if smoothing > 0:
                    temp_data = band_stretched.copy()
                    temp_data[mask] = 0
                    temp_data = gaussian_filter(temp_data, sigma=smoothing)
                    band_stretched[valid_mask] = temp_data[valid_mask]

                # Add offset to valid pixels
                if offset != 0:
                    band_stretched[valid_mask] += offset

                # Compute stretched stats (again only for valid pixels)
                stretched_valid = band_stretched[valid_mask]
                stretched_mean = float(np.nanmean(stretched_valid)) if len(stretched_valid) > 0 else np.nan
                stretched_std  = float(np.nanstd(stretched_valid))  if len(stretched_valid) > 0 else np.nan
                stretched_min  = float(np.nanmin(stretched_valid))  if len(stretched_valid) > 0 else np.nan
                stretched_max  = float(np.nanmax(stretched_valid))  if len(stretched_valid) > 0 else np.nan

                print(
                    f"Stats for band {bidx}: "
                    f"mean: {stretched_mean:.3f} vs {original_mean:.3f}, "
                    f"std: {stretched_std:.3f} vs {original_std:.3f}, "
                    f"min: {stretched_min:.3f} vs {original_min:.3f}, "
                    f"max: {stretched_max:.3f} vs {original_max:.3f}"
                )

                # Restore nodata
                if nodata_value is not None:
                    band_stretched[mask] = nodata_value

                # Cast back to final dtype
                band_stretched = band_stretched.astype(final_dtype)

                out_bands.append(band_stretched)

            # Update profile
            profile.update(dtype=final_dtype, count=src.count)

            # Output filename
            input_basename = os.path.basename(input_image_path)
            name_only, ext = os.path.splitext(input_basename)
            output_filename = f"{name_only}{output_basename}{ext}"
            output_path = os.path.join(output_image_folder, output_filename)

            # Write
            with rasterio.open(output_path, 'w', **profile) as dst:
                for i in range(src.count):
                    dst.write(out_bands[i], i+1)

        print(f"Saved to: {output_path}\n")

def get_image_largest_value(input_image_path, mask=None, override_mask_crs_epsg=None):
    """
    Get the largest value in a raster image, optionally masked by a shapefile.

    Args:
    input_image_path (str): Path to the raster image.
    mask (str, optional): Path to a geographic mask (shapefile). Defaults to None.
    override_mask_crs_epsg (int, optional): EPSG code to override the mask's coordinate system. Defaults to None.

    Returns:
    float: The largest value in the raster (masked if mask is provided).
    """
    # Open the raster dataset
    dataset = gdal.Open(input_image_path, gdal.GA_ReadOnly)
    if not dataset:
        raise FileNotFoundError(f"Raster file not found: {input_image_path}")

    # Get raster geotransform and nodata value
    geo_transform = dataset.GetGeoTransform()
    projection = dataset.GetProjection()
    x_size = dataset.RasterXSize
    y_size = dataset.RasterYSize

    band = dataset.GetRasterBand(1)
    nodata_value = band.GetNoDataValue()

    # Read the raster data as an array
    raster_array = band.ReadAsArray()
    if nodata_value is not None:
        raster_mask = (raster_array != nodata_value)
    else:
        raster_mask = np.ones_like(raster_array, dtype=bool)

    # If a mask is provided, check for overlaps in valid data
    mask_array = None
    if mask:
        mask_ds = ogr.Open(mask)
        if not mask_ds:
            raise FileNotFoundError(f"Mask file not found: {mask}")

        mask_layer = mask_ds.GetLayer()

        # Override the mask's CRS if requested
        if override_mask_crs_epsg:
            target_srs = osr.SpatialReference()
            target_srs.ImportFromEPSG(override_mask_crs_epsg)

            # Create a memory layer with the corrected CRS
            driver = ogr.GetDriverByName("Memory")
            mem_ds = driver.CreateDataSource("in_memory")
            mem_layer = mem_ds.CreateLayer("mask", srs=target_srs, geom_type=mask_layer.GetGeomType())

            # Copy features from the original mask layer
            for feature in mask_layer:
                mem_layer.CreateFeature(feature.Clone())

            mask_layer = mem_layer  # Replace original mask layer with corrected one

        # Create an in-memory raster for the mask
        driver = gdal.GetDriverByName("MEM")
        mask_raster = driver.Create("", x_size, y_size, 1, gdal.GDT_Byte)
        mask_raster.SetGeoTransform(geo_transform)
        mask_raster.SetProjection(projection)

        # Rasterize the shapefile to the mask raster
        gdal.RasterizeLayer(mask_raster, [1], mask_layer, burn_values=[1])
        mask_array = mask_raster.ReadAsArray()

        # Check if the mask has true values where the raster has no valid data
        overlap_check = np.logical_and(mask_array == 1, ~raster_mask)
        if np.any(overlap_check):
            print("Warning: The provided mask contains areas where there is no valid data in the input image.")

        # Clean up
        mask_ds = None
        mask_raster = None
        if override_mask_crs_epsg:
            mem_ds = None

    # Loop through all bands and find the largest value
    largest_value = None
    for band_index in range(1, dataset.RasterCount + 1):
        band = dataset.GetRasterBand(band_index)
        band_array = band.ReadAsArray()

        if nodata_value is not None:
            band_array = np.ma.masked_where(band_array == nodata_value, band_array)

        if mask_array is not None:
            # Apply the mask to the band data
            band_array = np.ma.masked_where(mask_array == 0, band_array)

        # Find the maximum value in the current band
        band_max = np.max(band_array)
        largest_value = max(largest_value, band_max) if largest_value is not None else band_max

    # Clean up
    dataset = None

    return float(largest_value)