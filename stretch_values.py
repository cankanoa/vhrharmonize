import os
import numpy as np
import rasterio
from rasterio import open as rio_open
from scipy.ndimage import gaussian_filter

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