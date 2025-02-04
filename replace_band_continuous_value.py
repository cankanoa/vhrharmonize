import rasterio
import numpy as np
from scipy import ndimage
import os
from tqdm import tqdm
def replace_band_continuous_values_in_largest_segment(
        input_image_path,
        output_image_path,
        search_value,
        replace_value,
        nodata_value=None,
        dtype=None,
        connectivity=8,
        minimum_contiguous_count_to_remove=20
):
    """
Replaces pixel values in contiguous regions of a multi-band raster
where ALL bands share the same search_value at a given pixel,
but only for segments with pixel counts >= minimum_contiguous_count_to_remove.

The connectivity (4- or 8-connected) defines how adjacent pixels are considered
part of the same region.

Parameters
----------
input_image_path : str
Path to the input raster.
output_image_path : str
Path to the output raster.
search_value : int or float
The pixel value to search for.
replace_value : int or float
The pixel value to replace with.
nodata_value : int or float, optional
The NoData value to set in the output raster.
Defaults to the input raster's NoData value if not provided.
dtype : str or rasterio dtype, optional
Data type of the output raster. Defaults to the data type of the input raster.
connectivity : int
Pixel connectivity to define contiguous regions.
- 4 for orthogonal neighbors (N, S, E, W).
- 8 for orthogonal + diagonal neighbors.
minimum_contiguous_count_to_remove : int
Minimum number of pixels in a segment for it to be removed (default is 1).
"""

    # ------------------
    # 1. READ RASTER
    # ------------------
    with rasterio.open(input_image_path) as src:
        data = src.read()  # Shape: (bands, rows, cols)
        profile = src.profile

        # Set nodata_value default from source if not given
        if nodata_value is None:
            nodata_value = src.nodata

        # If dtype is specified, cast data and update profile
        if dtype is not None:
            data = data.astype(dtype)
            profile.update(dtype=dtype)

    if not os.path.exists(output_image_path):
        os.makedirs(output_image_path)

    # ------------------
    # 2. CREATE MASK
    # ------------------
    # Create a 2D mask of all bands == search_value
    # shape = (rows, cols), True where all bands match search_value
    mask = np.all(data == search_value, axis=0)

    # If there are no pixels matching search_value at all, just write out unchanged data
    if not mask.any():
        profile.update({'nodata': nodata_value})
        with rasterio.open(output_image_path, 'w', **profile) as dst:
            dst.write(data)
        return

    # ------------------
    # 3. LABEL THE MASK (CONNECTED COMPONENTS)
    # ------------------
    # Define the connectivity structure
    # For 4-connected, cross shape; for 8-connected, a 3x3 block of True
    if connectivity == 4:
        structure = np.array([[0,1,0],
                              [1,1,1],
                              [0,1,0]], dtype=bool)
    else:
        structure = np.ones((3, 3), dtype=bool)

    labeled_mask, num_features = ndimage.label(mask, structure=structure)

    # If no features found, write out unchanged data
    if num_features == 0:
        profile.update({'nodata': nodata_value})
        with rasterio.open(output_image_path, 'w', **profile) as dst:
            dst.write(data)
        return

    # ------------------
    # 4. IDENTIFY SEGMENTS TO REMOVE
    # ------------------
    # label_sizes[i] = number of pixels in label i
    # label=0 is background, so we skip it
    label_sizes = np.bincount(labeled_mask.ravel())[1:]  # skip label=0
    # Find labels where the size is >= minimum_contiguous_count_to_remove
    labels_to_remove = np.where(label_sizes >= minimum_contiguous_count_to_remove)[0] + 1  # +1 to match label index

    # Create a boolean mask for all regions meeting the size criteria
    removal_mask = np.isin(labeled_mask, labels_to_remove)

    # ------------------
    # 5. REPLACE VALUES ONLY IN THOSE SEGMENTS
    # ------------------
    data[:, removal_mask] = replace_value

    # ------------------
    # 6. WRITE OUTPUT
    # ------------------
    profile.update({'nodata': nodata_value})

    with rasterio.open(output_image_path, 'w', **profile) as dst:
        dst.write(data)