import orthority as oty
import rasterio
import numpy as np
import os


def pansharpen_image(
    input_low_resolution_path,
    input_high_resolution_path,
    output_image_path,
    change_nodata_value = None,
    ):

    pan_sharp = oty.PanSharpen(input_high_resolution_path, input_low_resolution_path)
    pan_sharp.process(output_image_path, write_mask=False, overwrite=True)

    if change_nodata_value is not None:
        _change_nodata_value(output_image_path, change_nodata_value, -32768)


def _change_nodata_value(
    input_image_path,
    new_nodata_value,
    old_nodata
    ):

    """
    Replaces all pixels equal to `old_nodata` with `new_nodata_value`, and rewrites the file cleanly
    without any auxiliary metadata issues.

    :param input_image_path: Path to the original raster file.
    :param new_nodata_value: Value to use as the new NoData.
    :param old_nodata: The pixel value currently used as NoData.
    """
    # Load into memory
    with rasterio.open(input_image_path) as src:
        profile = src.profile.copy()
        data = src.read()
        dtype = src.dtypes[0]

        print(f"Replacing {old_nodata} → {new_nodata_value}")

        if np.issubdtype(data.dtype, np.floating):
            mask = np.isclose(data, old_nodata)
        else:
            mask = data == old_nodata

        if np.any(mask):
            data[mask] = new_nodata_value
        else:
            print("No pixels matched old nodata value")

    # Remove original file and .aux.xml if it exists
    os.remove(input_image_path)
    aux_path = input_image_path + ".aux.xml"
    if os.path.exists(aux_path):
        os.remove(aux_path)

    # Write clean raster
    profile.update(nodata=new_nodata_value)

    with rasterio.open(input_image_path, 'w', **profile) as dst:
        dst.write(data)