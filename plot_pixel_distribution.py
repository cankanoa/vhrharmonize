import rasterio
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
matplotlib.use('Agg')

def plot_pixel_distribution(input_image_path, bands=None, bins=256, save_plot_path=None):
    """
Plots the distribution (histogram) of pixel values for specified bands
in a raster file as line plots. Optionally saves the plot to a file.

Parameters
----------
input_image_path : str
Path to the raster image (e.g., GeoTIFF).

bands : set or list of int, optional
A collection of band indices (1-based) to plot. If not provided or None,
the function will plot all available bands in the image.
Example: {5, 3, 2}

bins : int, optional
Number of bins to use for the histogram. Default is 256.

save_plot_path : str, optional
Path to save the plot as an image file. If provided, the plot will
be saved to this path and not displayed.

Returns
-------
None (displays or saves a matplotlib figure).
    """
    # Open the raster
    with rasterio.open(input_image_path) as src:

        # If no bands specified, default to *all* bands
        if bands is None:
            bands = range(1, src.count + 1)

        plt.figure(figsize=(8, 5))

        # Loop over each requested band
        for b in sorted(bands):
            if b < 1 or b > src.count:
                print(f"Warning: Band {b} is out of range (1-{src.count}). Skipping.")
                continue

            # Read just the specified band
            band_data = src.read(b)

            # Flatten the band to 1D; ignore nodata by default (if masked via rasterio)
            # If the raster has an internal nodata mask, it's automatically
            # applied when reading with rasterio
            valid_data = band_data.compressed() if np.ma.isMaskedArray(band_data) else band_data.ravel()

            # Compute histogram
            hist, bin_edges = np.histogram(valid_data, bins=bins)

            # Plot as a line plot
            # We use the midpoints of the bin edges to plot a continuous line
            bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
            plt.plot(bin_centers, hist, label=f'Band {b}')

        # Set plot title to the file name
        file_name = input_image_path.split('/')[-1]
        plt.title(f"{file_name}", fontsize=7)
        plt.xlabel("Pixel Value")
        plt.ylabel("Frequency")
        plt.legend()
        plt.grid(True, linestyle='--', linewidth=0.4)
        plt.tight_layout()

        if save_plot_path:
            # Save the plot to the specified file path
            plt.savefig(save_plot_path, dpi=300)
            plt.close()
        else:
            # Show the plot
            plt.show()

# ----------------------------
# Example function call:
# ----------------------------
# plot_pixel_distribution(
#     "/path/to/your/image.tif",
#     bands={5, 3, 2},       # or just None to plot all
#     save_plot_path="/path/to/save/plot.png"
# )
