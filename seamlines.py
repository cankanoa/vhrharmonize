import os
import numpy as np
import cv2
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from shapely.geometry import LineString
from scipy.spatial import Voronoi
from osgeo import gdal


# ------------------------- Step 1: Generate Initial Mosaic -------------------------
def generate_initial_mosaic(image_paths):
    """
    Generates an initial mosaic seamline network using Voronoi diagrams.
    If fewer than 4 images are provided, fall back to a simple grid-based approach.
    """
    centers = []
    for path in image_paths:
        with rasterio.open(path) as img:
            h, w = img.shape
            centers.append((w // 2, h // 2))

    if len(centers) < 4:
        print(f"Warning: Only {len(centers)} images provided. Voronoi requires at least 4 points.")
        print("Using a simple midpoint connection strategy instead.")

        # Create a fallback seamline as a simple straight connection between midpoints
        return [((centers[i][0], centers[i][1]), (centers[i + 1][0], centers[i + 1][1]))
                for i in range(len(centers) - 1)]

    # Compute Voronoi diagram
    return Voronoi(centers)


# ------------------------- Step 2: Adjust Seamline Vertices -------------------------
def adjust_seamline_vertices(vor, images):
    """
Adjusts the position of seamline vertices to avoid obstacles.
    """
    adjusted_vertices = []
    for vertex in vor.vertices:
        if not is_on_obstacle(vertex, images):
            adjusted_vertices.append(vertex)
    return np.array(adjusted_vertices)


def is_on_obstacle(vertex, images):
    """
Checks if a vertex falls on an obstacle using image data.
    """
    for img in images:
        if img['mask'][int(vertex[1]), int(vertex[0])] == 1:  # Assuming a binary mask
            return True
    return False


# ------------------------- Step 3: Cost Function Computation -------------------------
def calculate_cost_function(overlap_region):
    """
Computes the cost function based on gray differences and gradient in the overlapping area.
    """
    gray_diff = np.abs(overlap_region[:, :, 0] - overlap_region[:, :, 1])
    gradient = np.gradient(overlap_region[:, :, 0]) + np.gradient(overlap_region[:, :, 1])
    cost = gray_diff + np.abs(gradient)
    return cost


# ------------------------- Step 4: Ant Colony Algorithm -------------------------
def ant_colony_algorithm(cost_map, iterations=100, num_ants=50):
    """
Uses the continuous space ant colony algorithm to find optimal seamlines.
    """
    rows, cols = cost_map.shape
    pheromones = np.ones_like(cost_map) * 0.1
    best_path = None
    best_cost = float('inf')

    for _ in range(iterations):
        paths, costs = [], []
        for _ in range(num_ants):
            path, cost = find_ant_path(cost_map, pheromones)
            paths.append(path)
            costs.append(cost)

        update_pheromones(pheromones, paths, costs)

        min_cost_idx = np.argmin(costs)
        if costs[min_cost_idx] < best_cost:
            best_cost = costs[min_cost_idx]
            best_path = paths[min_cost_idx]

    return best_path


def find_ant_path(cost_map, pheromones):
    """
Finds a path from start to end based on probability selection.
    """
    rows, cols = cost_map.shape
    start = (0, np.random.randint(0, cols))
    end = (rows - 1, np.random.randint(0, cols))

    path = [start]
    total_cost = 0
    current = start

    while current[0] < rows - 1:
        next_steps = [
            (current[0] + 1, current[1] - 1),
            (current[0] + 1, current[1]),
            (current[0] + 1, current[1] + 1)
        ]
        next_steps = [p for p in next_steps if 0 <= p[1] < cols]

        probabilities = np.array([pheromones[p] / (cost_map[p] + 1e-6) for p in next_steps])
        probabilities /= probabilities.sum()

        next_step = next_steps[np.random.choice(len(next_steps), p=probabilities)]
        path.append(next_step)
        total_cost += cost_map[next_step]
        current = next_step

    return path, total_cost


def update_pheromones(pheromones, paths, costs):
    """
Updates pheromone levels.
    """
    pheromones *= 0.8
    for path, cost in zip(paths, costs):
        pheromone_contribution = 1.0 / (cost + 1e-6)
        for p in path:
            pheromones[p] += pheromone_contribution


# ------------------------- Step 5: Save Seamlines to GeoPackage -------------------------
def save_seamlines_to_gpkg(seamlines, output_path, input_filenames):
    """
Saves seamlines as layers in a GeoPackage.
Ensures the number of seamlines and filenames match.
    """
    if len(seamlines) != len(input_filenames):
        print(f"Warning: Mismatch in number of seamlines ({len(seamlines)}) and input files ({len(input_filenames)}).")
        print("Adjusting filenames list to match the number of seamlines.")
        input_filenames = input_filenames[:len(seamlines)]  # Trim filenames list if too long

    gdf = gpd.GeoDataFrame({'filename': input_filenames, 'geometry': [LineString(line) for line in seamlines]}, crs="EPSG:4326")
    gdf.to_file(output_path, driver="GPKG")

    print(f"Seamlines saved to {output_path}")


# ------------------------- Step 6: Mask Input Images -------------------------
def mask_images(seamlines, image_paths, output_folder, basename):
    """
Masks images using seamlines and saves them.
    """
    os.makedirs(output_folder, exist_ok=True)

    for img_path, seamline in zip(image_paths, seamlines):
        with rasterio.open(img_path) as src:
            masked, _ = mask(src, [LineString(seamline)], crop=True)
            output_path = os.path.join(output_folder, os.path.basename(img_path).replace(".tif", f"_{basename}.tif"))

            with rasterio.open(output_path, 'w', **src.meta) as dest:
                dest.write(masked)

    return output_folder


# ------------------------- Step 7: Merge Masked Images -------------------------
def merge_masked_images(masked_folder, output_path):
    """
Merges all masked images into one.
    """
    tif_files = [os.path.join(masked_folder, f) for f in os.listdir(masked_folder) if f.endswith(".tif")]

    vrt = gdal.BuildVRT("/vsimem/temp.vrt", tif_files)
    gdal.Translate(output_path, vrt)
    vrt = None


# ------------------------- Main Function -------------------------
def process_images(
        input_images_paths,
        output_folder_path,
        output_seamline_basename,
        use_seamlines_to_create_masked_images=False,
        create_merged_image=False,
        iterations=100,
        num_ants=50
):
    """
Main function that executes all steps.
    """
    # Generate initial seamline network
    vor = generate_initial_mosaic(input_images_paths)

    # Check if Voronoi was successful or if it's using fallback
    if isinstance(vor, list):
        # Fallback strategy: Use predefined seamlines, no cost function required
        seamlines = vor
    else:
        # Compute cost functions only for valid Voronoi diagrams
        cost_maps = [calculate_cost_function(vor) for _ in input_images_paths]

        # Apply ant colony optimization
        seamlines = [ant_colony_algorithm(cost_map, iterations, num_ants) for cost_map in cost_maps]

    # Save seamlines to GeoPackage
    gpkg_path = os.path.join(output_folder_path, f"{output_seamline_basename}.gpkg")
    save_seamlines_to_gpkg(seamlines, gpkg_path, [os.path.basename(f) for f in input_images_paths])

    # Optional: Mask images using seamlines
    if use_seamlines_to_create_masked_images:
        masked_folder = os.path.join(output_folder_path, "masked_images")
        mask_images(seamlines, input_images_paths, masked_folder, output_seamline_basename)

        # Optional: Merge masked images
        if create_merged_image:
            merged_path = os.path.join(output_folder_path, f"merged_{output_seamline_basename}.tif")
            merge_masked_images(masked_folder, merged_path)

    return gpkg_path


# Example Usage
if __name__ == "__main__":
    process_images(
        ["/Users/kanoalindiwe/Downloads/temp/testDroneImages/GeoNodata/DJI2_Geo_No.tif",
         "/Users/kanoalindiwe/Downloads/temp/testDroneImages/GeoNodata/DJI1_Geo_No.tif"
         ],
        "/Users/kanoalindiwe/Downloads/temp/testDroneImages/seamlines",
        "seamlines",
        use_seamlines_to_create_masked_images=True,
        create_merged_image=True
    )