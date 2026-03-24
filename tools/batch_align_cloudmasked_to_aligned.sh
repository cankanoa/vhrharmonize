#!/usr/bin/env bash
set -euo pipefail

MOVING_DIR="${MOVING_DIR:-/mnt/s/Satellite_Imagery/Big_Island/Processed_cloudmasked}"
FIXED_IMAGE="${FIXED_IMAGE:-/mnt/x/PROJECTS_2/Big_Island/ChangeHI_Trees/Dry_Forest/Data/Raster/mean_intensity/mean_intensity_mosaic.tif}"
TEMP_DIR="${TEMP_DIR:-/home/manumea/tmp_align}"
MOVING_BAND_INDEX="${MOVING_BAND_INDEX:-6}"
FIXED_BAND_INDEX="${FIXED_BAND_INDEX:-0}"

mkdir -p "${TEMP_DIR}"
shopt -s nullglob

inputs=("${MOVING_DIR}"/*_cloud_masked.tif)

if [[ ${#inputs[@]} -eq 0 ]]; then
    echo "No *_cloud_masked.tif files found in ${MOVING_DIR}" >&2
    exit 1
fi

for moving_image in "${inputs[@]}"; do
    output_image="${moving_image%_cloud_masked.tif}_aligned.tif"

    if [[ -f "${output_image}" ]]; then
        echo "Skipping: ${moving_image}"
        echo "Reason:   aligned output already exists at ${output_image}"
        continue
    fi

    echo "Aligning: ${moving_image}"
    echo "Output:   ${output_image}"

    vhr-align-image-pair \
      --moving-image "${moving_image}" \
      --fixed-image "${FIXED_IMAGE}" \
      --output-image "${output_image}" \
      --registration-mode structural_wv3_lidar \
      --moving-band-index "${MOVING_BAND_INDEX}" \
      --fixed-band-index "${FIXED_BAND_INDEX}" \
      --clip-fixed-to-moving \
      --enforce-mutual-valid-mask \
      --output-on-moving-grid \
      --temp-dir "${TEMP_DIR}"
done
