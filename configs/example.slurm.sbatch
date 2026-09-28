#!/bin/bash
# Edit SBATCH lines for the cluster; vhr hpc-start passes the remote workflow config as $1


# Node settings
#SBATCH --account=uh
#SBATCH --partition=shared
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --job-name=vhr-workflow
#SBATCH --output=../run_{run_id}/logs/slurm-%j.out
#SBATCH --error=../run_{run_id}/logs/slurm-%j.err

# When dask is enabled
# #SBATCH --nodes=2
# #SBATCH --ntasks-per-node=1

# Setup conda
set -eo pipefail
module load lang/Anaconda3
eval "$(conda shell.bash hook)"


# Setup conda environment
ENV_NAME="vhrharmonize"
if [ ! -d "$(conda info --base)/envs/${ENV_NAME}" ]; then
    conda create -y -n "${ENV_NAME}" -c pytorch -c conda-forge gdal sixs pytorch torchvision cpuonly imageio python=3.11
fi
conda activate "${ENV_NAME}"
export PYTHONUNBUFFERED=1


# Must install via Pypi or from git source:
# Install from Pypi
python -m pip install --upgrade "vhrharmonize[defaults,dask]"

# Alternative: install from source
# SOURCE_REPO_URL="https://github.com/cankanoa/vhrharmonize.git"
# SOURCE_DIR="${HOME}/vhrharmonize_source"
# SOURCE_BRANCH="main"
# SOURCE_EXTRAS="[defaults,dask]"
# if [ ! -d "${SOURCE_DIR}/.git" ]; then
#     git clone --branch "${SOURCE_BRANCH}" "${SOURCE_REPO_URL}" "${SOURCE_DIR}"
# else
#     git -C "${SOURCE_DIR}" fetch origin
#     git -C "${SOURCE_DIR}" checkout "${SOURCE_BRANCH}"
#     git -C "${SOURCE_DIR}" pull --ff-only origin "${SOURCE_BRANCH}"
# fi
# python -m pip install --upgrade "${SOURCE_DIR}${SOURCE_EXTRAS}"


# Run the prepared YAML (top-level plugin names; param:/var:/const:/core: settings).
# jsonata-python is installed by the package dependency above.
vhr workflow --config "$1"

# When Dask is enabled, set shared.core:concurrent_processing_backend: dask
# and shared.core:concurrent_processing: 1 in the recipe; comment the command above out.
# SpectralMatch aggregate functions inherit ordinary parameters from shared:
# param:concurrent_processing_backend: dask
# param:image_threads: null
# param:dask_scheduler: [file, /remote/path/to/the/scheduler.json]
# VHR seamline_metadata also needs param:dask_scheduler_file pointing to that file.
# Use the same scheduler file below; core settings are not passed as function parameters.
# SCHEDULER_FILE="${HOME}/koa_scratch/vhrharmonize/run_{run_id}/logs/{run_id}_dask_scheduler.json"
# SCHEDULER_HOST="$(scontrol show hostnames "${SLURM_JOB_NODELIST}" | head -n 1)"
# mkdir -p "$(dirname "${SCHEDULER_FILE}")"
# rm -f "${SCHEDULER_FILE}"
# srun --nodes=1 --ntasks=1 -w "${SCHEDULER_HOST}" dask-scheduler --scheduler-file "${SCHEDULER_FILE}" &
# sleep 10
# srun --nodes="${SLURM_NNODES}" --ntasks="${SLURM_NNODES}" dask-worker --scheduler-file "${SCHEDULER_FILE}" --nthreads "${SLURM_CPUS_PER_TASK}" --local-directory "${SLURM_TMPDIR:-/tmp}" &
# sleep 10
# vhr workflow --config "$1"
