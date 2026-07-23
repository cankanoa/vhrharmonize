#!/bin/bash
# Edit SBATCH lines for the cluster; vhr-hpc passes the remote provider config as $1


# Node settings
#SBATCH --account=uh
#SBATCH --partition=shared
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --job-name=vhr-worldview
#SBATCH --output=../logs/slurm-%j.out
#SBATCH --error=../logs/slurm-%j.err


# Setup conda
set -eo pipefail
module load lang/Anaconda3
eval "$(conda shell.bash hook)"


# Setup conda environment
ENV_NAME="vhrharmonize"
if [ ! -d "$(conda info --base)/envs/${ENV_NAME}" ]; then
    conda create -y -n "${ENV_NAME}" -c pytorch -c conda-forge py6s sixs gdal pytorch torchvision cpuonly imageio python=3.11
fi
conda activate "${ENV_NAME}"


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


# Run with configs from yml
export PYTHONUNBUFFERED=1
vhr-worldview --config-yaml "$1"


# Use distributed workers with dask across nodes
# #SBATCH --nodes=2
# #SBATCH --ntasks-per-node=1
# SCHEDULER_FILE="${HOME}/koa_scratch/vhrharmonize/runs/{run_id}/logs/dask-scheduler_{run_id}.json"
# SCHEDULER_HOST="$(scontrol show hostnames "${SLURM_JOB_NODELIST}" | head -n 1)"
# srun --nodes=1 --ntasks=1 -w "${SCHEDULER_HOST}" dask-scheduler --scheduler-file "${SCHEDULER_FILE}" &
# sleep 10
# srun --nodes="${SLURM_NNODES}" --ntasks="${SLURM_NNODES}" dask-worker --scheduler-file "${SCHEDULER_FILE}" --nthreads "${SLURM_CPUS_PER_TASK}" --local-directory "${SLURM_TMPDIR:-/tmp}" &
# sleep 10
# vhr-worldview --config-yaml "$1" --dask-scheduler-file "${SCHEDULER_FILE}"
