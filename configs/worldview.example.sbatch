#!/bin/bash
# Edit SBATCH lines for the cluster; vhr-slurm passes the remote provider config as $1
#SBATCH --account=uh
#SBATCH --partition=shared
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=vhr-worldview
#SBATCH --output=../logs/slurm-%j.out
#SBATCH --error=../logs/slurm-%j.err

set -eo pipefail

# Setup conda
module load lang/Anaconda3
eval "$(conda shell.bash hook)"

# Setup conda environment
# conda create -n vhrharmonize -c conda-forge py6s sixs gdal python=3.11
conda activate vhrharmonize

# Default: install or upgrade the latest PyPI release
# python -m pip install --upgrade vhrharmonize

# Alternative: install from source instead of PyPI
# source_dir="${HOME}/vhrharmonize_source"
# if [ ! -d "${source_dir}/.git" ]; then
#   git clone https://github.com/your-org/vhrharmonize.git "${source_dir}"
# else
#   git -C "${source_dir}" pull --ff-only
# fi
# git -C "${source_dir}" checkout branch-name
# python -m pip install --upgrade "${source_dir}"[defaults]

vhr-worldview --config-yaml "$1"
