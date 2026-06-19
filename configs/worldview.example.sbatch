#!/bin/bash
# Edit SBATCH lines for the cluster; vhr-slurm passes the remote provider config as $1.
#SBATCH --account=uh
#SBATCH --partition=shared
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=vhr-worldview
#SBATCH --output=../logs/slurm-%j.out
#SBATCH --error=../logs/slurm-%j.err

set -euo pipefail

source ~/.bashrc
conda activate vhrharmonize

vhr-worldview --config-yaml "$1"
