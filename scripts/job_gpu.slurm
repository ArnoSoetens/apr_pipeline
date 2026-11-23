#!/bin/bash
# Set up and run APR windows on one GPU. Submit from the repository root:
#   sbatch scripts/job_gpu.slurm                         every point that needs setup, every unfinished window
#   sbatch scripts/job_gpu.slurm iso_00001 iso_00002     these data points
#   sbatch scripts/job_gpu.slurm iso_00001/a000-a005     these windows (main.py sets up iso_00001)
#   sbatch scripts/job_gpu.slurm --list todo.txt         entries from a file
# The arguments after the script name ("$@" below) go to main.py and run.py unchanged.
# pmemd comes from run.pmemd in data/inputs/config.toml: set pmemd = "cuda" there.
# The log goes to data/logs/gpu_<job id>.out (data/logs/ is created whenever a script reads the config;
# slurm itself does not create it, and without it the job fails without a log).
#SBATCH --job-name=apr_gpu
#SBATCH --partition=gpu-a100-small
#SBATCH --time=03:55:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --account=research-me-pe
#SBATCH --output=data/logs/gpu_%j.out
set -euo pipefail

# -------------------------------
# DelftBlue modules
# -------------------------------
module purge
module load 2026
module load gpu
module load compilers
module load .compiler-2024
module load gcc/12.3.0
module load openmpi/5.0.10
module load cuda/12.5
module load bzip2/1.0.8

# -------------------------------
# pAPRika conda environment and AmberTools26
# -------------------------------
source ~/opt/miniforge3/etc/profile.d/conda.sh
conda activate paprika
source ~/opt/amber26/install/amber.sh
# Make AmberTools26 Python packages (pytraj) visible
export PYTHONPATH="$HOME/opt/amber26/install/lib/python3.11/site-packages${PYTHONPATH:+:$PYTHONPATH}"

cd "$SLURM_SUBMIT_DIR"
echo "===== $(date '+%F %T') on $(hostname), GPU: ${CUDA_VISIBLE_DEVICES:-none} ====="
command -v python
command -v pmemd.cuda

echo "===== Setup (main.py) ====="
python scripts/main.py "$@"

echo "===== Run (run.py, one window at a time on the GPU) ====="
python scripts/run.py --cpus 1 "$@"
