#!/bin/bash
# Shared by all array jobs. There are no #SBATCH lines on purpose: submit_arrays.py passes
# partition, account, time, cores, memory and the array range from config.toml.
#   job_array.slurm setup|sim|analysis <worklist> <items_per_task> [offset]
#   job_array.slurm collect
set -e

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
# pAPRika conda environment
# -------------------------------

unset CONDA_SHLVL
source ~/opt/miniforge3/etc/profile.d/conda.sh
conda activate paprika


# -------------------------------
# AmberTools26
# -------------------------------
source ~/opt/amber26/install/amber.sh
export PYTHONPATH="$HOME/opt/amber26/install/lib/python3.11/site-packages${PYTHONPATH:+:$PYTHONPATH}"

# one thread per process: the parallelism comes from running many processes
export OMP_NUM_THREADS=1

echo "===== $(date '+%F %T') on $(hostname), $1 ====="
which python pmemd

# submit_arrays.py starts the job in the scripts folder (--chdir), so $PWD is the scripts folder
bash "$PWD/array_task.sh" "$@"
