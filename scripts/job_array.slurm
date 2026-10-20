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
# Many array tasks start at the same moment, and `conda activate` starts a python that reads
# a lot of files from the shared filesystem; some tasks failed there. So instead:
#   1. wait a random few seconds, so tasks do not all start together
#   2. put the environment's bin folder on PATH (conda itself is never started)
#   3. check that its python works; if not, wait and try again a few times
# Find ENV_DIR with `conda env list`.
ENV_DIR="$HOME/opt/miniforge3/envs/paprika"

sleep $(( RANDOM % 30 ))

if [ ! -d "$ENV_DIR" ]; then
    echo "ERROR: environment folder not found: $ENV_DIR" >&2
    exit 1
fi
export CONDA_PREFIX="$ENV_DIR"
export PATH="$ENV_DIR/bin:$PATH"

env_ok=0
for attempt in 1 2 3 4 5; do
    if python -c "import numpy, tomllib" 2>/dev/null; then
        env_ok=1
        break
    fi
    echo "python of the environment is not usable yet (attempt $attempt of 5), waiting" >&2
    sleep $(( RANDOM % 20 + 5 ))
done
if [ "$env_ok" -ne 1 ]; then
    echo "ERROR: the environment $ENV_DIR is not usable from $(hostname)" >&2
    exit 1
fi

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
