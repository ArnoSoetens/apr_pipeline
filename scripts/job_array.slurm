#!/bin/bash
# Shared by all array jobs. submit_arrays.py passes Slurm resources and stage arguments.
#   job_array.slurm setup|sim|analysis|implicit <worklist> <items_per_task> <offset> [...]
#   job_array.slurm collect
set -euo pipefail

ENV_DIR="$HOME/opt/miniforge3/envs/paprika"
AMBERTOOLS_HOME="$HOME/opt/amber26/install"

# Spread array-task startup to reduce simultaneous shared-filesystem load.
sleep $(( RANDOM % 30 ))

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
# pAPRika environment (PATH approach)
# -------------------------------
if [ ! -d "$ENV_DIR" ]; then
    echo "ERROR: environment folder not found: $ENV_DIR" >&2
    exit 1
fi
export PATH="$ENV_DIR/bin:$PATH"

# Retry the Python environment check in case the shared filesystem is slow.
wait_for_python() {
    local attempt
    for attempt in 1 2 3 4 5; do
        if python -c "import numpy, tomllib" 2>/dev/null; then
            return 0
        fi
        echo "Python environment not ready (attempt $attempt of 5); retrying" >&2
        sleep $(( RANDOM % 20 + 5 ))
    done

    echo "ERROR: Python environment is not usable on $(hostname)" >&2
    return 1
}
# the simulation tasks run only bash and pmemd (array_task.sh sim): no Python to wait for
if [ "$1" != "sim" ]; then
    wait_for_python
fi

# -------------------------------
# AmberTools26
# -------------------------------
source "$AMBERTOOLS_HOME/amber.sh"

# The antechamber wrapper needs AMBERHOME to locate wrapped_progs.
# The successful conda-activation test showed this should be the paprika env.
export AMBERHOME="$ENV_DIR"
export PATH="$AMBERTOOLS_HOME/bin:$ENV_DIR/bin:$PATH"
export PYTHONPATH="$AMBERTOOLS_HOME/lib/python3.11/site-packages${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1

if [ ! -x "$AMBERHOME/bin/wrapped_progs/antechamber" ]; then
    echo "ERROR: expected antechamber program not found at:" >&2
    echo "  $AMBERHOME/bin/wrapped_progs/antechamber" >&2
    exit 1
fi

echo "===== $(date '+%F %T') on $(hostname), stage: $1 ====="
echo "AMBERHOME=$AMBERHOME"
if [ "$1" != "sim" ]; then
    command -v python
    command -v antechamber
    command -v sqm
fi
command -v pmemd || true

# submit_arrays.py sets the working directory to scripts/.
bash "$PWD/array_task.sh" "$@"