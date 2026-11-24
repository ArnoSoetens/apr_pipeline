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
# array_task.sh recognizes shared-filesystem errors by these paths
export ENV_DIR AMBERTOOLS_HOME

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

# -------------------------------
# Python environment check
# -------------------------------
# Under load the shared filesystem sometimes reports files of the environment as missing
# (another file each time; a retry works). Import what this stage's Python will import, retrying
# a few times, so most of those errors happen here instead of in the middle of the stage.
# array_task.sh also reruns a script that still runs into one.
case "$1" in
    setup)    CHECK_MODULES=(setup_runner structure placement parametrize apr_setup) ;;
    analysis) CHECK_MODULES=(apr_analysis) ;;
    implicit) CHECK_MODULES=(implicit) ;;
    collect)  CHECK_MODULES=(collect_results) ;;
    *)        CHECK_MODULES=() ;;          # sim runs only bash and pmemd; cleanup_logs is light
esac

wait_for_python() {
    local attempt error
    for attempt in 1 2 3 4 5; do
        # sys.path as the scripts set it up: apr_pipeline/ and scripts/ (the working directory)
        if error=$(python -c '
import importlib, sys
sys.path[:0] = [sys.argv[1] + "/../apr_pipeline", sys.argv[1]]
import numpy, tomllib
for name in sys.argv[2:]:
    importlib.import_module(name)
' "$PWD" "${CHECK_MODULES[@]}" 2>&1); then
            return 0
        fi
        echo "Python environment not ready (attempt $attempt of 5): $(tail -n 1 <<< "$error")" >&2
        sleep $(( RANDOM % 20 + 5 ))
    done

    echo "ERROR: Python environment is not usable on $(hostname); last error:" >&2
    echo "$error" >&2
    return 1
}
if [ "$1" != "sim" ]; then
    wait_for_python
fi

# submit_arrays.py sets the working directory to scripts/.
bash "$PWD/array_task.sh" "$@"