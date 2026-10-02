#!/bin/bash
#SBATCH --job-name="APR-test"
#SBATCH --partition=gpu-a100-small
#SBATCH --time=03:55:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --account=research-me-pe

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
source ~/opt/miniforge3/etc/profile.d/conda.sh
conda activate paprika

# -------------------------------
# AmberTools26
# -------------------------------
source ~/opt/amber26/install/amber.sh

# Make AmberTools26 Python packages (pytraj) visible
export PYTHONPATH="$HOME/opt/amber26/install/lib/python3.11/site-packages${PYTHONPATH:+:$PYTHONPATH}"

# -------------------------------
# Check environment
# -------------------------------
echo "===== Environment ====="
echo "AMBERHOME=$AMBERHOME"
which python
python --version

echo "pAPRika:"
python -c "import paprika; print(paprika.__file__)"
python -c "import importlib.metadata as m; print(m.version('paprika'))"

echo "Amber Python packages:"
python -c "import pytraj; print('pytraj:', pytraj.__file__)"
python -c "import parmed; print('ParmEd:', parmed.__file__)"

echo "pAPRika build:"
python -c "from paprika.build import align, dummy; print('paprika.build: OK')"

echo "Amber executables:"
which sander
which pmemd
which pmemd.MPI

# -------------------------------
# Run
# -------------------------------
echo "===== Setup (main.py) ====="
python "$SLURM_SUBMIT_DIR/main.py"

export PMEMD=pmemd.cuda
export MAX_JOBS=2

echo "===== Running run.sh ====="
bash "$SLURM_SUBMIT_DIR/run.sh"
