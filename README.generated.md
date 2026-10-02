# apr_pipeline

Automated Attach-Pull-Release (APR) pipeline for calculating binding free energies of host-guest complexes using `paprIka` and `AmberTools`.

## Folder Layout

```text
.
├── apr_pipeline/              # Core library for system preparation
│   ├── apr_analysis.py        # Analysis scripts (TI-block bootstrap)
│   ├── apr_setup.py           # Paprika window and restraint setup
│   ├── create_data_directories.py  # XML/SMILES directory generation
│   ├── enumerate              # Isomer enumeration tool (external)
│   ├── parametrize.py         # Antechamber/GAFF2 parametrization
│   ├── placement.py           # Guest alignment and host-cavity placement
│   └── structure.py           # SMILES to 3D PDB conversion (ETKDG + UFF)
├── scripts/                   # Execution and orchestration scripts
│   ├── main.py                # Setup driver (from XML to APR windows)
│   ├── analysis.py            # Post-simulation data gathering
│   ├── run.sh / run_gpu.sh    # Simulation execution scripts (CPU/GPU)
│   └── *.in                   # Amber MD input templates (min, prod)
├── .gitignore                 # Excludes data, logs, and SLURM scripts
└── results.csv                # Final calculated binding affinities
```

## Excluded Files

To keep the repository clean, the following are excluded via `.gitignore`:
- **`data/`**: Contains all generated structures, Amber topology files, and simulation trajectories.
- **`**/*.slurm`**: SLURM submission scripts and log files (`*.out`, `*.err`).
- **`__pycache__`**: Python compiled files.

## Getting Started

### 1. Prerequisites
Ensure you have the following installed and in your PATH:
- **AmberTools** (specifically `antechamber`, `parmchk2`, `tleap`, and `pmemd`)
- **Python 3.x** with the following packages:
  - `rdkit`
  - `paprika`
  - `numpy`, `scipy`, `networkx`, `matplotlib`

### 2. Setup Phase
Place your host files (`BCD.pdb`, `BCD.mol2`, `BCD.frcmod`) and your input XML (isomers/SMILES) in a directory named `data/`.

Run the setup script to generate the structures and simulation windows:
```bash
python scripts/main.py
```
This script will:
1. Parse the input XML and create directories for each isomer.
2. Generate 3D coordinates from SMILES.
3. Place guests inside the host cavity and check for clashes.
4. Parametrize guests with GAFF2.
5. Create APR windows (Attach and Pull phases) using `paprIka`.

### 3. Simulation Phase
Execute the MD simulations for all generated windows. You can use the provided bash scripts or wrap them in your own SLURM submission logic:
```bash
# For CPU execution
./scripts/run.sh

# For GPU execution
export PMEMD=pmemd.cuda
./scripts/run_gpu.sh
```

### 4. Analysis Phase
Once the production trajectories (`production.nc`) are generated for all windows, calculate the binding free energies:
```bash
python scripts/analysis.py
```
This will produce a `results.json` in each data directory and a global `results.csv` in the root folder containing the calculated $dG_{bind}$ and SEM.