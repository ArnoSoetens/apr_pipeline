"""Set up data points: guest structure, placement in the host, parametrization, APR windows.

    python main.py                  serial: create data point folders from the configured SMILES list, then set up
                                    every data point whose windows are not complete
    python main.py --force          ... and redo the complete ones too
    python main.py iso_00001 ...    set up just these data points (this is what the slurm array tasks do)

A data point that fails gets a setup_failed.txt with the reason; one that succeeds loses it.
"""
import argparse
import os
import subprocess
import sys

script_dir = os.path.dirname(os.path.abspath(__file__))

# bootstrap: just enough to import apr_config (apr_pipeline/ sits next to scripts/)
sys.path.insert(0, os.path.abspath(os.path.join(script_dir, "..", "apr_pipeline")))

from apr_config import CFG

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("names", nargs="*", metavar="NAME",
                    help="data points to set up (default: all that still need it)")
parser.add_argument("--force", action="store_true",
                    help="serial mode: also redo data points whose windows are complete")
args = parser.parse_args()

from rdkit import Chem

data_dir = CFG["directories"]["data"]
host_dir = CFG["directories"]["host"]
pipeline_dir = CFG["directories"]["pipeline"]
HOST = CFG["naming"]["host_prefix"]
GUEST = CFG["naming"]["guest_prefix"]

sys.path.insert(0, pipeline_dir)   # lets Python find the "pipeline" package

from structure import *
from parametrize import *
from placement import *
from apr_setup import *
from datapoints import (
    SETUP_STATES,
    clear_failed,
    mark_failed,
    prepared_data_points,
    scan,
    setup_data_points,
)


serial = not args.names
if serial:
    smiles_list = CFG["directories"]["smiles_list"]
    if not prepared_data_points() and (
        os.path.isfile(smiles_list) or not setup_data_points()
    ):
        subprocess.run(
            [
                sys.executable,
                "create_data_directories.py",
                smiles_list,
                data_dir,
            ],
            cwd=pipeline_dir,
            check=True,
        )
    names = [dp.name for dp in setup_data_points()]
else:
    names = args.names

HOST_PDB = os.path.join(host_dir, f"{HOST}.pdb")
HOST_mol2 = os.path.join(host_dir, f"{HOST}.mol2")
HOST_frcmod = os.path.join(host_dir, f"{HOST}.frcmod")
host = Chem.MolFromPDBFile(HOST_PDB, removeHs=False)


n_done = n_failed = n_skipped = 0
for i, name in enumerate(names, 1):
    data_point_dir = os.path.join(data_dir, name)      # full path to this isomer's folder
    print(f"[{i}/{len(names)}] {name}", flush=True)

    smi_file = os.path.join(data_point_dir, f"{name}.smi")
    prepared_files = (
        os.path.join(data_point_dir, f"{GUEST}.pdb"),
        os.path.join(data_point_dir, f"{GUEST}.mol2"),
        os.path.join(data_point_dir, f"{GUEST}.frcmod"),
        os.path.join(data_point_dir, f'{CFG["naming"]["complex_prefix"]}.pdb'),
    )
    use_prepared_files = all(os.path.isfile(path) for path in prepared_files)
    has_smiles = os.path.isfile(smi_file)
    if not has_smiles and not use_prepared_files:
        print(
            name,
            f"SKIPPED: needs {name}.smi or all prepared guest files "
            f"({GUEST}.pdb, {GUEST}.mol2, {GUEST}.frcmod, "
            f'{CFG["naming"]["complex_prefix"]}.pdb)',
        )
        n_skipped += 1
        continue
    if serial and not args.force and scan(data_point_dir)[0] not in SETUP_STATES:
        print(name, "SKIPPED: windows are complete (--force redoes them)")
        n_skipped += 1
        continue

    try:
        if use_prepared_files:
            mol = Chem.MolFromPDBFile(prepared_files[0], removeHs=False)
            if mol is None:
                raise ValueError(f"could not read prepared guest PDB: {prepared_files[0]}")
        else:
            with open(smi_file) as f:
                smiles = f.read().strip()
            mol = smiles_to_pdb(smiles, prepared_files[0])

        chain = find_carbon_backbone_chain(mol)
        if len(chain) < 2:
            raise ValueError(f"guest needs at least two connected carbon atoms: {prepared_files[0]}")
        G1, G2 = pick_backbone_alignment_atoms(chain)

        if not use_prepared_files:
            guest_aligned = place_guest_in_cavity(mol, G1, G2)[0]
            generate_pdb(host, guest_aligned, out_dir=data_point_dir)

            parametrize(data_point_dir, net_charge=Chem.GetFormalCharge(mol))

        setup_datapoint(data_point_dir, host_dir, G1, G2)
    except Exception as e:
        print(name, "FAILED:", e)
        mark_failed(data_point_dir, e)
        n_failed += 1
        continue

    clear_failed(data_point_dir)
    n_done += 1

print(f"setup finished: {n_done} done, {n_failed} failed, {n_skipped} skipped", flush=True)
