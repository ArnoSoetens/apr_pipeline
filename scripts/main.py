"""Set up data points: guest structure, placement in the host, parametrization, APR windows.

    python main.py                  serial: create the data point folders from the XML, then set up
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
from datapoints import SETUP_STATES, clear_failed, data_points, mark_failed, scan


serial = not args.names
if serial:
    # no names given: make the folders from the XML first, then take every data point
    subprocess.run(["python3", "create_data_directories.py", CFG["directories"]["input_xml"], data_dir], cwd=pipeline_dir)
    names = [dp.name for dp in data_points()]
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
    if not os.path.isfile(smi_file):                   # not a data point folder
        print(name, "SKIPPED: no .smi file")
        n_skipped += 1
        continue
    if serial and not args.force and scan(data_point_dir)[0] not in SETUP_STATES:
        print(name, "SKIPPED: windows are complete (--force redoes them)")
        n_skipped += 1
        continue
    with open(smi_file) as f:
        smiles = f.read().strip()
                  # the real SMILES

    try:
        mol = smiles_to_pdb(smiles, os.path.join(data_point_dir, f"{GUEST}.pdb"))
    except ValueError as e:
        print(name, "FAILED:", e)
        mark_failed(data_point_dir, e)
        n_failed += 1
        continue

    try:
        chain = find_carbon_backbone_chain(mol)
        G1, G2 = pick_backbone_alignment_atoms(chain)

        guest_aligned = place_guest_in_cavity(mol, G1, G2)[0]
        generate_pdb(host,guest_aligned,out_dir=data_point_dir)

        parametrize(data_point_dir)

        setup_datapoint(data_point_dir, host_dir, G1, G2)
    except Exception as e:
        print(name, "FAILED:", e)
        mark_failed(data_point_dir, e)
        n_failed += 1
        continue

    clear_failed(data_point_dir)
    n_done += 1

print(f"setup finished: {n_done} done, {n_failed} failed, {n_skipped} skipped", flush=True)
