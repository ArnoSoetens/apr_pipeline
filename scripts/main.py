import os
import subprocess
import sys
from rdkit import Chem

script_dir = os.path.dirname(os.path.abspath(__file__))

# bootstrap: just enough to import apr_config (apr_pipeline/ sits next to scripts/)
sys.path.insert(0, os.path.abspath(os.path.join(script_dir, "..", "apr_pipeline")))

from apr_config import CFG

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




subprocess.run(["python3", "create_data_directories.py", CFG["directories"]["input_xml"], data_dir], cwd=pipeline_dir)

HOST_PDB = os.path.join(host_dir, f"{HOST}.pdb")
HOST_mol2 = os.path.join(host_dir, f"{HOST}.mol2")
HOST_frcmod = os.path.join(host_dir, f"{HOST}.frcmod")
host = Chem.MolFromPDBFile(HOST_PDB, removeHs=False)


for name in sorted(os.listdir(data_dir)):
    data_point_dir = os.path.join(data_dir, name)      # full path to this isomer's folder
    if not os.path.isdir(data_point_dir):              # skip data_input.xml and other files
        continue

    smi_file = os.path.join(data_point_dir, f"{name}.smi")
    if not os.path.isfile(smi_file):                   # not a data point folder
        continue
    with open(smi_file) as f:
        smiles = f.read().strip()
                  # the real SMILES

    try:
        mol = smiles_to_pdb(smiles, os.path.join(data_point_dir, f"{GUEST}.pdb"))
    except ValueError as e:
        print(name, "FAILED:", e)
        continue

    chain = find_carbon_backbone_chain(mol)
    G1, G2 = pick_backbone_alignment_atoms(chain)

    try:
        guest_aligned = place_guest_in_cavity(mol, G1, G2)[0]
        generate_pdb(host,guest_aligned,out_dir=data_point_dir)

        parametrize(data_point_dir)

        setup_datapoint(data_point_dir, host_dir, G1, G2)
    except Exception as e:
        print(name, "FAILED:", e)
        continue





    

    

  

    
