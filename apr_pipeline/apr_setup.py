#!/usr/bin/env python
# Setup APR windows for BCD-FOA, orientation 1, implicit solvent (GB).
# Produces complex/ and windows/ which are renamed by run.sh.

import os
import subprocess

import numpy as np
import parmed as pmd

from paprika import restraints
from paprika.build import align, dummy
from paprika.build.system import TLeap
from paprika.io import save_restraints
from paprika.restraints.amber import amber_restraint_line
from paprika.restraints.utils import create_window_list

import shutil


##host = "BCD"
##guest = "UNL"


# --- Build vacuum complex with tleap ---

def build_system( data_dir_host, data_dir_guest, host_prefix = "BCD", guest_prefix = "UNL"):
    os.makedirs("complex", exist_ok=True)
    system = TLeap()
    system.output_path = "complex"
    system.pbc_type = None
    system.neutralize = False
    system.template_lines = [
        "source leaprc.gaff",
        "loadamberparams frcmod.ions1lm_126_tip3p",
        f"loadamberparams {data_dir_host}/{host_prefix}.frcmod",
        f"{host_prefix} = loadmol2 {data_dir_host}/{host_prefix}.mol2",
        f"loadamberparams {data_dir_guest}/{guest_prefix}.frcmod",
        f"{guest_prefix} = loadmol2 {data_dir_guest}/{guest_prefix}.mol2",
        f"model = loadpdb {data_dir_guest}/{host_prefix}-{guest_prefix}.pdb",
        "check model",
        "savepdb model vac.pdb",
        "saveamberparm model vac.prmtop vac.rst7",
    ]
    system.build(clean_files=False)
    return

# --- Align the pulling axis and add dummy atoms ---
# Orientation 1: C3 leads (faces away from host opening during pull).

def guest_masks(b1_idx, b2_idx, guest_prefix="UNL"):
    s = pmd.load_file("complex/vac.prmtop", "complex/vac.rst7", structure=True)
    guest_atoms = [a for a in s.atoms if a.residue.name == guest_prefix]
    G1 = f":{guest_prefix}@{guest_atoms[b1_idx].idx + 1}"
    G2 = f":{guest_prefix}@{guest_atoms[b2_idx].idx + 1}"
    return G1, G2


def add_dummy(G1, G2, data_dir_host, data_dir_guest, host_prefix="BCD", guest_prefix="UNL"): 
    structure = pmd.load_file("complex/vac.prmtop", "complex/vac.rst7", structure=True)
    aligned = align.zalign(structure, G1, G2)
    aligned.save("complex/aligned.prmtop", overwrite=True)
    aligned.save("complex/aligned.rst7", overwrite=True)

    structure = pmd.load_file("complex/aligned.prmtop", "complex/aligned.rst7", structure=True)
    structure = dummy.add_dummy(structure, residue_name="DM1", z=-6.0)
    structure = dummy.add_dummy(structure, residue_name="DM2", z=-9.0)
    structure = dummy.add_dummy(structure, residue_name="DM3", z=-11.2, y=2.2)
    structure.save("complex/aligned_with_dummy.prmtop", overwrite=True)
    structure.save("complex/aligned_with_dummy.rst7", overwrite=True)
    structure.save("complex/aligned_with_dummy.pdb", overwrite=True)

    dummy.write_dummy_frcmod(filepath="complex/dummy.frcmod")
    dummy.write_dummy_mol2(residue_name="DM1", filepath="complex/dm1.mol2")
    dummy.write_dummy_mol2(residue_name="DM2", filepath="complex/dm2.mol2")
    dummy.write_dummy_mol2(residue_name="DM3", filepath="complex/dm3.mol2")

# --- Final vacuum system with dummy atoms ---
    system = TLeap()
    system.output_path = "complex"
    system.pbc_type = None
    system.neutralize = False
    system.template_lines = [
        "source leaprc.gaff",
        "loadamberparams frcmod.ions1lm_126_tip3p",
        f"loadamberparams {data_dir_host}/{host_prefix}.frcmod",
        f"{host_prefix} = loadmol2 {data_dir_host}/{host_prefix}.mol2",
        f"loadamberparams {data_dir_guest}/{guest_prefix}.frcmod",
        f"{guest_prefix} = loadmol2 {data_dir_guest}/{guest_prefix}.mol2",
        "loadamberparams dummy.frcmod",
        "DM1 = loadmol2 dm1.mol2",
        "DM2 = loadmol2 dm2.mol2",
        "DM3 = loadmol2 dm3.mol2",
        "model = loadpdb aligned_with_dummy.pdb",
        "check model",
        f"savepdb model {host_prefix}-{guest_prefix}-dum.pdb",
        f"saveamberparm model {host_prefix}-{guest_prefix}-dum.prmtop {host_prefix}-{guest_prefix}-dum.rst7",
    ]
    system.build()
    return

def window_setup():
    attach_string = "0.00 0.01 0.0375 0.40 0.80 1.60 2.40 4.00 5.50 8.65 11.80 18.10 24.40 37.00 49.60 74.80 100.00"
    attach_fractions = [float(i) / 100 for i in attach_string.split()]
    initial_distance = 6.0
    pull_distances = np.arange(initial_distance, 14.5 + initial_distance, 0.5)
    release_fractions = []
    windows = [len(attach_fractions), len(pull_distances), len(release_fractions)]
    return windows, attach_fractions

def set_restraints(G1, G2, host_prefix= "BCD", guest_prefix="UNL"):

    windows, attach_fractions= window_setup()

    H = [f":{host_prefix}@75", f":{host_prefix}@33", f":{host_prefix}@138"]
    H1, H2 = [], []
    for i in range(7):
        if i == 6:
            H1.append([f":{host_prefix}@142", f":{host_prefix}@127", f":{host_prefix}@129", f":{host_prefix}@12"])
            H2.append([f":{host_prefix}@127", f":{host_prefix}@129", f":{host_prefix}@12", f":{host_prefix}@14"])
        else:
            H1.append([f":{host_prefix}@{16+i*21}", f":{host_prefix}@{1+i*21}", f":{host_prefix}@{3+i*21}", f":{host_prefix}@{33+i*21}"])
            H2.append([f":{host_prefix}@{1+i*21}", f":{host_prefix}@{3+i*21}", f":{host_prefix}@{33+i*21}", f":{host_prefix}@{35+i*21}"])

    D1, D2, D3 = ":DM1", ":DM2", ":DM3"

    structure = pmd.load_file(f"complex/{host_prefix}-{guest_prefix}-dum.prmtop", f"complex/{host_prefix}-{guest_prefix}-dum.rst7")

    # Static restraints keep host in a fixed position/orientation throughout all windows.
    static_restraints = []
    for masks, fc in [
        ([D1, H[0]], 10.0),
        ([D2, D1, H[0]], 100.0),
        ([D3, D2, D1, H[0]], 100.0),
        ([D1, H[0], H[1]], 100.0),
        ([D2, D1, H[0], H[1]], 100.0),
        ([D1, H[0], H[1], H[2]], 100.0),
    ]:
        static_restraints.append(
            restraints.static_DAT_restraint(
                restraint_mask_list=masks,
                num_window_list=windows,
                ref_structure=structure,
                force_constant=fc,
                amber_index=True,
            )
        )
    for masks in H1 + H2:
        static_restraints.append(
            restraints.static_DAT_restraint(
                restraint_mask_list=masks,
                num_window_list=windows,
                ref_structure=structure,
                force_constant=10.0,
                amber_index=True,
            )
        )

    # Guest restraints drive the APR coordinate: attach raises force constants,
    # pull translates the guest, release lowers force constants analytically.
    guest_restraints = []
    for mask1, mask2, mask3, attach_tgt, attach_fc, pull_tgt in [
        (D1,  G1,   None,  6.0,   10.0,  20.0),   # distance (Å, kcal/mol/Å²)
        (D2,  D1,   G1,    180.0, 100.0, 180.0),   # angle (deg, kcal/mol/rad²)
        (D1,  G1,   G2,    180.0, 100.0, 180.0),   # torsion
    ]:
        r = restraints.DAT_restraint()
        r.mask1 = mask1
        r.mask2 = mask2
        if mask3:
            r.mask3 = mask3
        r.topology = structure
        r.auto_apr = True
        r.continuous_apr = True
        r.amber_index = True
        r.attach["target"] = attach_tgt
        r.attach["fraction_list"] = attach_fractions
        r.attach["fc_final"] = attach_fc
        r.pull["target_final"] = pull_tgt
        r.pull["num_windows"] = windows[1]
        r.initialize()
        guest_restraints.append(r)

    # Create directories and write disang.rest for each window.
    window_list = create_window_list(guest_restraints)
    for window in window_list:
        if not os.path.isdir(f"windows/{window}"):
            os.makedirs(f"windows/{window}")

    all_restraints = static_restraints + guest_restraints
    for window in window_list:
        with open(f"windows/{window}/disang.rest", "w") as f:
            for r in all_restraints:
                line = amber_restraint_line(r, window)
                if line is not None:
                    f.write(line)

    save_restraints(guest_restraints, filepath="windows/restraints.json")

    return guest_restraints

def create_windows(G1, G2, host_prefix='BCD', guest_prefix='UNL'):

    guest_restraints = set_restraints(G1,G2, host_prefix, guest_prefix)

    window_list = create_window_list(guest_restraints)

    for window in window_list:
        if window[0] == "a":
            shutil.copy(f"complex/{host_prefix}-{guest_prefix}-dum.prmtop", f"windows/{window}/{host_prefix}-{guest_prefix}-dum.prmtop")
            shutil.copy(f"complex/{host_prefix}-{guest_prefix}-dum.rst7", f"windows/{window}/{host_prefix}-{guest_prefix}-dum.rst7")
        elif window[0] == "p":
            s = pmd.load_file(f"complex/{host_prefix}-{guest_prefix}-dum.prmtop", f"complex/{host_prefix}-{guest_prefix}-dum.rst7", structure=True)
            target_diff = (
                guest_restraints[0].phase["pull"]["targets"][int(window[1:])]
                - guest_restraints[0].pull["target_initial"]
            )
            print(f"Window {window}: translating guest {target_diff.magnitude:0.1f} Å.")
            for atom in s.atoms:
                if atom.residue.name == guest_prefix.upper():
                    atom.xz += target_diff.magnitude
            s.save(f"windows/{window}/{host_prefix}-{guest_prefix}-dum.prmtop", overwrite=True)
            s.save(f"windows/{window}/{host_prefix}-{guest_prefix}-dum.rst7", overwrite=True)

    return

def setup_datapoint(dp, data_dir_host, b1_idx, b2_idx, host_prefix="BCD", guest_prefix="UNL"):
    old = os.getcwd()
    os.chdir(dp)
    try:
        build_system(data_dir_host, dp, host_prefix, guest_prefix)
        G1, G2 = guest_masks(b1_idx, b2_idx, guest_prefix)   # needs complex/vac.prmtop
        print("masks:", G1, G2)
        add_dummy(G1, G2, data_dir_host, dp, host_prefix, guest_prefix)
        create_windows(G1, G2, host_prefix, guest_prefix)
    finally:
        os.chdir(old)
