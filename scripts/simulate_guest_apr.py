#!/usr/bin/env python3
"""Build and run the guest-only implicit-solvent p000 APR window."""

import argparse
import shlex
import subprocess
import sys
from pathlib import Path
from string import Template

import parmed as pmd
from rdkit import Chem

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "apr_pipeline"))

from apr_config import (
    CFG,
    get_attach_fractions,
    get_pmemd_command,
    get_pull_distances,
    get_radii_set,
)
from simulate_guest import validate_guest_inputs
from structure import find_carbon_backbone_chain, pick_backbone_alignment_atoms
from paprika import restraints
from paprika.build import align, dummy
from paprika.io import save_restraints
from paprika.restraints.amber import amber_restraint_line
from paprika.restraints.utils import create_window_list

GUEST = CFG["naming"]["guest_prefix"]
ATOM_TYPE = CFG["parametrize"]["atom_type"]
RUN = CFG["run"]
APR = CFG["apr"]
APR_RESTRAINTS = APR["restraints"]
TOPOLOGY_NAME = "guest_apr_p000"


def _run(command, cwd, log_path):
    result = subprocess.run(
        command, cwd=cwd, capture_output=True, text=True, check=False
    )
    log_path.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError(
            f"Command failed with exit status {result.returncode}: "
            f"{shlex.join(command)}\n{detail[-4000:]}"
        )


def _atom_masks(pdb_path, structure):
    molecule = Chem.MolFromPDBFile(str(pdb_path), removeHs=False)
    if molecule is None:
        raise ValueError(f"Could not read guest PDB with RDKit: {pdb_path}")

    backbone = find_carbon_backbone_chain(molecule)
    if len(backbone) < 2:
        raise ValueError("Guest needs at least two connected carbon atoms for APR alignment")
    first_idx, second_idx = pick_backbone_alignment_atoms(backbone)

    atom_names = []
    for atom_idx in (first_idx, second_idx):
        atom_info = molecule.GetAtomWithIdx(atom_idx).GetPDBResidueInfo()
        if atom_info is None:
            raise ValueError("Could not read PDB atom names for APR backbone atoms")
        atom_names.append(atom_info.GetName().strip())

    guest_atoms = [
        atom for atom in structure.atoms if atom.residue.name == GUEST.upper()
    ]
    for name in atom_names:
        matches = [atom for atom in guest_atoms if atom.name.strip() == name]
        if len(matches) != 1:
            raise ValueError(
                f"Expected exactly one {GUEST} atom named {name!r} in the Amber structure"
            )

    return tuple(f":{GUEST}@{name}" for name in atom_names)


def _leap_script(pdb_path, mol2_path, frcmod_path):
    paths = (pdb_path, mol2_path, frcmod_path)
    if any(any(char in str(path) for char in ('"', "\n", "\r")) for path in paths):
        raise ValueError("Input file paths must not contain quotes or newline characters")

    return "\n".join(
        [
            f"source leaprc.{ATOM_TYPE}",
            "loadamberparams frcmod.ions1lm_126_tip3p",
            f"set default PBRadii {get_radii_set()}",
            f'loadamberparams "{frcmod_path}"',
            f'{GUEST} = loadmol2 "{mol2_path}"',
            'loadamberparams "dummy.frcmod"',
            'DM1 = loadmol2 "dm1.mol2"',
            'DM2 = loadmol2 "dm2.mol2"',
            'DM3 = loadmol2 "dm3.mol2"',
            f'guest = loadpdb "{pdb_path}"',
            "check guest",
            f'saveamberparm guest {TOPOLOGY_NAME}.prmtop {TOPOLOGY_NAME}.rst7',
            "quit",
            "",
        ]
    )


def _build_topology(pdb_path, mol2_path, frcmod_path, output_dir):
    structure = pmd.load_file(str(pdb_path))
    g1, g2 = _atom_masks(pdb_path, structure)
    structure = align.zalign(structure, g1, g2)

    for residue_name, x, y, z in (
        ("DM1", 0.0, 0.0, APR["dummy_z"][0]),
        ("DM2", 0.0, 0.0, APR["dummy_z"][1]),
        ("DM3", 0.0, APR["dummy3_y"], APR["dummy_z"][2]),
    ):
        structure = dummy.add_dummy(
            structure, residue_name=residue_name, x=x, y=y, z=z
        )
    pdb_with_dummy = output_dir / "guest_with_dummy.pdb"
    structure.save(str(pdb_with_dummy), overwrite=True)

    dummy.write_dummy_frcmod(filepath=str(output_dir / "dummy.frcmod"))
    for residue_name in ("DM1", "DM2", "DM3"):
        dummy.write_dummy_mol2(
            residue_name=residue_name,
            filepath=str(output_dir / f"{residue_name.lower()}.mol2"),
        )

    leap = _leap_script(pdb_with_dummy.name, mol2_path, frcmod_path)
    (output_dir / "leap.in").write_text(leap, encoding="utf-8")
    _run(["tleap", "-f", "leap.in"], output_dir, output_dir / "tleap.log")


def _create_p000_restraints(topology, g1, g2):
    windows = [len(get_attach_fractions()), len(get_pull_distances()), 0]
    d1, d2 = ":DM1", ":DM2"
    restraints_for_guest = []

    # Keep the existing APR p000 mask definitions unchanged.
    for mask1, mask2, mask3, target, force_constant, pull_target in (
        (d1, g1, None, APR["initial_distance"], APR_RESTRAINTS["distance_fc"],
         APR["final_distance"]),
        (d2, d1, g1, APR_RESTRAINTS["angle_target"], APR_RESTRAINTS["angle_fc"],
         APR_RESTRAINTS["angle_target"]),
        (d1, g1, g2, APR_RESTRAINTS["torsion_target"],
         APR_RESTRAINTS["torsion_fc"], APR_RESTRAINTS["torsion_target"]),
    ):
        restraint = restraints.DAT_restraint()
        restraint.mask1 = mask1
        restraint.mask2 = mask2
        if mask3:
            restraint.mask3 = mask3
        restraint.topology = topology
        restraint.auto_apr = True
        restraint.continuous_apr = True
        restraint.amber_index = True
        restraint.attach["target"] = target
        restraint.attach["fraction_list"] = get_attach_fractions()
        restraint.attach["fc_final"] = force_constant
        restraint.pull["target_final"] = pull_target
        restraint.pull["num_windows"] = windows[1]
        restraint.initialize()
        restraints_for_guest.append(restraint)

    if "p000" not in create_window_list(restraints_for_guest):
        raise ValueError("pAPRika did not generate the expected p000 pull window")
    return restraints_for_guest


def _write_p000_restraints(output_dir, g1, g2):
    topology = pmd.load_file(
        str(output_dir / f"{TOPOLOGY_NAME}.prmtop"),
        str(output_dir / f"{TOPOLOGY_NAME}.rst7"),
    )
    guest_restraints = _create_p000_restraints(topology, g1, g2)
    lines = [
        amber_restraint_line(restraint, "p000")
        for restraint in guest_restraints
    ]
    (output_dir / "disang.rest").write_text(
        "".join(line for line in lines if line), encoding="utf-8"
    )
    save_restraints(
        guest_restraints, filepath=str(output_dir / "restraints.json")
    )


def _render_inputs(output_dir):
    template_dir = Path(CFG["directories"]["templates"])
    values = {
        "igb": RUN["igb"],
        "temperature": RUN["temperature"],
        "restraint_wt": RUN["restraint_wt"],
        "maxcyc": RUN["minimize_maxcyc"],
        "ncyc": RUN["minimize_ncyc"],
        "dt": RUN["dt"],
        "gamma_ln": RUN["gamma_ln"],
        "write_every": RUN["write_every"],
        "restart_every": RUN["restart_every"],
    }
    inputs = (
        ("minimize_nopbc.in", "minimize.in", {}),
        ("production_nopbc1.in", "equilibrate.in", {"nstlim": RUN["steps_equil"]}),
        ("production_nopbc2.in", "production.in", {"nstlim": RUN["steps_production"]}),
    )
    for template_name, output_name, extra in inputs:
        template = Template(
            (template_dir / template_name).read_text(encoding="utf-8")
        )
        (output_dir / output_name).write_text(
            template.substitute({**values, **extra}), encoding="utf-8"
        )


def _run_stages(output_dir, pmemd_command):
    pmemd = shlex.split(pmemd_command)
    stages = (
        ("minimize", "minimize.in", f"{TOPOLOGY_NAME}.rst7", "minimized.rst7"),
        ("equilibrate", "equilibrate.in", "minimized.rst7", "equilibrated.rst7"),
        ("production", "production.in", "equilibrated.rst7", "production.rst7"),
    )
    for name, input_file, start_file, restart_file in stages:
        command = [
            *pmemd,
            "-O",
            "-i", input_file,
            "-o", f"{name}.out",
            "-p", f"{TOPOLOGY_NAME}.prmtop",
            "-c", start_file,
            "-r", restart_file,
            "-ref", f"{TOPOLOGY_NAME}.rst7",
            "-inf", f"{name}.mdinfo",
        ]
        if name != "minimize":
            command.extend(
                [
                    "-x", f"{name}.nc",
                    "-v", f"{name}.mdvel",
                    "-e", f"{name}.mden",
                ]
            )
        _run(command, output_dir, output_dir / f"{name}.command.log")


def run_guest_apr_simulation(
    pdb_path, mol2_path, frcmod_path, output_dir, pmemd_command=None
):
    pdb_path = Path(pdb_path).expanduser().resolve()
    mol2_path = Path(mol2_path).expanduser().resolve()
    frcmod_path = Path(frcmod_path).expanduser().resolve()
    output_dir = Path(output_dir).expanduser().resolve()

    for path, label in (
        (pdb_path, "PDB"),
        (mol2_path, "MOL2"),
        (frcmod_path, "FRCMOD"),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"{label} file does not exist: {path}")
    validate_guest_inputs(pdb_path, mol2_path)

    output_dir.mkdir(parents=True, exist_ok=True)
    _build_topology(pdb_path, mol2_path, frcmod_path, output_dir)
    topology = pmd.load_file(
        str(output_dir / f"{TOPOLOGY_NAME}.prmtop"),
        str(output_dir / f"{TOPOLOGY_NAME}.rst7"),
    )
    g1, g2 = _atom_masks(pdb_path, topology)
    _write_p000_restraints(output_dir, g1, g2)
    _render_inputs(output_dir)
    _run_stages(output_dir, pmemd_command or get_pmemd_command())
    print(f"Guest p000 simulation complete. Output: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Build and run a guest-only implicit-solvent p000 APR window."
    )
    parser.add_argument("pdb", type=Path, help="Guest-only PDB coordinates")
    parser.add_argument(
        "mol2", type=Path, help="Parameterized guest MOL2 (GAFF atom types/charges)"
    )
    parser.add_argument(
        "--frcmod",
        type=Path,
        help="Guest FRCMOD (default: <PDB directory>/<guest_prefix>.frcmod)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=f"Output directory (default: <PDB directory>/{TOPOLOGY_NAME})",
    )
    parser.add_argument(
        "--pmemd-command",
        help="Override the configured PMEMD command (for example, 'pmemd')",
    )
    args = parser.parse_args()
    frcmod = args.frcmod or args.pdb.parent / f"{GUEST}.frcmod"
    output_dir = args.output_dir or args.pdb.parent / TOPOLOGY_NAME
    run_guest_apr_simulation(
        args.pdb, args.mol2, frcmod, output_dir, args.pmemd_command
    )


if __name__ == "__main__":
    main()
