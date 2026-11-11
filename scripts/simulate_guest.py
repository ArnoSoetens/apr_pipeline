#!/usr/bin/env python3
"""Build and run a guest-only Amber simulation in GB implicit solvent."""

import argparse
import shlex
import subprocess
import sys
from pathlib import Path
from string import Template

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "apr_pipeline"))

from apr_config import CFG, get_pmemd_command, get_radii_set


def _mol2_atoms(path):
    atoms = []
    in_atom_section = False
    with path.open(encoding="utf-8") as mol2_file:
        for line in mol2_file:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
            elif line.startswith("@<TRIPOS>"):
                in_atom_section = False
            elif in_atom_section and line.strip():
                fields = line.split()
                if len(fields) >= 8:
                    atoms.append((fields[1], fields[7]))
    return atoms


def _pdb_atoms(path):
    atoms = []
    residues = set()
    with path.open(encoding="utf-8") as pdb_file:
        for line in pdb_file:
            if line.startswith(("ATOM  ", "HETATM")):
                atoms.append(line[12:16].strip())
                residues.add((line[17:20].strip(), line[21:22], line[22:27]))
    return atoms, residues


def validate_guest_inputs(pdb_path, mol2_path):
    """Check that the PDB coordinates map to one residue in the MOL2 template."""
    pdb_atoms, residues = _pdb_atoms(pdb_path)
    mol2_atoms = _mol2_atoms(mol2_path)

    if not pdb_atoms:
        raise ValueError(f"No ATOM/HETATM records found in {pdb_path}")
    if len(residues) != 1:
        raise ValueError(f"{pdb_path} must contain exactly one guest residue")
    if not mol2_atoms:
        raise ValueError(f"No atoms found in the MOL2 atom section of {mol2_path}")

    pdb_resname = next(iter(residues))[0]
    mol2_resnames = {resname for _, resname in mol2_atoms}
    if mol2_resnames != {pdb_resname}:
        raise ValueError(
            f"PDB residue name {pdb_resname!r} must match the MOL2 substructure "
            f"name(s) {sorted(mol2_resnames)!r}"
        )

    mol2_names = [name for name, _ in mol2_atoms]
    if len(pdb_atoms) != len(mol2_names) or sorted(pdb_atoms) != sorted(mol2_names):
        raise ValueError(
            "PDB and MOL2 must contain the same atom names and atom count; "
            "check for missing hydrogens or mismatched atom naming"
        )
    if len(set(pdb_atoms)) != len(pdb_atoms):
        raise ValueError("Guest atom names must be unique within the PDB residue")


def _amber_value(value):
    if isinstance(value, bool):
        return ".true." if value else ".false."
    return str(value)


def render_inputs(output_dir):
    run = CFG["run"]
    values = {
        "igb": run["igb"],
        "temperature": run["temperature"],
        "dt": run["dt"],
        "gamma_ln": run["gamma_ln"],
        "write_every": run["write_every"],
        "restart_every": run["restart_every"],
        "minimize_maxcyc": run["minimize_maxcyc"],
        "minimize_ncyc": run["minimize_ncyc"],
        "steps_equil": run["steps_equil"],
        "steps_production": run["steps_production"],
    }
    templates = {
        "minimize.in": "guest_minimize.in",
        "equilibrate.in": "guest_equilibrate.in",
        "production.in": "guest_production.in",
    }
    template_dir = SCRIPT_DIR / "templates"
    for output_name, template_name in templates.items():
        template = Template((template_dir / template_name).read_text(encoding="utf-8"))
        (output_dir / output_name).write_text(
            template.substitute({key: _amber_value(value) for key, value in values.items()}),
            encoding="utf-8",
        )


def _run(command, cwd, log_path):
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True, check=False)
    log_path.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError(
            f"Command failed with exit status {result.returncode}: "
            f"{shlex.join(command)}\n{detail[-4000:]}"
        )


def _leap_script(pdb_path, mol2_path, frcmod_path):
    paths = (pdb_path, mol2_path, frcmod_path)
    if any(any(char in str(path) for char in ('"', "\n", "\r")) for path in paths):
        raise ValueError("Input file paths must not contain quotes or newline characters")

    return "\n".join(
        [
            f"source leaprc.{CFG['parametrize']['atom_type']}",
            "loadamberparams frcmod.ions1lm_126_tip3p",
            f"set default PBRadii {get_radii_set()}",
            f'loadamberparams "{frcmod_path}"',
            f'{CFG["naming"]["guest_prefix"]} = loadmol2 "{mol2_path}"',
            f'guest = loadpdb "{pdb_path}"',
            "check guest",
            "savepdb guest prepared_guest.pdb",
            "saveamberparm guest guest.prmtop guest.rst7",
            "quit",
            "",
        ]
    )


def run_guest_simulation(
    pdb_path, mol2_path, output_dir, frcmod_path=None, pmemd_command=None
):
    pdb_path = Path(pdb_path).expanduser().resolve()
    mol2_path = Path(mol2_path).expanduser().resolve()
    output_dir = Path(output_dir).expanduser().resolve()
    frcmod_path = Path(frcmod_path).expanduser().resolve() if frcmod_path else None

    for path, label in ((pdb_path, "PDB"), (mol2_path, "MOL2")):
        if not path.is_file():
            raise FileNotFoundError(f"{label} file does not exist: {path}")
    if frcmod_path is not None and not frcmod_path.is_file():
        raise FileNotFoundError(f"FRCMOD file does not exist: {frcmod_path}")
    validate_guest_inputs(pdb_path, mol2_path)

    output_dir.mkdir(parents=True, exist_ok=True)
    if frcmod_path is None:
        frcmod_path = output_dir / "guest.frcmod"
        _run(
            [
                "parmchk2",
                "-i",
                str(mol2_path),
                "-f",
                "mol2",
                "-o",
                str(frcmod_path),
                "-s",
                str(CFG["parametrize"]["parmchk_set"]),
            ],
            output_dir,
            output_dir / "parmchk2.log",
        )

    (output_dir / "leap.in").write_text(
        _leap_script(pdb_path, mol2_path, frcmod_path), encoding="utf-8"
    )
    _run(["tleap", "-f", "leap.in"], output_dir, output_dir / "tleap.log")
    render_inputs(output_dir)

    pmemd = shlex.split(pmemd_command or get_pmemd_command())
    stages = [
        (
            "minimize",
            "minimize.in",
            "guest.rst7",
            "minimized.rst7",
            None,
        ),
        (
            "equilibrate",
            "equilibrate.in",
            "minimized.rst7",
            "equilibrated.rst7",
            "equilibration",
        ),
        (
            "production",
            "production.in",
            "equilibrated.rst7",
            "production.rst7",
            "production",
        ),
    ]

    for name, input_file, start_file, restart_file, trajectory in stages:
        command = [
            *pmemd,
            "-O",
            "-i",
            input_file,
            "-o",
            f"{name}.out",
            "-p",
            "guest.prmtop",
            "-c",
            start_file,
            "-r",
            restart_file,
            "-inf",
            f"{name}.mdinfo",
        ]
        if trajectory is not None:
            command.extend(
                [
                    "-x",
                    f"{trajectory}.nc",
                    "-v",
                    f"{trajectory}.mdvel",
                    "-e",
                    f"{trajectory}.mden",
                ]
            )
        _run(command, output_dir, output_dir / f"{name}.command.log")

    print(f"Guest simulation complete. Output: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Build and simulate a single guest with Amber GB implicit solvent."
    )
    parser.add_argument("pdb", type=Path, help="Guest-only PDB coordinates")
    parser.add_argument("mol2", type=Path, help="Parameterized guest MOL2 (GAFF atom types/charges)")
    parser.add_argument(
        "--frcmod",
        type=Path,
        help="Optional guest FRCMOD; parmchk2 generates one from the MOL2 if omitted",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory (default: <PDB stem>_implicit in the current directory)",
    )
    parser.add_argument(
        "--pmemd-command",
        help="Override the configured PMEMD command (for example, 'pmemd' for serial runs)",
    )
    args = parser.parse_args()
    output_dir = args.output_dir or Path(f"{args.pdb.stem}_implicit")
    run_guest_simulation(
        args.pdb, args.mol2, output_dir, args.frcmod, args.pmemd_command
    )


if __name__ == "__main__":
    main()
