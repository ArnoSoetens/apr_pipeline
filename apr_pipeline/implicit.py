"""Guest-only GB simulations of one data point, in <dp>/windows_implicit/<run>/.

    implicit_only   the free guest: tleap (guest.prmtop), minimize -> equilibrate -> production,
                    with data/inputs/implicit_only/*.in (no restraints)
    implicit_p000   the guest alone in the p000 APR window: dummy atoms, the guest's APR restraints
                    at p000 (disang.rest, restraints.json), minimize -> equilibrate -> production with
                    the APR window .in files in data/inputs/

Both read the data point's <guest>.pdb, .mol2 and .frcmod; the pmemd command is run.pmemd from the
config. A run that finishes writes the done marker (it is never redone); a step that fails writes
run_failed.txt (step, exit code, last line of its output). scripts/run_implicit.py picks the points.
The tleap scripts, dummy atoms, restraints and pmemd commands are those of the former
scripts/simulate_guest.py and simulate_guest_apr.py.
"""
import shlex
import shutil
import subprocess
import time
import socket
from pathlib import Path

from apr_config import CFG, get_attach_fractions, get_pmemd_command, get_pull_distances, get_radii_set
from datapoints import (DONE_MARKER, GUEST_FRCMOD, GUEST_MOL2, GUEST_PDB, IMPLICIT_DIR, IMPLICIT_ONLY,
                        IMPLICIT_P000, RUN_FAILED_FILE, TRAJECTORY)

RUNS = (IMPLICIT_ONLY, IMPLICIT_P000)
GUEST = CFG["naming"]["guest_prefix"]
ATOM_TYPE = CFG["parametrize"]["atom_type"]
APR = CFG["apr"]
APR_RESTRAINTS = APR["restraints"]
TOPOLOGY = {IMPLICIT_ONLY: "guest", IMPLICIT_P000: "guest_apr_p000"}   # <name>.prmtop / .rst7


class StepFailed(Exception):
    def __init__(self, step, code, output, detail=""):
        super().__init__(f"{step} failed with exit code {code}" + (f": {detail}" if detail else ""))
        self.step, self.code, self.output = step, code, output


# --------------------------------------------------------------------------
# where things are
# --------------------------------------------------------------------------

def run_dir(dp, run):
    return Path(dp) / IMPLICIT_DIR / run


def finished(dp, run):
    """Done marker + production.nc in the run folder."""
    folder = run_dir(dp, run)
    return (folder / DONE_MARKER).is_file() and (folder / TRAJECTORY).is_file()


def guest_files(dp):
    dp = Path(dp)
    return dp / GUEST_PDB, dp / GUEST_MOL2, dp / GUEST_FRCMOD


def missing_inputs(dp):
    """The guest files a run needs that the data point does not have (empty: ready)."""
    return [path.name for path in guest_files(dp) if not path.is_file()]


def clean_plan(folder):
    """Files a finished run folder does not need any more: everything but .prmtop files,
    production.nc and the done marker (subfolders included)."""
    folder = Path(folder)
    keep = {TRAJECTORY, DONE_MARKER}
    remove = []
    for path in sorted(folder.rglob("*")):
        if path.is_dir() and not path.is_symlink():
            continue
        if path.parent == folder and (path.suffix == ".prmtop" or path.name in keep):
            continue
        remove.append(path)
    return remove


# --------------------------------------------------------------------------
# the two runs (tleap and inputs as in the former simulate_guest*.py)
# --------------------------------------------------------------------------

def _mol2_atoms(path):
    atoms, in_atoms = [], False
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
        elif line.startswith("@<TRIPOS>"):
            in_atoms = False
        elif in_atoms and line.strip():
            fields = line.split()
            if len(fields) >= 8:
                atoms.append((fields[1], fields[7]))
    return atoms


def _pdb_atoms(path):
    atoms, residues = [], set()
    for line in Path(path).read_text(encoding="utf-8").splitlines():
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
        raise ValueError(f"PDB residue name {pdb_resname!r} must match the MOL2 substructure "
                         f"name(s) {sorted(mol2_resnames)!r}")
    mol2_names = [name for name, _ in mol2_atoms]
    if len(pdb_atoms) != len(mol2_names) or sorted(pdb_atoms) != sorted(mol2_names):
        raise ValueError("PDB and MOL2 must contain the same atom names and atom count; "
                         "check for missing hydrogens or mismatched atom naming")
    if len(set(pdb_atoms)) != len(pdb_atoms):
        raise ValueError("Guest atom names must be unique within the PDB residue")


def _check_paths(*paths):
    if any(any(char in str(path) for char in ('"', "\n", "\r")) for path in paths):
        raise ValueError("Input file paths must not contain quotes or newline characters")


def _leap_only(pdb_path, mol2_path, frcmod_path):
    _check_paths(pdb_path, mol2_path, frcmod_path)
    return "\n".join([
        f"source leaprc.{ATOM_TYPE}",
        "loadamberparams frcmod.ions1lm_126_tip3p",
        f"set default PBRadii {get_radii_set()}",
        f'loadamberparams "{frcmod_path}"',
        f'{GUEST} = loadmol2 "{mol2_path}"',
        f'guest = loadpdb "{pdb_path}"',
        "check guest",
        "savepdb guest prepared_guest.pdb",
        "saveamberparm guest guest.prmtop guest.rst7",
        "quit",
        "",
    ])


def _leap_p000(pdb_path, mol2_path, frcmod_path):
    _check_paths(pdb_path, mol2_path, frcmod_path)
    name = TOPOLOGY[IMPLICIT_P000]
    return "\n".join([
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
        f"saveamberparm guest {name}.prmtop {name}.rst7",
        "quit",
        "",
    ])


def _tool(command, cwd, log_name, step):
    """Run tleap/parmchk2; output to <log_name>; StepFailed on a non-zero exit."""
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True, check=False)
    (Path(cwd) / log_name).write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        detail = (result.stderr or result.stdout).strip().splitlines()
        raise StepFailed(step, result.returncode, log_name, detail[-1] if detail else "")


def anchor_masks(dp, structure):
    """Amber masks of G1/G2: the anchors in guest_input.json (or the chain ends, lower one first,
    as placement picks them), by atom name in the guest PDB."""
    from rdkit import Chem

    from guest_input import anchors
    from placement import order_anchors
    from structure import backbone_anchors

    pdb_path = Path(dp) / GUEST_PDB
    molecule = Chem.MolFromPDBFile(str(pdb_path), removeHs=False)
    if molecule is None:
        raise ValueError(f"Could not read guest PDB with RDKit: {pdb_path}")
    chosen = anchors(dp)
    if chosen is None:
        chosen = order_anchors(molecule.GetConformer().GetPositions(),
                               *backbone_anchors(molecule, source=str(pdb_path)))
    names = []
    for index in chosen:
        info = molecule.GetAtomWithIdx(index).GetPDBResidueInfo()
        if info is None:
            raise ValueError("Could not read PDB atom names for the APR anchor atoms")
        names.append(info.GetName().strip())
    guest_atoms = [atom for atom in structure.atoms if atom.residue.name == GUEST.upper()]
    for name in names:
        if len([atom for atom in guest_atoms if atom.name.strip() == name]) != 1:
            raise ValueError(f"Expected exactly one {GUEST} atom named {name!r} in the Amber structure")
    return tuple(f":{GUEST}@{name}" for name in names)


def _p000_restraints(topology, g1, g2):
    from paprika import restraints
    from paprika.restraints.utils import create_window_list

    windows = [len(get_attach_fractions()), len(get_pull_distances()), 0]
    d1, d2 = ":DM1", ":DM2"
    guest_restraints = []
    for mask1, mask2, mask3, target, force_constant, pull_target in (
        (d1, g1, None, APR["initial_distance"], APR_RESTRAINTS["distance_fc"], APR["final_distance"]),
        (d2, d1, g1, APR_RESTRAINTS["angle_target"], APR_RESTRAINTS["angle_fc"],
         APR_RESTRAINTS["angle_target"]),
        (d1, g1, g2, APR_RESTRAINTS["torsion_target"], APR_RESTRAINTS["torsion_fc"],
         APR_RESTRAINTS["torsion_target"]),
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
        guest_restraints.append(restraint)
    if "p000" not in create_window_list(guest_restraints):
        raise ValueError("pAPRika did not generate the expected p000 pull window")
    return guest_restraints


def _build_only(dp, folder):
    pdb, mol2, frcmod = guest_files(dp)
    (folder / "leap.in").write_text(_leap_only(pdb, mol2, frcmod), encoding="utf-8")
    _tool(["tleap", "-f", "leap.in"], folder, "tleap.log", "tleap")


def _build_p000(dp, folder):
    import parmed as pmd
    from paprika.build import align, dummy
    from paprika.io import save_restraints
    from paprika.restraints.amber import amber_restraint_line

    pdb, mol2, frcmod = guest_files(dp)
    structure = pmd.load_file(str(pdb))
    g1, g2 = anchor_masks(dp, structure)
    structure = align.zalign(structure, g1, g2)
    for residue_name, x, y, z in (("DM1", 0.0, 0.0, APR["dummy_z"][0]),
                                  ("DM2", 0.0, 0.0, APR["dummy_z"][1]),
                                  ("DM3", 0.0, APR["dummy3_y"], APR["dummy_z"][2])):
        structure = dummy.add_dummy(structure, residue_name=residue_name, x=x, y=y, z=z)
    with_dummy = folder / "guest_with_dummy.pdb"
    structure.save(str(with_dummy), overwrite=True)
    dummy.write_dummy_frcmod(filepath=str(folder / "dummy.frcmod"))
    for residue_name in ("DM1", "DM2", "DM3"):
        dummy.write_dummy_mol2(residue_name=residue_name, filepath=str(folder / f"{residue_name.lower()}.mol2"))
    (folder / "leap.in").write_text(_leap_p000(with_dummy.name, mol2, frcmod), encoding="utf-8")
    _tool(["tleap", "-f", "leap.in"], folder, "tleap.log", "tleap")

    name = TOPOLOGY[IMPLICIT_P000]
    topology = pmd.load_file(str(folder / f"{name}.prmtop"), str(folder / f"{name}.rst7"))
    g1, g2 = anchor_masks(dp, topology)
    guest_restraints = _p000_restraints(topology, g1, g2)
    lines = [amber_restraint_line(restraint, "p000") for restraint in guest_restraints]
    (folder / "disang.rest").write_text("".join(line for line in lines if line), encoding="utf-8")
    save_restraints(guest_restraints, filepath=str(folder / "restraints.json"))


def _stages(run):
    """(step, .in file, start coordinates, restart file, trajectory name or None, use -ref)."""
    inputs = Path(CFG["directories"]["inputs"])
    if run == IMPLICIT_ONLY:
        only = inputs / "implicit_only"
        return (
            ("minimize", only / "minimize.in", "guest.rst7", "minimized.rst7", None, False),
            ("equilibrate", only / "equilibrate.in", "minimized.rst7", "equilibrated.rst7", "equilibration", False),
            ("production", only / "production.in", "equilibrated.rst7", "production.rst7", "production", False),
        )
    top = TOPOLOGY[IMPLICIT_P000]
    return (
        ("minimize", inputs / "minimize_nopbc.in", f"{top}.rst7", "minimized.rst7", None, True),
        ("equilibrate", inputs / "production_nopbc1.in", "minimized.rst7", "equilibrated.rst7", "equilibrate", True),
        ("production", inputs / "production_nopbc2.in", "equilibrated.rst7", "production.rst7", "production", True),
    )


def _pmemd(run, folder):
    pmemd = shlex.split(get_pmemd_command())
    top = TOPOLOGY[run]
    for step, input_file, start, restart, trajectory, reference in _stages(run):
        command = [*pmemd, "-O", "-i", str(input_file), "-o", f"{step}.out", "-p", f"{top}.prmtop",
                   "-c", start, "-r", restart]
        if reference:
            command += ["-ref", f"{top}.rst7"]
        command += ["-inf", f"{step}.mdinfo"]
        if trajectory:
            command += ["-x", f"{trajectory}.nc", "-v", f"{trajectory}.mdvel", "-e", f"{trajectory}.mden"]
        result = subprocess.run(command, cwd=folder, capture_output=True, text=True, check=False)
        (folder / f"{step}.command.log").write_text(result.stdout + result.stderr, encoding="utf-8")
        if result.returncode:
            raise StepFailed(step, result.returncode, f"{step}.out")


def _record_failure(folder, step, code, last):
    (folder / RUN_FAILED_FILE).write_text(
        f"date: {time.strftime('%F %T')}\nhost: {socket.gethostname()}\nstep: {step}\n"
        f"exit code: {code}\n{last}\n")


def _last_line(path):
    lines = [line for line in path.read_text(errors="replace").splitlines() if line.strip()] \
        if path.is_file() else []
    return lines[-1] if lines else "(no output)"


def run(dp, run_name, clean=False):
    """Run one guest-only simulation of data point dp from scratch.

    The run folder is emptied first (a finished run is never passed here: run_implicit.py skips it).
    Raises StepFailed (after writing run_failed.txt) or ValueError/FileNotFoundError for bad inputs.
    """
    if run_name not in RUNS:
        raise ValueError(f"unknown guest-only run {run_name!r} (one of {RUNS})")
    dp = Path(dp).resolve()
    missing = missing_inputs(dp)
    if missing:
        raise FileNotFoundError(f"{dp.name} has no {', '.join(missing)} (parametrize first)")
    pdb, mol2, _ = guest_files(dp)
    validate_guest_inputs(pdb, mol2)

    folder = run_dir(dp, run_name)
    if folder.exists():
        shutil.rmtree(folder)                 # an unfinished earlier attempt: start again
    folder.mkdir(parents=True)
    try:
        if run_name == IMPLICIT_ONLY:
            _build_only(dp, folder)
        else:
            _build_p000(dp, folder)
        _pmemd(run_name, folder)
    except StepFailed as failure:
        _record_failure(folder, failure.step, failure.code,
                        f"last line of {failure.output}: {_last_line(folder / failure.output)}")
        raise
    except Exception as error:                # e.g. bad anchors: no command ran
        _record_failure(folder, "build", "-", f"error: {error}")
        raise
    (folder / DONE_MARKER).write_text("")
    if clean:
        for path in clean_plan(folder):
            path.unlink()
        for sub in sorted((p for p in folder.rglob("*") if p.is_dir()), reverse=True):
            if not any(sub.iterdir()):
                sub.rmdir()
