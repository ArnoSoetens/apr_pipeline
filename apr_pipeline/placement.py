"""
guest_placement.py

Place a guest molecule in a host cavity, check for clashes, and write the
host+guest PDB files in the format the downstream steps expect.

Typical use
-----------
    from guest_placement import place_guest_in_cavity, generate_pdb, StageFailure

    guest_aligned, R, midpoint = place_guest_in_cavity(guest, b1_idx=0, b2_idx=3)
    try:
        out = generate_pdb(host, guest_aligned, out_dir="isomers/CCCCF")
    except StageFailure as e:
        print("skipped:", e.reason)   # e.g. clash -> record in status.json, move on
    else:
        print(out["guest_pdb"], out["complex_pdb"])

Files written into out_dir by generate_pdb():
    <guest_prefix>.pdb     guest only:  UNL residue set to 2, no TER
    <complex_prefix>.pdb   host+guest:  TER inserted before the UNL block,
                                        UNL residue set to 2
Everything else in the PDB text is left exactly as RDKit wrote it.
"""
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.spatial.distance import cdist

__all__ = [
    "StageFailure",
    "rotation_matrix_from_vectors",
    "place_guest_in_cavity",
    "clash_report",
    "radial_containment",
    "Mol_to_PDB",
    "convert",
    "renumber_guest",
    "generate_pdb",
]


class StageFailure(Exception):
    """An expected, recordable failure of a pipeline stage (clash, bad PDB, ...).

    A driver loop should catch this, store `reason` in the datapoint's
    status.json, and move on to the next datapoint instead of crashing.
    """

    def __init__(self, reason, stage=None, details=None):
        super().__init__(reason)
        self.reason = reason
        self.stage = stage
        self.details = details or {}


# --------------------------------------------------------------------------
# Geometry: placement, clash and containment checks
# --------------------------------------------------------------------------

def rotation_matrix_from_vectors(a, b):
    '''Rotation matrix R such that R @ a is parallel to b (both unit vectors).'''
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    if s < 1e-8:
        if c > 0:
            return np.eye(3)
        perp = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        axis = np.cross(a, perp)
        axis = axis / np.linalg.norm(axis)
        return 2 * np.outer(axis, axis) - np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s ** 2))


def place_guest_in_cavity(guest_mol, b1_idx, b2_idx, cavity_center=(0, 0, 0), cavity_axis=(0, 0, 1)):
    coords = guest_mol.GetConformer().GetPositions()

    backbone_axis = coords[b2_idx] - coords[b1_idx]
    backbone_axis = backbone_axis / np.linalg.norm(backbone_axis)
    backbone_midpoint = (coords[b1_idx] + coords[b2_idx]) / 2.0

    R = rotation_matrix_from_vectors(backbone_axis, cavity_axis)
    new_coords = (coords - backbone_midpoint) @ R.T + cavity_center

    new_mol = Chem.Mol(guest_mol)
    conf = new_mol.GetConformer()
    for i in range(new_mol.GetNumAtoms()):
        x, y, z = new_coords[i]
        conf.SetAtomPosition(i, (float(x), float(y), float(z)))
    return new_mol, R, backbone_midpoint


def clash_report(host_mol, guest_mol, clash_cutoff=1.5):
    hc = host_mol.GetConformer().GetPositions()
    gc = guest_mol.GetConformer().GetPositions()
    dmat = cdist(hc, gc)
    i, j = np.unravel_index(dmat.argmin(), dmat.shape)
    return {
        'min_distance': float(dmat.min()),
        'closest_host_atom': int(i),
        'closest_guest_atom': int(j),
        'n_pairs_under_cutoff': int((dmat < clash_cutoff).sum()),
        'cutoff_used': clash_cutoff,
    }


def radial_containment(guest_mol, cavity_center, cavity_axis):
    coords = guest_mol.GetConformer().GetPositions()
    rel = coords - cavity_center
    along = np.outer(rel @ cavity_axis, cavity_axis)
    perp = rel - along
    radial = np.linalg.norm(perp, axis=1)
    return {'min': float(radial.min()), 'max': float(radial.max()), 'mean': float(radial.mean())}


# --------------------------------------------------------------------------
# PDB writing
# --------------------------------------------------------------------------

def Mol_to_PDB(host, guest_aligned, complex_prefix="BCD-UNL", guest_prefix="UNL",
               out_dir=".", flavor=2):
    """Write the guest alone and the host+guest complex as PDB files.

    flavor is passed to RDKit; 2 means "write no CONECT records".
    Returns (guest_path, complex_path) as Path objects.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    guest_path = out_dir / f"{guest_prefix}.pdb"
    complex_path = out_dir / f"{complex_prefix}.pdb"

    Chem.MolToPDBFile(guest_aligned, str(guest_path), flavor=flavor)
    combined = Chem.CombineMols(host, guest_aligned)
    Chem.MolToPDBFile(combined, str(complex_path), flavor=flavor)
    return guest_path, complex_path


# --------------------------------------------------------------------------
# PDB post-processing: residue number (and TER) for the guest
# --------------------------------------------------------------------------



def eol_of(line: str) -> str:
    """Return the line ending of `line` ('\\n', '\\r\\n' or '' if none)."""
    return line[len(line.rstrip("\r\n")):]


def convert(lines, resname="UNL", resnum=2, add_ter=True):
    """Return (new_lines, n_ter_added, n_renumbered).

    A TER is only added at the start of a resname block when add_ter is True,
    some ATOM/HETATM record already came earlier in the file (so a guest-only
    file never gets one), and the previous line isn't already a TER.
    """
    RESNAME_COLS = slice(17, 20)   # PDB cols 18-20
    RESNUM_COLS = slice(22, 26)    # PDB cols 23-26 (right-justified)

    if not 0 <= resnum <= 9999:
        raise ValueError("resnum must fit in the 4-character PDB residue field (0-9999)")

    out = []
    in_block = False      # inside a run of matching HETATM lines?
    seen_atoms = False    # any ATOM/HETATM line earlier in the file?
    n_ter = 0
    n_renum = 0

    for line in lines:
        is_target = (
            line.startswith("HETATM")
            and len(line.rstrip("\r\n")) >= RESNUM_COLS.stop
            and line[RESNAME_COLS].strip() == resname
        )

        if is_target:
            if (add_ter and not in_block and seen_atoms
                    and not (out and out[-1].startswith("TER"))):
                out.append("TER  " + (eol_of(line) or "\n"))
                n_ter += 1
            line = line[:RESNUM_COLS.start] + f"{resnum:>4d}" + line[RESNUM_COLS.stop:]
            n_renum += 1
            in_block = True
        else:
            in_block = False

        if line.startswith(("ATOM", "HETATM")):
            seen_atoms = True
        out.append(line)

    return out, n_ter, n_renum


def renumber_guest(input_path, output_path=None, resnum=2, resname="UNL", add_ter=False):
    """Set the residue number of the `resname` lines in a PDB file.

    Defaults are for a guest-only file: residue 2, no TER. Pass add_ter=True
    for a host+guest file. Writes to output_path, or overwrites input_path if
    output_path is None. Returns (n_renumbered, n_ter_added).
    Raises ValueError if no matching HETATM lines exist (nothing is written).
    """
    input_path = Path(input_path)
    output_path = Path(output_path) if output_path else input_path

    # newline="" -> no newline translation, so line endings survive untouched
    with open(input_path, newline="") as f:
        lines = f.readlines()

    new_lines, n_ter, n_renum = convert(lines, resname, resnum, add_ter)
    if n_renum == 0:
        raise ValueError(f"no HETATM lines with residue name {resname!r} in {input_path}")

    with open(output_path, "w", newline="") as f:
        f.writelines(new_lines)
    return n_renum, n_ter


# --------------------------------------------------------------------------
# The stage: clash check -> write PDBs -> fix residue numbering
# --------------------------------------------------------------------------

def generate_pdb(host, guest_aligned, out_dir=".", complex_prefix="BCD-UNL", guest_prefix="UNL",
                 resnum=2, resname="UNL", flavor=2, clash_cutoff=1.5):
    """Write <guest_prefix>.pdb and <complex_prefix>.pdb into out_dir.

    1. Clash check: if any host-guest atom pair is closer than clash_cutoff
       (Angstrom), raise StageFailure BEFORE anything is written.
       Pass clash_cutoff=None to skip the check.
    2. Write both files with RDKit (see Mol_to_PDB).
    3. Set the guest's residue number to `resnum` in both files, and insert a
       TER before the guest in the complex file (the guest-only file gets none).

    Returns {"guest_pdb": Path, "complex_pdb": Path, "clash": report or None}.
    Raises StageFailure on a clash or if the guest residue can't be found.
    """
    report = None
    if clash_cutoff is not None:
        report = clash_report(host, guest_aligned, clash_cutoff)
        if report["n_pairs_under_cutoff"] > 0:
            raise StageFailure(
                f"clash: {report['n_pairs_under_cutoff']} host-guest atom pair(s) closer than "
                f"{clash_cutoff} A (min {report['min_distance']:.2f} A between host atom "
                f"{report['closest_host_atom']} and guest atom {report['closest_guest_atom']})",
                stage="generate_pdb",
                details=report,
            )

    guest_pdb, complex_pdb = Mol_to_PDB(host, guest_aligned, complex_prefix, guest_prefix,
                                        out_dir=out_dir, flavor=flavor)
    try:
        renumber_guest(guest_pdb, resnum=resnum, resname=resname, add_ter=False)
        renumber_guest(complex_pdb, resnum=resnum, resname=resname, add_ter=True)
    except ValueError as e:
        raise StageFailure(f"PDB post-processing failed: {e}", stage="generate_pdb") from e

    return {"guest_pdb": guest_pdb, "complex_pdb": complex_pdb, "clash": report}, print("PDB succesfully created")