"""guest_input.json: what a data point's setup stages hand to each other.

    {
      "smiles": "[O-]C(=O)C(F)(F)C(F)(F)C(F)(F)F",  from <name>.smi; null when the guest started as a PDB
      "net_charge": -1,                             the charge parametrize uses. Filled in while null:
                                                    with config parametrize.net_charge = "auto" from the
                                                    SMILES (a PDB-only guest stays null: set it by hand),
                                                    with an integer there that integer
      "anchors": [1, 5]                             G1, G2: the restrained guest atoms, as the atom serial
                                                    numbers in <guest_prefix>.pdb (column 2). null:
                                                    placement picks the carbon chain ends, G1 the lower
                                                    one, and writes them here. Set by hand: placement
                                                    turns G1->G2 onto +cavity_axis and keeps them.
    }

Every setup stage starts with ensure_guest_input(), so the file exists as soon as any stage has
touched the data point. Every value can be edited by hand; the next stage reads what is in the file,
and smiles/net_charge/anchors are only filled in automatically while they are still null, so an
edit is never overwritten. To apply a changed config net_charge to a point that already has a
charge, set its net_charge back to null.
"""
import json
from pathlib import Path

from apr_config import CFG

GUEST_INPUT = "guest_input.json"
GUEST_PDB = f'{CFG["naming"]["guest_prefix"]}.pdb'
DEFAULTS = {"smiles": None, "net_charge": None, "anchors": None}


def read_guest_input(dp):
    """The data point's guest_input.json as a dict (all keys present; missing file -> defaults)."""
    path = Path(dp) / GUEST_INPUT
    if not path.is_file():
        return dict(DEFAULTS)
    try:
        values = json.loads(path.read_text())
    except json.JSONDecodeError as e:
        raise ValueError(f"{path} is not valid JSON: {e}") from e
    if not isinstance(values, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return {**DEFAULTS, **values}


def update_guest_input(dp, **values):
    """Set some values in guest_input.json (creating it if needed), keep the others."""
    unknown = set(values) - set(DEFAULTS)
    if unknown:
        raise KeyError(f"unknown guest_input.json key(s): {', '.join(sorted(unknown))}")
    merged = {**read_guest_input(dp), **values}
    (Path(dp) / GUEST_INPUT).write_text(json.dumps(merged, indent=2) + "\n")
    return merged


def smiles_charge(smiles):
    """Formal charge of a SMILES, or None when RDKit cannot read it."""
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    return None if mol is None else Chem.GetFormalCharge(mol)


def ensure_guest_input(dp):
    """Create guest_input.json if the data point has none; an existing file is left alone.

    smiles comes from <name>.smi when the folder has one (otherwise null), net_charge from
    default_charge(); anchors stay null until placement.
    """
    dp = Path(dp)
    if (dp / GUEST_INPUT).is_file():
        return read_guest_input(dp)
    smi = dp / f"{dp.name}.smi"
    smiles = smi.read_text().strip() if smi.is_file() else None
    return update_guest_input(dp, smiles=smiles, net_charge=default_charge(smiles))


def default_charge(smiles):
    """What an empty net_charge is filled with: the config integer (parametrize.net_charge), or with
    "auto" the SMILES' formal charge (None without a SMILES)."""
    configured = CFG["parametrize"]["net_charge"]
    if configured != "auto":
        return configured
    return smiles_charge(smiles) if smiles else None


def net_charge(dp):
    """The guest charge parametrize uses: this point's guest_input.json net_charge.

    While it is null it is filled in first (default_charge) and saved, so the file always shows
    the charge that was used. A value in the file is never replaced, not even by a config integer.
    """
    path = Path(dp) / GUEST_INPUT
    values = read_guest_input(dp)
    charge = values["net_charge"]
    if charge is None:
        charge = default_charge(values["smiles"])
        if charge is None:
            raise ValueError(
                f"set net_charge (an integer) in {path}: this guest has no SMILES to take it from "
                f"(or set parametrize.net_charge in config.toml to an integer for every empty one)")
        update_guest_input(dp, net_charge=charge)
    if isinstance(charge, bool) or not isinstance(charge, int):
        raise ValueError(f"net_charge in {path} must be an integer, got {charge!r}")
    return charge


def save_anchors(dp, b1_idx, b2_idx):
    """Store the G1/G2 anchors (0-based RDKit guest atom indices) as guest PDB serial numbers."""
    return update_guest_input(dp, anchors=[int(b1_idx) + 1, int(b2_idx) + 1])


def anchors(dp):
    """The G1/G2 anchors as 0-based guest atom indices, or None when not stored yet."""
    serials = read_guest_input(dp)["anchors"]
    if serials is None:
        return None
    if (not isinstance(serials, list) or len(serials) != 2
            or not all(isinstance(s, int) and not isinstance(s, bool) and s >= 1 for s in serials)
            or serials[0] == serials[1]):
        raise ValueError(
            f"anchors in {Path(dp) / GUEST_INPUT} must be two different atom serial numbers "
            f"(column 2 of {GUEST_PDB}), e.g. [1, 5]; got {serials!r}")
    return serials[0] - 1, serials[1] - 1
