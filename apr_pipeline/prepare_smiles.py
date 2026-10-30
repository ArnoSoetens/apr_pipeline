#!/usr/bin/env python3
"""Convert the enumerator XML to plain SMILES and make configured RDKit products."""
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from rdkit import Chem

from apr_config import CFG
from structure import find_carbon_backbone_chain


def extract_isomers(xml_path):
    root = ET.parse(xml_path).getroot()
    records = []
    skipped = 0
    for isomer in root.iter("isomer"):
        isomer_id = isomer.get("isomer_id", "unknown")
        smiles_element = isomer.find("constitutional_SMILES")
        smiles = (smiles_element.text or "").strip() if smiles_element is not None else ""
        if not smiles:
            print(f"Warning: isomer {isomer_id} has no SMILES, skipping.", file=sys.stderr)
            skipped += 1
            continue
        records.append({"isomer_id": isomer_id, "smiles": smiles})
    return records, skipped


def write_smiles(smiles, smiles_path):
    Path(smiles_path).parent.mkdir(parents=True, exist_ok=True)
    with open(smiles_path, "w") as smiles_file:
        smiles_file.writelines(f"{value}\n" for value in smiles)


def terminal_carbons(mol):
    return [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetSymbol() == "C"
        and sum(neighbor.GetSymbol() == "C" for neighbor in atom.GetNeighbors()) <= 1
    ]


def candidate_sites(mol, strategy):
    if strategy == "backbone_ends":
        chain = find_carbon_backbone_chain(mol)
        return sorted({chain[0], chain[-1]}) if chain else []
    if strategy == "all_terminal":
        return terminal_carbons(mol)
    raise ValueError(f"Unsupported functionalize.site_strategy: {strategy}")


def f_to_replace(atom):
    fluorines = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == "F"]
    carbon_neighbors = sum(neighbor.GetSymbol() == "C" for neighbor in atom.GetNeighbors())
    if atom.GetTotalNumHs() == 0 and len(fluorines) == 3 and carbon_neighbors == 1:
        return fluorines[0].GetIdx()
    return None


def attach_group(mol, site, group):
    atom = mol.GetAtomWithIdx(site)
    if atom.GetTotalNumHs() > 0:
        atom_to_remove = None
    else:
        atom_to_remove = f_to_replace(atom)
        if atom_to_remove is None:
            return None

    combined = Chem.CombineMols(mol, group)
    editable = Chem.RWMol(combined)
    editable.AddBond(site, mol.GetNumAtoms(), Chem.BondType.SINGLE)
    if atom_to_remove is not None:
        editable.RemoveAtom(atom_to_remove)
    product = editable.GetMol()
    Chem.SanitizeMol(product)
    return product


def parse_group(name, smiles):
    group = Chem.MolFromSmiles(smiles)
    if group is None:
        raise ValueError(f"Invalid SMILES for functional group {name!r}: {smiles!r}")
    if len(Chem.GetMolFrags(group)) != 1:
        raise ValueError(f"Functional group {name!r} must be a single connected fragment")
    return group


def make_products(records, groups, site_strategy):
    parsed_groups = {name: parse_group(name, smiles) for name, smiles in groups.items()}
    products = {}

    for record in records:
        mol = Chem.MolFromSmiles(record["smiles"])
        if mol is None:
            raise ValueError(
                f"Invalid enumerated SMILES for isomer {record['isomer_id']!r}: "
                f"{record['smiles']!r}"
            )
        for site in candidate_sites(mol, site_strategy):
            for group_name, group in parsed_groups.items():
                product = attach_group(mol, site, group)
                if product is None:
                    continue
                product_smiles = Chem.MolToSmiles(
                    product, canonical=True, isomericSmiles=True
                )
                products.setdefault(product_smiles, None)
    return list(products)


def main():
    directories = CFG["directories"]
    records, skipped = extract_isomers(directories["input_xml"])
    if CFG["functionalize"]["enabled"]:
        smiles = make_products(
            records,
            CFG["functionalize"]["groups"],
            CFG["functionalize"]["site_strategy"],
        )
        if not smiles:
            raise ValueError("Functionalization generated no products; check the input SMILES and configured groups")
    else:
        smiles = [record["smiles"] for record in records]

    write_smiles(smiles, directories["smiles_list"])
    if not CFG["functionalize"]["enabled"]:
        print(
            f"Wrote {len(smiles)} enumerated SMILES to {directories['smiles_list']}"
            + (f" ({skipped} source isomers skipped)" if skipped else "")
        )
        return

    print(
        f"Wrote {len(smiles)} unique functionalized SMILES to {directories['smiles_list']}"
        + (f" ({skipped} source isomers skipped)" if skipped else "")
    )


if __name__ == "__main__":
    main()
