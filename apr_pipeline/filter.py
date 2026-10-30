from pathlib import Path
import xml.etree.ElementTree as ET
import random

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors
from rdkit.SimDivFilters import rdSimDivPickers

# ---------- Configuration ----------
INPUT_XML = Path("1000.xml")
OUTPUT_CSV = Path("diverse_1000.csv")
AUDIT_CSV = Path("smiles_validation_audit.csv")

TARGET_N = 2000
SEED = 42
FP_RADIUS = 2
FP_BITS = 2048

# ---------- 1. Parse XML and validate structures ----------
root = ET.parse(INPUT_XML).getroot()

records = []
audit = []

for group in root.findall(".//isomer_list"):
    declared_formula = group.get("formula")

    for item in group.findall("./isomer"):
        mol_id = item.get("isomer_id")
        smiles_node = item.find("constitutional_SMILES")
        smiles = (
            smiles_node.text.strip()
            if smiles_node is not None and smiles_node.text
            else ""
        )

        mol = Chem.MolFromSmiles(smiles) if smiles else None

        if mol is None:
            audit.append({
                "isomer_id": mol_id,
                "original_smiles": smiles,
                "status": "invalid_smiles",
            })
            continue

        canonical = Chem.MolToSmiles(
            mol, canonical=True, isomericSmiles=True
        )
        calculated_formula = rdMolDescriptors.CalcMolFormula(mol)

        if calculated_formula != declared_formula:
            audit.append({
                "isomer_id": mol_id,
                "original_smiles": smiles,
                "canonical_smiles": canonical,
                "declared_formula": declared_formula,
                "calculated_formula": calculated_formula,
                "status": "formula_mismatch",
            })
            # Keep it for now, but flag the discrepancy.
            # Review this record before trusting the selection.

        records.append({
            "isomer_id": mol_id,
            "original_smiles": smiles,
            "canonical_smiles": canonical,
            "formula": declared_formula,
            "calculated_formula": calculated_formula,
            "mol": mol,
        })

df = pd.DataFrame(records)

if df.empty:
    raise ValueError("No valid SMILES were parsed.")

# Deduplicate by canonical molecular structure.
# Keep the first original record and audit later occurrences.
duplicate_mask = df.duplicated("canonical_smiles", keep="first")

for _, row in df.loc[duplicate_mask].iterrows():
    audit.append({
        "isomer_id": row["isomer_id"],
        "original_smiles": row["original_smiles"],
        "canonical_smiles": row["canonical_smiles"],
        "status": "duplicate_canonical_structure",
    })

df = df.loc[~duplicate_mask].reset_index(drop=True)

# Stop rather than silently return fewer than 1000 molecules.
if len(df) < TARGET_N:
    raise ValueError(
        f"Only {len(df)} unique valid structures remain; "
        f"cannot select {TARGET_N}."
    )

# Do not proceed silently if the XML formula labels disagree with RDKit.
mismatches = df[
    df["formula"] != df["calculated_formula"]
]
if not mismatches.empty:
    raise ValueError(
        f"{len(mismatches)} formula mismatches found. "
        "Inspect the source labels and audit before continuing."
    )

print("Valid unique structures:", len(df))
print("Formula groups:", df["formula"].nunique())

# ---------- 2. Generate Morgan fingerprints ----------
generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=FP_RADIUS,
    fpSize=FP_BITS,
)

molecules = df["mol"].tolist()
fingerprints = [generator.GetFingerprint(m) for m in molecules]

# ---------- 3. Choose one reproducible seed per formula ----------
# This guarantees coverage of each formula group.
rng = random.Random(SEED)
first_picks = []

for formula, indices in df.groupby("formula").groups.items():
    first_picks.append(rng.choice(list(indices)))

if len(first_picks) > TARGET_N:
    raise ValueError(
        "There are more formula groups than available selection slots."
    )

# ---------- 4. MaxMin selection ----------
# Tanimoto distance = 1 - Tanimoto similarity.
# MaxMin repeatedly selects structures far from the selected set.
def tanimoto_distance(i, j):
    similarity = DataStructs.TanimotoSimilarity(
        fingerprints[i], fingerprints[j]
    )
    return 1.0 - similarity

picker = rdSimDivPickers.MaxMinPicker()

selected_indices = list(
    picker.LazyBitVectorPick(
        fingerprints,
        len(fingerprints),
        TARGET_N,
        firstPicks=first_picks,
        seed=SEED,
    )
)

if len(selected_indices) != TARGET_N:
    raise RuntimeError(
        f"Expected {TARGET_N} picks; got {len(selected_indices)}."
    )

selected = df.iloc[selected_indices].copy()

# ---------- 5. Validate selection and save ----------
assert len(selected) == TARGET_N
assert selected["canonical_smiles"].nunique() == TARGET_N
assert selected["formula"].nunique() == df["formula"].nunique()

# Do not export RDKit Mol objects to CSV.
selected.drop(columns=["mol"]).to_csv(OUTPUT_CSV, index=False)

pd.DataFrame(audit).to_csv(AUDIT_CSV, index=False)

# Compare formula representation before and after selection.
before = df["formula"].value_counts().rename("before")
after = selected["formula"].value_counts().rename("selected")
summary = pd.concat([before, after], axis=1).fillna(0).astype(int)
summary["selected_fraction"] = summary["selected"] / summary["before"]
summary = summary.sort_index()

print(f"\nSaved {len(selected)} molecules to {OUTPUT_CSV}")
print(f"Saved validation audit to {AUDIT_CSV}")
print("\nFormula representation:")
print(summary.to_string())