#!/usr/bin/env python3
"""
Parses an isomer_enumeration XML file and, for every <isomer>, creates:

    <isomers_dir>/iso_00001/iso_00001.smi      (one line: the SMILES)

and writes <isomers_dir>/isomers.csv mapping folder -> isomer_id -> SMILES.

Usage:
    python3 create_data_directories.py <path_to_xml> <isomers_output_dir>
"""
import csv
import sys
import os
import xml.etree.ElementTree as ET


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <path_to_xml> <isomers_output_dir>", file=sys.stderr)
        sys.exit(1)

    xml_path = sys.argv[1]
    isomers_dir = sys.argv[2]

    if not os.path.isfile(xml_path):
        print(f"Error: XML file not found: {xml_path}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(isomers_dir, exist_ok=True)

    root = ET.parse(xml_path).getroot()

    rows = []
    skipped = 0

    for number, isomer in enumerate(root.iter("isomer"), start=1):
        isomer_id = isomer.get("isomer_id", "unknown")
        smiles_el = isomer.find("constitutional_SMILES")

        if smiles_el is None or not (smiles_el.text or "").strip():
            print(f"Warning: isomer {isomer_id} has no SMILES, skipping.", file=sys.stderr)
            skipped += 1
            continue

        smiles = smiles_el.text.strip()
        dir_name = f"iso_{number:05d}"

        smiles_dir = os.path.join(isomers_dir, dir_name)
        os.makedirs(smiles_dir, exist_ok=True)

        with open(os.path.join(smiles_dir, f"{dir_name}.smi"), "w") as f:
            f.write(smiles + "\n")

        rows.append([dir_name, isomer_id, smiles])

    with open(os.path.join(isomers_dir, "isomers.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["folder", "isomer_id", "smiles"])
        writer.writerows(rows)

    print(f"Wrote {len(rows)} isomer directories to {isomers_dir}"
          + (f" ({skipped} skipped)" if skipped else ""))


if __name__ == "__main__":
    main()