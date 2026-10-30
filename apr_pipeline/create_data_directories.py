#!/usr/bin/env python3
"""Create one data-point directory for each SMILES in the configured text file.

Usage:
    python3 create_data_directories.py [<path_to_smiles_txt> <isomers_output_dir>]
"""
import csv
import os
import sys

from apr_config import CFG


def main():
    if len(sys.argv) == 1:
        smiles_path = CFG["directories"]["smiles_list"]
        isomers_dir = CFG["directories"]["data"]
    elif len(sys.argv) == 3:
        smiles_path = sys.argv[1]
        isomers_dir = sys.argv[2]
    else:
        print(f"Usage: {sys.argv[0]} [<path_to_smiles_txt> <isomers_output_dir>]", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(smiles_path):
        print(f"Error: SMILES file not found: {smiles_path}", file=sys.stderr)
        sys.exit(1)
    os.makedirs(isomers_dir, exist_ok=True)
    rows = []
    product_count = 0
    with open(smiles_path) as smiles_file:
        for line_number, line in enumerate(smiles_file, start=1):
            smiles = line.strip()
            if not smiles:
                continue
            product_count += 1
            dir_name = (
                f"{CFG['naming']['folder_prefix']}"
                f"{product_count:0{CFG['naming']['folder_digits']}d}"
            )
            smiles_dir = os.path.join(isomers_dir, dir_name)
            os.makedirs(smiles_dir, exist_ok=True)

            with open(os.path.join(smiles_dir, f"{dir_name}.smi"), "w") as smiles_output:
                smiles_output.write(smiles + "\n")

            rows.append([dir_name, smiles])

    with open(os.path.join(isomers_dir, "isomers.csv"), "w", newline="") as mapping_file:
        writer = csv.writer(mapping_file)
        writer.writerow(["folder", "smiles"])
        writer.writerows(rows)

    print(f"Wrote {product_count} isomer directories to {isomers_dir}")


if __name__ == "__main__":
    main()
