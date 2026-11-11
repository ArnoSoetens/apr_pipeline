#!/usr/bin/env python3
"""Build results.csv from the results.json of every data point.

Run it after the analysis array (submit_arrays.py does this for you), or by hand any time.
"""
import csv
import json
import os
import sys
from pathlib import Path

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "..", "apr_pipeline"))

from apr_config import CFG
from datapoints import has_complete_run


def collect_result_rows(data_dir):
    rows = []
    missing = 0
    for data_point in sorted(Path(data_dir).iterdir()):
        if not data_point.is_dir():
            continue

        result_file = data_point / "results.json"
        if not result_file.is_file():
            if has_complete_run(data_point):
                missing += 1
            continue

        # A saved result is durable completion evidence after cleanup removes trajectories.
        smiles_file = data_point / f"{data_point.name}.smi"
        smiles = smiles_file.read_text().strip() if smiles_file.is_file() else ""
        result = json.loads(result_file.read_text())
        rows.append([
            data_point.name,
            smiles,
            result["dG_bind"],
            result["dG_sem"],
        ])
    return rows, missing


def main():
    rows, missing = collect_result_rows(CFG["directories"]["data"])
    with open(CFG["analysis"]["results_csv"], "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["smiles_folder", "smiles", "dG_bind_kcal_per_mol", "sem"])
        writer.writerows(rows)

    print(f"wrote {len(rows)} results to {CFG['analysis']['results_csv']} "
          f"({missing} complete data points have no results.json yet)")


if __name__ == "__main__":
    main()
