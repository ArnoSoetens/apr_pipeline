#!/usr/bin/env python3
"""Build results.csv from the results.json of every data point.

Run it after the analysis array (submit_arrays.py does this for you), or by hand any time.
"""
import csv
import json
import os
import sys

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "..", "apr_pipeline"))

from apr_config import CFG

data_dir = CFG["directories"]["data"]

rows = []
missing = 0
for name in sorted(os.listdir(data_dir)):
    data_point_dir = os.path.join(data_dir, name)
    smi_file = os.path.join(data_point_dir, f"{name}.smi")
    if not os.path.isfile(smi_file):              # not a data point folder
        continue

    result_file = os.path.join(data_point_dir, "results.json")
    if not os.path.isfile(result_file):
        missing += 1
        continue

    with open(smi_file) as f:
        smiles = f.read().strip()
    with open(result_file) as f:
        result = json.load(f)
    rows.append([name, smiles, result["dG_bind"], result["dG_sem"]])

with open(CFG["analysis"]["results_csv"], "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["smiles_folder", "smiles", "dG_bind_kcal_per_mol", "sem"])
    writer.writerows(rows)

print(f"wrote {len(rows)} results to {CFG['analysis']['results_csv']} "
      f"({missing} data points have no results.json yet)")
