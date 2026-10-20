import csv
import os
import sys

script_dir = os.path.dirname(os.path.abspath(__file__))

# bootstrap: just enough to import apr_config (apr_pipeline/ sits next to scripts/)
sys.path.insert(0, os.path.abspath(os.path.join(script_dir, "..", "apr_pipeline")))

from apr_config import CFG

data_dir = CFG["directories"]["data"]
sys.path.insert(0, CFG["directories"]["pipeline"])

from apr_analysis import run

# Data points to analyze: the names given as arguments (e.g. by a slurm array task),
# otherwise every folder in data/. With arguments results.csv is not written here:
# collect_results.py builds it afterwards from all results.json files.
names = sys.argv[1:] or sorted(os.listdir(data_dir))

summary = []
for name in names:
    data_point_dir = os.path.join(data_dir, name)

    if not os.path.isdir(os.path.join(data_point_dir, "windows")):
        continue

    # the real SMILES, read from the .smi file in this folder
    smi_file = os.path.join(data_point_dir, f"{name}.smi")
    if os.path.exists(smi_file):
        with open(smi_file) as f:
            smiles = f.read().strip()
    else:
        smiles = ""          # no .smi file found, leave the column empty

    try:
        out = run(data_point_dir)
    except Exception as e:
        print(f"analysis {name} FAILED: {e}")
        continue

    summary.append([name, smiles, out["dG_bind"], out["dG_sem"]])
    print(f"analysis {name} done: {out['dG_bind']:.2f} +/- {out['dG_sem']:.2f} kcal/mol")

if not sys.argv[1:]:
    with open(CFG["analysis"]["results_csv"], "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["smiles_folder", "smiles", "dG_bind_kcal_per_mol", "sem"])
        writer.writerows(summary)