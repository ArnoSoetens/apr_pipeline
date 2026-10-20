"""Loads config.toml once and exposes it as CFG (a plain dict) plus a few derived values.

Every module does:   from apr_config import CFG
The config file is  <project>/config.toml  unless the environment variable APR_CONFIG points elsewhere.
"""
import os
from pathlib import Path

try:
    import tomllib                    # python >= 3.11
except ModuleNotFoundError:           # older python: pip install tomli
    import tomli as tomllib

import numpy as np

PIPELINE_DIR = Path(__file__).resolve().parent
CONFIG_PATH = Path(os.environ.get("APR_CONFIG", PIPELINE_DIR.parent / "config.toml"))

# tleap radii recommended for each GB model (Amber manual table)
RADII_BY_IGB = {1: "mbondi", 2: "mbondi2", 5: "mbondi2", 7: "bondi", 8: "mbondi3"}
# config name -> pAPRika method name
METHODS = {"ti": "ti-block", "mbar": "mbar-autoc"}
PMEMD_KINDS = ("serial", "cuda", "mpi")


def _radii_for(cfg):
    run = cfg["run"]
    if run["radii_set"] != "auto":
        return run["radii_set"]
    if run["igb"] not in RADII_BY_IGB:
        raise ValueError(f"config.toml: no radii set known for igb = {run['igb']}; "
                         f"set run.radii_set explicitly")
    return RADII_BY_IGB[run["igb"]]


def _validate(cfg):
    def need(ok, msg):
        if not ok:
            raise ValueError(f"config.toml: {msg}")

    p, a, r, an = cfg["placement"], cfg["apr"], cfg["run"], cfg["analysis"]
    need(len(cfg["naming"]["guest_prefix"]) <= 3, "naming.guest_prefix can be at most 3 characters (PDB residue name)")
    need(p["clash_cutoff"] > 0, "placement.clash_cutoff must be > 0")
    need(len(p["cavity_center"]) == 3, "placement.cavity_center needs 3 numbers")
    need(len(p["cavity_axis"]) == 3 and any(p["cavity_axis"]), "placement.cavity_axis needs 3 numbers, not all 0")

    fractions = [float(x) for x in a["attach_string"].split()]
    need(fractions == sorted(fractions), "apr.attach_string must be in ascending order")
    need(a["step"] > 0, "apr.step must be > 0")
    need(a["initial_distance"] < a["final_distance"], "apr.initial_distance must be < apr.final_distance")
    need(len(a["dummy_z"]) == 3, "apr.dummy_z needs 3 numbers")
    need(len(a["host_anchor_atoms"]) == 3, "apr.host_anchor_atoms needs 3 atom numbers")
    need(len(a["restraints"]["static_fc"]) == 6, "apr.restraints.static_fc needs 6 numbers")

    s = cfg["slurm"]
    need(s["max_array_size"] <= 1001, "slurm.max_array_size is above the cluster maximum (1001)")
    need(s["sim_cpus"] >= 1 and s["sim_windows_per_cpu"] >= 1, "slurm.sim_cpus and sim_windows_per_cpu must be >= 1")
    need(s["sim_task_hours"] > 0 and s["sim_minutes_per_window"] > 0, "slurm.sim_task_hours and sim_minutes_per_window must be > 0")
    need(s["analysis_data_points_per_task"] >= 1, "slurm.analysis_data_points_per_task must be >= 1")
    need(an["method"] in METHODS, f"analysis.method must be one of {list(METHODS)}")
    need(r["pmemd"] in PMEMD_KINDS, f"run.pmemd must be one of {list(PMEMD_KINDS)}")
    _radii_for(cfg)


def load_config(path=CONFIG_PATH):
    with open(path, "rb") as f:
        cfg = tomllib.load(f)

    base = Path(path).resolve().parent                 # relative paths start here
    for key, value in cfg["directories"].items():
        cfg["directories"][key] = str((base / value).resolve())
    csv = cfg["analysis"]["results_csv"]
    cfg["analysis"]["results_csv"] = str((base / csv).resolve())

    naming = cfg["naming"]
    if not naming.get("complex_prefix"):
        naming["complex_prefix"] = f'{naming["host_prefix"]}-{naming["guest_prefix"]}'

    _validate(cfg)
    return cfg


CFG = load_config()


# ---- values derived from the config (shared by setup and analysis) ----------

def get_attach_fractions():
    return [float(i) / 100 for i in CFG["apr"]["attach_string"].split()]


def get_pull_distances():
    a = CFG["apr"]
    # "+ step / 2" keeps the final distance in despite float rounding
    return np.arange(a["initial_distance"], a["final_distance"] + a["step"] / 2, a["step"])


def get_radii_set():
    return _radii_for(CFG)


def get_method():
    return METHODS[CFG["analysis"]["method"]]


def get_pmemd_command():
    run = CFG["run"]
    if run["pmemd"] == "cuda":
        return "pmemd.cuda"
    if run["pmemd"] == "mpi":
        return f"mpirun -np {run['mpi_ranks']} pmemd.MPI"
    return "pmemd"
