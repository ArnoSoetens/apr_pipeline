"""Loads config.toml once and exposes it as CFG (a plain dict) plus a few derived values.

Every module does:   from apr_config import CFG
The config file is  <data>/inputs/config.toml. The data folder is <repository>/data, or the folder in
the environment variable APR_DATA. CFG["directories"] is not read from the config but derived:
    data       the data folder (data point folders: iso_00001, ...)
    inputs     <data>/inputs: config.toml, host files, SMILES list, Amber .in files, results.csv
    host       <data>/inputs (BCD.pdb, BCD.mol2, BCD.frcmod)
    logs       <data>/logs: slurm output and worklists
    smiles_list <data>/inputs/<input.smiles_list>
    scripts, pipeline, templates: the code's own folders
"""
import math
import os
import sys
from pathlib import Path

try:
    import tomllib                    # python >= 3.11
except ModuleNotFoundError:           # older python: pip install tomli
    import tomli as tomllib

import numpy as np

PIPELINE_DIR = Path(__file__).resolve().parent
REPO_DIR = PIPELINE_DIR.parent
DATA_DIR = Path(os.environ.get("APR_DATA") or REPO_DIR / "data").expanduser().resolve()
INPUTS = "inputs"
LOGS = "logs"
RESERVED = (INPUTS, LOGS)              # folders in the data folder that are never a data point
CONFIG_PATH = DATA_DIR / INPUTS / "config.toml"

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

    wall = a["restraints"]["walls"]["distance"]
    need(
        isinstance(wall["enabled"], bool),
        "apr.restraints.walls.distance.enabled must be true or false",
    )
    need(
        wall["mode"] in ("manual", "rim"),
        "apr.restraints.walls.distance.mode must be 'manual' or 'rim'",
    )
    if wall["mode"] == "rim":
        rim_groups = wall.get("rim_atom_groups")
        need(
            isinstance(rim_groups, list)
            and len(rim_groups) == 2
            and all(isinstance(group, list) and len(group) == 7 for group in rim_groups),
            "apr.restraints.walls.distance.rim_atom_groups must contain two groups of 7 host atom numbers in rim mode",
        )
        rim_atom_numbers = [number for group in rim_groups for number in group]
        need(
            all(isinstance(number, int) and not isinstance(number, bool) and number > 0
                for number in rim_atom_numbers)
            and len(set(rim_atom_numbers)) == len(rim_atom_numbers),
            "apr.restraints.walls.distance.rim_atom_groups must contain unique positive integers",
        )

    bounds = [wall[key] for key in ("r1", "r2", "r3", "r4")]
    need(
        all(isinstance(value, (int, float))
            and not isinstance(value, bool)
            and math.isfinite(value)
            for value in bounds),
        "apr.restraints.walls.distance.r1-r4 must be finite numbers",
    )
    if wall["mode"] == "manual":
        need(
            bounds[0] < bounds[1] < bounds[2] < bounds[3],
            "manual apr.restraints.walls.distance bounds must satisfy r1 < r2 < r3 < r4",
        )
    for key in ("rk2", "rk3"):
        need(
            isinstance(wall[key], (int, float))
            and not isinstance(wall[key], bool)
            and math.isfinite(wall[key])
            and wall[key] > 0,
            f"apr.restraints.walls.distance.{key} must be a finite positive number",
        )
    need(
        wall["guest"] in ("g1", "centroid"),
        "apr.restraints.walls.distance.guest must be 'g1' or 'centroid'",
    )
    if wall["mode"] == "rim" and wall["guest"] == "g1":
        rim_offset = wall.get("rim_offset")
        need(
            isinstance(rim_offset, (int, float))
            and not isinstance(rim_offset, bool)
            and math.isfinite(rim_offset)
            and rim_offset >= 0,
            "apr.restraints.walls.distance.rim_offset must be a finite "
            "nonnegative distance in rim mode",
        )
    if wall["mode"] == "rim" and wall["guest"] == "centroid":
        radius_offset = wall.get("radius_offset")
        need(
            isinstance(radius_offset, (int, float))
            and not isinstance(radius_offset, bool)
            and math.isfinite(radius_offset),
            "apr.restraints.walls.distance.radius_offset must be a finite distance "
            "(Angstrom, may be negative) in rim mode with guest = 'centroid'",
        )
        parabolic_width = wall.get("parabolic_width")
        need(
            isinstance(parabolic_width, (int, float))
            and not isinstance(parabolic_width, bool)
            and math.isfinite(parabolic_width)
            and parabolic_width > 0,
            "apr.restraints.walls.distance.parabolic_width must be a finite "
            "positive distance in rim mode",
        )

    def finite_number(value):
        return (isinstance(value, (int, float)) and not isinstance(value, bool)
                and math.isfinite(value))

    cone = a["restraints"]["walls"]["cone"]
    need(isinstance(cone["enabled"], bool),
         "apr.restraints.walls.cone.enabled must be true or false")
    need(cone["guest"] in ("g1", "centroid"),
         "apr.restraints.walls.cone.guest must be 'g1' or 'centroid'")
    need(finite_number(cone["angle_offset"]),
         "apr.restraints.walls.cone.angle_offset must be a finite number of degrees")
    for key in ("parabolic_width", "rk2"):
        need(finite_number(cone[key]) and cone[key] > 0,
             f"apr.restraints.walls.cone.{key} must be a finite positive number")
    if cone["enabled"]:
        rim_groups = wall.get("rim_atom_groups")
        need(
            isinstance(rim_groups, list)
            and len(rim_groups) == 2
            and all(isinstance(group, list) and len(group) == 7
                    and all(isinstance(n, int) and not isinstance(n, bool) and n > 0 for n in group)
                    for group in rim_groups),
            "apr.restraints.walls.cone needs apr.restraints.walls.distance.rim_atom_groups "
            "(two groups of 7 host oxygen numbers), also when the distance wall is manual or off",
        )

    s = cfg["slurm"]
    need(s["max_array_size"] <= 1001, "slurm.max_array_size is above the cluster maximum (1001)")
    need(s["sim_cpus"] >= 1, "slurm.sim_cpus must be >= 1")
    need(s["sim_task_hours"] > 0 and s["sim_minutes_per_window"] > 0, "slurm.sim_task_hours and sim_minutes_per_window must be > 0")
    need(s["analysis_data_points_per_task"] >= 1, "slurm.analysis_data_points_per_task must be >= 1")
    need(s["setup_data_points_per_task"] >= 1 and s["setup_cpus"] >= 1, "slurm.setup_data_points_per_task and setup_cpus must be >= 1")
    need(an["method"] in METHODS, f"analysis.method must be one of {list(METHODS)}")
    need(isinstance(an["enhanced"], bool), "analysis.enhanced must be true or false")
    need(r["pmemd"] in PMEMD_KINDS, f"run.pmemd must be one of {list(PMEMD_KINDS)}")
    functionalize = cfg["functionalize"]
    need(
        isinstance(functionalize["enabled"], bool),
        "functionalize.enabled must be true or false",
    )
    need(
        functionalize["site_strategy"] in ("backbone_ends", "all_terminal"),
        "functionalize.site_strategy must be 'backbone_ends' or 'all_terminal'",
    )
    need(
        functionalize["replace_policy"] == "h_first_cf3_only",
        "functionalize.replace_policy must be 'h_first_cf3_only'",
    )
    need(bool(functionalize["groups"]), "functionalize.groups must contain at least one group")
    net_charge = cfg["parametrize"]["net_charge"]
    need(
        net_charge == "auto"
        or (isinstance(net_charge, int) and not isinstance(net_charge, bool)),
        "parametrize.net_charge must be 'auto' or an integer",
    )
    output_pdb = cfg["structure"]["output_pdb"]
    need(
        isinstance(output_pdb, str) and output_pdb.endswith(".pdb")
        and Path(output_pdb).name == output_pdb
        and output_pdb not in (f'{cfg["naming"]["complex_prefix"]}.pdb', "sqm.pdb"),
        "structure.output_pdb must be a file name ending in .pdb (not the complex PDB or sqm.pdb)",
    )
    _radii_for(cfg)


def directories(data_dir, smiles_list):
    """The folders, derived from the data folder and the code's location."""
    data_dir = Path(data_dir)
    inputs = data_dir / INPUTS
    return {
        "data": str(data_dir),
        "inputs": str(inputs),
        "host": str(inputs),
        "logs": str(data_dir / LOGS),
        "smiles_list": str(inputs / smiles_list),
        "scripts": str(REPO_DIR / "scripts"),
        "pipeline": str(PIPELINE_DIR),
        "templates": str(REPO_DIR / "scripts" / "templates"),
    }


def load_config(path=CONFIG_PATH):
    """Read <data>/inputs/config.toml; the data folder is the folder above inputs/."""
    path = Path(path).resolve()
    if not path.is_file():
        raise FileNotFoundError(
            f"no settings found: expected {path} (the data folder is {path.parent.parent}; "
            f"set APR_DATA to use another data folder)")
    with open(path, "rb") as f:
        cfg = tomllib.load(f)

    if "directories" in cfg:
        print(f"note: [directories] in {path} is ignored; the folders follow from the data folder",
              file=sys.stderr)
    cfg["directories"] = directories(path.parent.parent,
                                     cfg.get("input", {}).get("smiles_list", "smiles.txt"))

    # analysis.convergence_plot: part of analysis.enhanced now; ignored if present
    cfg["analysis"].pop("convergence_plot", None)
    cfg["analysis"].setdefault("enhanced", False)
    for key, value in (("implicit_cpus", 1), ("implicit_data_points_per_task", 1),
                       ("implicit_minutes_per_data_point", 10), ("implicit_mem", "1G")):
        cfg["slurm"].setdefault(key, value)
    cfg["analysis"].setdefault("results_csv", "inputs/results.csv")
    cfg["structure"].setdefault("output_pdb", f'{cfg["naming"]["guest_prefix"]}.pdb')
    walls = cfg["apr"]["restraints"].setdefault("walls", {})
    distance_wall = walls.setdefault("distance", {})
    distance_wall.setdefault("enabled", False)
    distance_wall.setdefault("mode", "manual")
    distance_wall.setdefault("r1", 0.0)
    distance_wall.setdefault("r2", 1.0)
    distance_wall.setdefault("r3", 12.0)
    distance_wall.setdefault("r4", 13.0)
    distance_wall.setdefault("rk2", 100.0)
    distance_wall.setdefault("rk3", 100.0)
    distance_wall.setdefault("guest", "g1")
    if distance_wall["mode"] == "rim":
        distance_wall.setdefault("rim_offset", 1.0)
        distance_wall.setdefault("radius_offset", -1.0)
        distance_wall.setdefault("parabolic_width", 1.0)
    cone_wall = walls.setdefault("cone", {})
    cone_wall.setdefault("enabled", False)
    cone_wall.setdefault("guest", "g1")
    cone_wall.setdefault("angle_offset", 0.0)
    cone_wall.setdefault("parabolic_width", 5.0)
    cone_wall.setdefault("rk2", 100.0)
    csv = cfg["analysis"]["results_csv"]              # relative to the data directory, so that
    cfg["analysis"]["results_csv"] = str((Path(cfg["directories"]["data"]) / csv).resolve())  # data/ holds it

    naming = cfg["naming"]
    if not naming.get("complex_prefix"):
        naming["complex_prefix"] = f'{naming["host_prefix"]}-{naming["guest_prefix"]}'

    _validate(cfg)
    return cfg


CFG = load_config()
Path(CFG["directories"]["logs"]).mkdir(exist_ok=True)   # slurm needs it before a job starts (--output)


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
