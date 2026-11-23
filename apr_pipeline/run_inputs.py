"""Write the Amber .in files and run.env into <data>/inputs/, filled from config.toml.

    <data>/inputs/minimize_nopbc.in, production_nopbc1.in, production_nopbc2.in
                                 APR windows (scripts/run.sh) and the guest-only p000 window
    <data>/inputs/implicit_only/minimize.in, equilibrate.in, production.in
                                 the free guest (scripts/run_implicit.py --implicit-only)
    <data>/inputs/run.env        the settings run.sh reads (plain KEY="value" lines, no Python needed)

scripts/run.py and scripts/run_implicit.py call write_inputs() before they start pmemd.
"""
import os
from pathlib import Path
from string import Template

from apr_config import CFG, get_pmemd_command

IN_FILES = ("minimize_nopbc.in", "production_nopbc1.in", "production_nopbc2.in")
RUN_ENV = "run.env"
IMPLICIT_ONLY_DIR = "implicit_only"
# implicit-only .in file -> its template (no dummy atoms, no restraints, production continues
# equilibration with its velocities)
IMPLICIT_ONLY_IN = {"minimize.in": "guest_minimize.in", "equilibrate.in": "guest_equilibrate.in",
                    "production.in": "guest_production.in"}


def in_file_values():
    """The ${placeholder} values of each .in template."""
    run = CFG["run"]
    common = {"igb": run["igb"], "temperature": run["temperature"], "restraint_wt": run["restraint_wt"]}
    production = {"dt": run["dt"], "gamma_ln": run["gamma_ln"],
                  "write_every": run["write_every"], "restart_every": run["restart_every"]}
    return {
        "minimize_nopbc.in": {**common, "maxcyc": run["minimize_maxcyc"], "ncyc": run["minimize_ncyc"]},
        "production_nopbc1.in": {**common, **production, "nstlim": run["steps_equil"]},
        "production_nopbc2.in": {**common, **production, "nstlim": run["steps_production"]},
    }


def implicit_only_values():
    """The ${placeholder} values of the implicit-only templates (one set for all three)."""
    run = CFG["run"]
    values = {key: run[key] for key in ("igb", "temperature", "dt", "gamma_ln", "write_every",
                                        "restart_every", "minimize_maxcyc", "minimize_ncyc",
                                        "steps_equil", "steps_production")}
    return {key: (".true." if v is True else ".false." if v is False else v) for key, v in values.items()}


def run_settings():
    """The values run.sh reads from run.env."""
    run = CFG["run"]
    return {
        "DATA_DIR": CFG["directories"]["data"],
        "PMEMD": get_pmemd_command(),
        "INPUT_DIR": CFG["directories"]["inputs"],
        "TOPOLOGY": f'{CFG["naming"]["complex_prefix"]}-dum',
        "DONE_MARKER": run["done_marker"],
        "RANKS": run["mpi_ranks"] if run["pmemd"] == "mpi" else 1,   # CPUs one window uses
    }


def _bash_quoted(value):
    text = str(value)
    for ch in ("\\", '"', "$", "`"):
        text = text.replace(ch, "\\" + ch)
    return f'"{text}"'


def _write_atomic(path, text):
    """Write via a temporary file and a rename, so a reader never sees half a file."""
    tmp = path.with_name(f"{path.name}.tmp{os.getpid()}")
    tmp.write_text(text)
    os.replace(tmp, path)


def write_inputs():
    """Fill <data>/inputs/ with all .in files and run.env. Returns the run.env path."""
    templates = Path(CFG["directories"]["templates"])
    inputs = Path(CFG["directories"]["inputs"])
    inputs.mkdir(parents=True, exist_ok=True)
    for name, values in in_file_values().items():
        _write_atomic(inputs / name, Template((templates / name).read_text()).substitute(values))
    implicit_only = inputs / IMPLICIT_ONLY_DIR
    implicit_only.mkdir(exist_ok=True)
    for name, template in IMPLICIT_ONLY_IN.items():
        _write_atomic(implicit_only / name,
                      Template((templates / template).read_text()).substitute(implicit_only_values()))
    lines = ["# written by apr_pipeline/run_inputs.py from config.toml; read by scripts/run.sh"]
    lines += [f"{key}={_bash_quoted(value)}" for key, value in run_settings().items()]
    env = inputs / RUN_ENV
    _write_atomic(env, "\n".join(lines) + "\n")
    return env
