"""Write what scripts/run.sh needs into <data>/inputs/: the Amber .in files and run.env.

    scripts/templates/*.in  (with ${placeholders})  ->  <data>/inputs/*.in   filled from config.toml
    <data>/inputs/run.env    the settings run.sh reads (plain KEY="value" lines, no Python needed)

scripts/run.py calls write_inputs() before it starts run.sh.
"""
import os
from pathlib import Path
from string import Template

from apr_config import CFG, get_pmemd_command

IN_FILES = ("minimize_nopbc.in", "production_nopbc1.in", "production_nopbc2.in")
RUN_ENV = "run.env"


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


def run_settings():
    """The values run.sh reads from run.env."""
    run = CFG["run"]
    return {
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
    """Fill <data>/inputs/ with the three .in files and run.env. Returns the run.env path."""
    templates = Path(CFG["directories"]["templates"])
    inputs = Path(CFG["directories"]["inputs"])
    inputs.mkdir(parents=True, exist_ok=True)
    for name, values in in_file_values().items():
        _write_atomic(inputs / name, Template((templates / name).read_text()).substitute(values))
    lines = ["# written by apr_pipeline/run_inputs.py from config.toml; read by scripts/run.sh"]
    lines += [f"{key}={_bash_quoted(value)}" for key, value in run_settings().items()]
    env = inputs / RUN_ENV
    _write_atomic(env, "\n".join(lines) + "\n")
    return env
