#!/usr/bin/env python3
"""Print `export` lines for bash, built from <data>/inputs/config.toml.

In a bash script:
    CONFIG_ENV="$(python config_env.py)"
    eval "$CONFIG_ENV"

The config decides: a variable with the same name in the environment is overwritten.
"""
import os
import sys
from pathlib import Path

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "..", "apr_pipeline"))

from apr_config import CFG, get_pmemd_command


def dq(value):
    """Escape a value so it is safe inside bash double quotes."""
    s = str(value)
    for ch in ("\\", '"', "$", "`"):
        s = s.replace(ch, "\\" + ch)
    return s


def export(name, value, raw=False):
    text = value if raw else dq(value)
    print(f'export {name}="{text}"')


D = CFG["directories"]
export("DATA_DIR", D["data"])
export("HOST_PDB", str(Path(D["host"]) / f'{CFG["naming"]["host_prefix"]}.pdb'))
export("HOST_PREFIX", CFG["naming"]["host_prefix"])
export("GUEST_PREFIX", CFG["naming"]["guest_prefix"])
export("INPUT_DIR", D["inputs"])
export("LOG_DIR", D["logs"])
export("COMPLEX_PREFIX", CFG["naming"]["complex_prefix"])
export("TRAJECTORY", CFG["analysis"]["trajectory"])
export("DONE_MARKER", CFG["run"]["done_marker"])
export("PMEMD", get_pmemd_command())

if CFG["run"]["max_jobs"] > 0:
    export("MAX_JOBS", CFG["run"]["max_jobs"])
else:   # 0 = as many windows as CPUs: the slurm request, or all cores outside slurm
    export("MAX_JOBS", "${SLURM_CPUS_PER_TASK:-$(nproc)}", raw=True)
