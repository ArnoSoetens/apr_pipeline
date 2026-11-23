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

from apr_config import CFG


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
export("LOG_DIR", D["logs"])
export("COMPLEX_PREFIX", CFG["naming"]["complex_prefix"])
export("TRAJECTORY", CFG["analysis"]["trajectory"])
