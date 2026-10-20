#!/usr/bin/env python3
"""Fill the Amber .in templates with the values from config.toml.

templates/*.in  (with ${placeholders})  ->  inputs/*.in
run.sh calls this before starting pmemd; you can also run it by hand to look at the result.
"""
import os
import sys
from pathlib import Path
from string import Template

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "..", "apr_pipeline"))

from apr_config import CFG

run = CFG["run"]
templates_dir = Path(CFG["directories"]["templates"])
inputs_dir = Path(CFG["directories"]["inputs"])

# values shared by all files
common = {
    "igb": run["igb"],
    "temperature": run["temperature"],
    "restraint_wt": run["restraint_wt"],
}
# values for the production files
production = {
    "dt": run["dt"],
    "gamma_ln": run["gamma_ln"],
    "write_every": run["write_every"],
    "restart_every": run["restart_every"],
}

files = {
    "minimize_nopbc.in": {"maxcyc": run["minimize_maxcyc"], "ncyc": run["minimize_ncyc"]},
    "production_nopbc1.in": {**production, "nstlim": run["steps_equil"]},
    "production_nopbc2.in": {**production, "nstlim": run["steps_production"]},
}

inputs_dir.mkdir(parents=True, exist_ok=True)
for name, extra in files.items():
    text = Template((templates_dir / name).read_text()).substitute({**common, **extra})
    tmp = inputs_dir / f"{name}.tmp{os.getpid()}"
    tmp.write_text(text)
    os.replace(tmp, inputs_dir / name)    # atomic: many array tasks may render at the same moment
    print(f"wrote {inputs_dir / name}")
