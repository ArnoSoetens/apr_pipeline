"""Set up data points: guest structure, placement in the host, parametrization, APR windows.

    python main.py                          every data point folder, each from the stage it still needs
                                            (folders are created from the configured SMILES list only
                                            when there are none)
    python main.py iso_00001 ...            just these data points (this is what the slurm array tasks do)
    python main.py --from placement [...]   redo placement and everything after it, whatever exists

Stages, in order: structure -> placement -> parametrize -> apr_setup. Each data point starts at the
first stage whose output files are missing (datapoints.next_stage), so a folder with only a PDB starts
at placement, a placed guest at parametrize, and so on; existing outputs are trusted.
A data point that fails gets a setup_failed.txt with the stage and reason, and the loop goes on;
one that succeeds loses it. A data point with simulation output (any windows/*/production.nc) is
always skipped, even with --from or by name: setup would overwrite its windows.
"""
import argparse
import subprocess
import sys
from pathlib import Path

script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir.parent / "apr_pipeline"))   # apr_pipeline/ sits next to scripts/

from apr_config import CFG
from datapoints import (
    DONE,
    STAGES,
    clear_failed,
    ensure_data_directories,
    has_simulation_output,
    mark_failed,
    next_stage,
)

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("names", nargs="*", metavar="NAME",
                    help="data points to set up (default: all that still need it)")
parser.add_argument("--from", dest="from_stage", metavar="STAGE",
                    choices=STAGES,
                    help="start every selected data point at this stage, even if its outputs exist "
                         f"({', '.join(STAGES)})")
args = parser.parse_args()

import apr_setup
import parametrize
import placement
import structure
from guest_input import ensure_guest_input

RUN_STAGE = {
    "structure": structure.run,
    "placement": placement.run,
    "parametrize": parametrize.run,
    "apr_setup": apr_setup.run,
}
assert set(RUN_STAGE) == set(STAGES)


def describe(error):
    """One-line reason, including the end of a failed command's output."""
    if isinstance(error, subprocess.CalledProcessError):
        output = (error.stderr or error.stdout or "").strip().splitlines()
        return f"{error}" + (f" -- {output[-1]}" if output else "")
    return str(error) or type(error).__name__


if args.names:
    names = args.names
else:
    try:
        names = [dp.name for dp in ensure_data_directories()]
    except FileNotFoundError as e:
        sys.exit(f"main.py: {e}")

data_dir = Path(CFG["directories"]["data"])
n_done = n_failed = n_skipped = 0
for i, name in enumerate(names, 1):
    dp = data_dir / name
    print(f"[{i}/{len(names)}] {name}", flush=True)

    if has_simulation_output(dp):
        print(name, f'SKIPPED: has simulation output (windows/*/{CFG["analysis"]["trajectory"]}); '
                    "setup never touches it", flush=True)
        n_skipped += 1
        continue
    if not dp.is_dir():
        print(name, f"SKIPPED: no such data point folder: {dp}", flush=True)
        n_skipped += 1
        continue
    start = args.from_stage or next_stage(dp)
    if start is None:
        print(name, f"SKIPPED: no input any stage recognizes ({name}.smi or a guest PDB)", flush=True)
        n_skipped += 1
        continue
    if start == DONE:
        print(name, "SKIPPED: setup is complete (--from apr_setup redoes it)", flush=True)
        n_skipped += 1
        continue
    stages = STAGES[STAGES.index(start):]

    stage = "guest_input"
    try:
        ensure_guest_input(dp)
        for stage in stages:
            print(f"  {stage}", flush=True)
            RUN_STAGE[stage](dp)
    except Exception as e:
        reason = f"{stage}: {describe(e)}"
        print(name, "FAILED in", reason, flush=True)
        mark_failed(dp, reason)
        n_failed += 1
        continue

    clear_failed(dp)
    n_done += 1

print(f"setup finished: {n_done} done, {n_failed} failed, {n_skipped} skipped", flush=True)
