"""Set up data points: guest structure, placement in the host, parametrization, APR windows.

    python main.py                  serial: set up every data point folder whose windows are not complete
                                    (folders are created from the configured SMILES list only when there are none)
    python main.py --force          ... and redo the complete ones too
    python main.py iso_00001 ...    set up just these data points (this is what the slurm array tasks do)

Each data point runs its stages in order (structure -> placement -> parametrize -> apr_setup):
all four from <name>.smi, or only apr_setup when the guest PDB, MOL2, FRCMOD and complex PDB exist.
A data point that fails gets a setup_failed.txt with the stage and reason, and the loop goes on;
one that succeeds loses it. A data point with simulation output (any windows/*/production.nc) is
always skipped, even with --force or by name: setup would overwrite its windows.
"""
import argparse
import subprocess
import sys
from pathlib import Path

script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir.parent / "apr_pipeline"))   # apr_pipeline/ sits next to scripts/

from apr_config import CFG

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("names", nargs="*", metavar="NAME",
                    help="data points to set up (default: all that still need it)")
parser.add_argument("--force", action="store_true",
                    help="serial mode: also redo data points whose windows are complete")
args = parser.parse_args()

import apr_setup
import parametrize
import placement
import structure
from datapoints import (
    COMPLEX_PDB,
    GUEST_FRCMOD,
    GUEST_MOL2,
    GUEST_PDB,
    SETUP_STATES,
    STAGES,
    clear_failed,
    ensure_data_directories,
    has_simulation_output,
    mark_failed,
    scan,
)
from guest_input import ensure_guest_input

RUN_STAGE = {
    "structure": structure.run,
    "placement": placement.run,
    "parametrize": parametrize.run,
    "apr_setup": apr_setup.run,
}
assert set(RUN_STAGE) == set(STAGES)


def stages_to_run(dp):
    """All stages from the SMILES, or only APR setup for a fully prepared guest; None: no input."""
    if all((dp / name).is_file() for name in (GUEST_PDB, GUEST_MOL2, GUEST_FRCMOD, COMPLEX_PDB)):
        return ["apr_setup"]
    if (dp / f"{dp.name}.smi").is_file():
        return list(STAGES)
    return None


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
    stages = stages_to_run(dp)
    if stages is None:
        print(name, f"SKIPPED: needs {name}.smi or all prepared guest files "
                    f"({GUEST_PDB}, {GUEST_MOL2}, {GUEST_FRCMOD}, {COMPLEX_PDB})", flush=True)
        n_skipped += 1
        continue
    if not args.names and not args.force and scan(dp)[0] not in SETUP_STATES:
        print(name, "SKIPPED: windows are complete (--force redoes them)", flush=True)
        n_skipped += 1
        continue

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
