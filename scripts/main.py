"""Set up data points: guest structure, placement in the host, parametrization, APR windows.

    python main.py                          every data point folder, each from the stage it still needs
                                            (folders are created from the configured SMILES list only
                                            when there are none)
    python main.py iso_00001 ...            just these data points (this is what the slurm array tasks do)
    python main.py iso_00001-iso_00010      a range of data points
    python main.py --list todo.txt          the data points in a file (one or more per line, # comments)
    python main.py iso_00001/a000-a005      window entries (as run.py takes them) set up their data point
    python main.py --from placement [...]   redo placement and everything after it, whatever exists
    python main.py --force [--from STAGE] NAME ...
                                            redo setup of points that have simulation output: first
                                            delete that output (asks first; --yes without a terminal)

Stages, in order: structure -> placement -> parametrize -> apr_setup. Each data point starts at the
first stage whose output files are missing (datapoints.next_stage), so a folder with only a PDB starts
at placement, a placed guest at parametrize, and so on; existing outputs are trusted.
A data point that fails gets a setup_failed.txt with the stage and reason, and the loop goes on;
one that succeeds loses it. A data point with simulation output (any windows/*/production.nc) is
always skipped, even with --from or by name: setup would overwrite its windows. Only --force
changes that: it deletes, per point and just before its setup, what the new setup makes stale
(datapoints.force_plan: windows/, complex/, results.json, plots/, ...; from parametrize or earlier
also windows_implicit/), then sets the point up from --from or from the start (structure with a
.smi, otherwise placement). Inputs (.smi, input PDBs, guest_input.json) are never deleted.
The stages can also be run one at a time: python apr_pipeline/<stage>.py [NAME ...].
"""
import argparse
import sys
from pathlib import Path

script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir.parent / "apr_pipeline"))   # apr_pipeline/ sits next to scripts/

from apr_config import CFG
from datapoints import (DONE, IMPLICIT_DIR, RESULTS, STAGES, TRAJECTORY, ensure_data_directories,
                        force_plan, next_stage, start_stage)
from selection import add_selection_args, selected

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
add_selection_args(parser)
parser.add_argument("--from", dest="from_stage", metavar="STAGE",
                    choices=STAGES,
                    help="start every selected data point at this stage, even if its outputs exist "
                         f"({', '.join(STAGES)})")
parser.add_argument("--force", action="store_true",
                    help="delete the simulation output of the selected points and redo their setup "
                         "(from --from, or from the start); asks for confirmation first")
parser.add_argument("--yes", action="store_true", help="with --force: do not ask (needed without a terminal)")
args = parser.parse_args()
if args.yes and not args.force:
    parser.error("--yes only goes with --force")

from setup_runner import confirm, disk_usage, remove, run_data_points


def forced_start(dp):
    return args.from_stage or start_stage(dp)


def plan(dp):
    """The stages to run for dp, or why it is skipped."""
    if args.force:
        return STAGES[STAGES.index(forced_start(dp)):]
    start = args.from_stage or next_stage(dp)
    if start is None:
        return f"no input any stage recognizes ({dp.name}.smi or a guest PDB)"
    if start == DONE:
        return "setup is complete (--from apr_setup redoes it)"
    return STAGES[STAGES.index(start):]


def every_data_point():
    """Every folder a stage recognizes; created from the SMILES list when there are none."""
    return [dp.name for dp in ensure_data_directories()]


try:
    names = selected(args, todo=every_data_point, strip_windows=True)   # iso_00001/a000 -> iso_00001
except (FileNotFoundError, ValueError) as e:
    sys.exit(f"main.py: {e}")

prepare = None
if args.force:
    data_dir = Path(CFG["directories"]["data"])
    plans = {name: force_plan(data_dir / name, forced_start(data_dir / name))
             for name in names if (data_dir / name).is_dir()}
    plans = {name: paths for name, paths in plans.items() if paths}
    if plans:
        trajectories = sum(len(list(path.glob(f"*/{TRAJECTORY}"))) for paths in plans.values()
                           for path in paths if path.name == "windows")
        n_results = sum(any(path.name == RESULTS for path in paths) for paths in plans.values())
        n_implicit = sum(any(path.name == IMPLICIT_DIR for path in paths) for paths in plans.values())
        size = sum(disk_usage(paths) for paths in plans.values()) / 1024 ** 3
        shown = ", ".join(list(plans)[:10]) + (", ..." if len(plans) > 10 else "")
        print(f"--force deletes the simulation output of {len(plans)} data point(s): {trajectories} "
              f"trajectories ({size:.1f} GB), results of {n_results}, guest-only runs of {n_implicit}:\n"
              f"  {shown}\nCancel running or queued Slurm jobs for these points first.", flush=True)
        if not confirm("Continue?", assume_yes=args.yes):
            sys.exit("main.py: nothing deleted")

    def prepare(dp):
        paths = plans.get(dp.name)
        if paths:
            remove(paths)
            print("  force: removed " + ", ".join(p.name + ("/" if not p.suffix else "") for p in paths),
                  flush=True)

run_data_points(names, plan, prepare=prepare)

if args.force and plans:          # the deleted results leave results.csv as well
    sys.path.insert(0, str(script_dir))
    import collect_results
    collect_results.main()
