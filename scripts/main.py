"""Set up data points: guest structure, placement in the host, parametrization, APR windows.

    python main.py                          every data point folder, each from the stage it still needs
                                            (folders are created from the configured SMILES list only
                                            when there are none)
    python main.py iso_00001 ...            just these data points (this is what the slurm array tasks do)
    python main.py iso_00001-iso_00010      a range of data points
    python main.py --list todo.txt          the data points in a file (one or more per line, # comments)
    python main.py iso_00001/a000-a005      window entries (as run.py takes them) set up their data point
    python main.py --from placement [...]   redo placement and everything after it, whatever exists

Stages, in order: structure -> placement -> parametrize -> apr_setup. Each data point starts at the
first stage whose output files are missing (datapoints.next_stage), so a folder with only a PDB starts
at placement, a placed guest at parametrize, and so on; existing outputs are trusted.
A data point that fails gets a setup_failed.txt with the stage and reason, and the loop goes on;
one that succeeds loses it. A data point with simulation output (any windows/*/production.nc) is
always skipped, even with --from or by name: setup would overwrite its windows.
The stages can also be run one at a time: python apr_pipeline/<stage>.py [NAME ...].
"""
import argparse
import sys
from pathlib import Path

script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir.parent / "apr_pipeline"))   # apr_pipeline/ sits next to scripts/

from datapoints import DONE, STAGES, ensure_data_directories, next_stage
from selection import add_selection_args, selected

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
add_selection_args(parser)
parser.add_argument("--from", dest="from_stage", metavar="STAGE",
                    choices=STAGES,
                    help="start every selected data point at this stage, even if its outputs exist "
                         f"({', '.join(STAGES)})")
args = parser.parse_args()

from setup_runner import run_data_points


def plan(dp):
    """The stages to run for dp, or why it is skipped."""
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

run_data_points(names, plan)
