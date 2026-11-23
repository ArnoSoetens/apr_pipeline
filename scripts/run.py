"""Run the APR windows of set-up data points with pmemd (minimize -> production part 1 -> part 2).

    python run.py                          every unfinished window of every set-up data point
    python run.py iso_00001 ...            the unfinished windows of these data points
    python run.py iso_00001/a000           one window;  iso_00001/a000-a005  a range of windows
    python run.py iso_00001-iso_00010      a range of data points;  --list todo.txt  entries from a file
    python run.py --cpus 4 ...             4 CPUs (default: run.max_jobs; 0 = all CPUs); --cpus 1 = serial
    python run.py --dry-run ...            only list the windows that would run

A finished window (done marker + production.nc) is never run again. A data point whose setup is not
finished is skipped. The .in files and run.env are written into <data>/inputs/ from config.toml,
then scripts/run.sh runs the windows (see there; a failed window gets run_failed.txt).
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir.parent / "apr_pipeline"))   # apr_pipeline/ sits next to scripts/

from apr_config import CFG
from datapoints import DONE, WINDOWS, all_data_points, next_stage, window_finished
from selection import add_selection_args, selected

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
add_selection_args(parser, windows=True)
parser.add_argument("--cpus", type=int, metavar="N",
                    help="CPUs to use (default: run.max_jobs in config.toml; 0 = all CPUs)")
parser.add_argument("--dry-run", action="store_true", help="only list the windows that would run")
args = parser.parse_args()
if args.cpus is not None and args.cpus < 1:
    parser.error("--cpus must be at least 1")

data_dir = Path(CFG["directories"]["data"])


def every_set_up_point():
    return [(dp.name, None) for dp in all_data_points() if next_stage(dp) == DONE]


try:
    entries = selected(args, todo=every_set_up_point, windows=True)
except (FileNotFoundError, ValueError) as e:
    sys.exit(f"run.py: {e}")
unknown = sorted({w for _, w in entries if w is not None and w not in WINDOWS})
if unknown:
    sys.exit(f"run.py: not a configured window: {', '.join(unknown)} (windows are {WINDOWS[0]} ... {WINDOWS[-1]})")

named = bool(args.names or args.list_file)   # chosen by hand (not the default selection)
to_run = {}                                  # data point -> its window folders to run, in order
skipped = set()


def skip(name, reason):
    if name not in skipped:                  # once per data point, not once per window
        print(f"{name} SKIPPED: {reason}", flush=True)
        skipped.add(name)


for name, window in entries:
    dp = data_dir / name
    if not dp.is_dir():
        skip(name, "no such data point folder")
        continue
    stage = next_stage(dp)
    if stage != DONE:
        skip(name, f"setup not finished (next stage: {stage})")
        continue
    folders = [dp / "windows" / w for w in ([window] if window else WINDOWS)]
    done = [f for f in folders if window_finished(f)]
    if window and done:
        print(f"{name}/{window} SKIPPED: already done", flush=True)
    elif named and not window and len(done) == len(folders):
        print(f"{name} SKIPPED: all windows already done", flush=True)
    queue = to_run.setdefault(name, [])
    queue += [f for f in folders if f not in done and f not in queue]

window_dirs = [f for folders in to_run.values() for f in folders]
if not window_dirs:
    print("nothing to run", flush=True)
    sys.exit(0)
print(f"{len(window_dirs)} window(s) in {sum(1 for f in to_run.values() if f)} data point(s)", flush=True)
if args.dry_run:
    for name, folders in to_run.items():
        if folders:
            print(f"  {name}: {len(folders)} window(s): {' '.join(f.name for f in folders)}")
    sys.exit(0)

cpus = args.cpus or CFG["run"]["max_jobs"] or int(os.environ.get("SLURM_CPUS_PER_TASK") or os.cpu_count() or 1)

from run_inputs import write_inputs

env_file = write_inputs()
result = subprocess.run(["bash", str(script_dir / "run.sh"), "--env", str(env_file), "--cpus", str(cpus),
                         *map(str, window_dirs)])
sys.exit(result.returncode)
