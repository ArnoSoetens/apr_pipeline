#!/usr/bin/env python3
"""Submit the simulations and the analysis as slurm arrays.

    python submit_arrays.py status                 # what is finished and what is left
    python submit_arrays.py all                    # show what would be submitted (nothing is submitted)
    python submit_arrays.py all --submit           # simulations -> analysis -> results.csv
    python submit_arrays.py sim --submit           # only the simulations
    python submit_arrays.py analysis --submit      # only the analysis

Simulation
    Every window that still has to be simulated goes into ONE list (lines like iso_00001/a000).
    The list is cut into work packages of  sim_cpus x sim_windows_per_cpu  windows. One array task
    runs one package on sim_cpus CPUs, one window per CPU at a time.
    More tasks than max_array_size are split over several arrays.
Analysis
    The data points that are (or will be) fully simulated, a few per task. The analysis
    starts when all simulation arrays have ended; data points with unfinished windows are skipped.

Options:
    --windows-per-task N    size of a work package (default: sim_cpus x sim_windows_per_cpu)
    --points-per-task N     data points per analysis task (default: from config.toml)
    --only NAME ...         only these data points, e.g. --only iso_00001 iso_00002

Only unfinished work is put in the list, so after a failure or a time-out just run it again.
"""
import argparse
import math
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "..", "apr_pipeline"))

from apr_config import CFG, get_attach_fractions, get_pull_distances

DATA = Path(CFG["directories"]["data"])
LOGS = Path(CFG["directories"]["logs"])
SL = CFG["slurm"]
MARKER = CFG["run"]["done_marker"]
TOPOLOGY = f'{CFG["naming"]["complex_prefix"]}-dum.prmtop'

# window folder names: a000 ... (attach) and p000 ... (pull)
WINDOWS = ([f"a{i:03d}" for i in range(len(get_attach_fractions())-1)] # the last attach window is the first pull window
           + [f"p{i:03d}" for i in range(len(get_pull_distances()))])


# ---------------------------------------------------------------------------
# What state is each data point in, and which of its windows are left?
# ---------------------------------------------------------------------------

def scan(dp):
    """(state, windows still to simulate)"""
    windows = dp / "windows"
    if not windows.is_dir():
        return "no_windows", []
    if not all((windows / w / TOPOLOGY).exists() for w in WINDOWS):
        return "setup_incomplete", []
    left = [w for w in WINDOWS if not (windows / w / MARKER).exists()]
    if left:
        return "needs_simulation", left
    return ("analyzed" if (dp / "results.json").exists() else "needs_analysis"), []


def data_points(only=None):
    dps = sorted(p for p in DATA.iterdir() if p.is_dir() and (p / f"{p.name}.smi").exists())
    if only:
        dps = [p for p in dps if p.name in only]
    return dps


# ---------------------------------------------------------------------------
# Work packages, arrays and time limits
# ---------------------------------------------------------------------------

def hhmmss(minutes):
    minutes = int(math.ceil(minutes))
    return f"{minutes // 60:02d}:{minutes % 60:02d}:00"


def default_package():
    return SL["sim_cpus"] * SL["sim_windows_per_cpu"]


def waves(package):
    """How many windows each CPU runs one after another."""
    return math.ceil(package / SL["sim_cpus"])


def sim_time(package):
    """Time limit: the package needs waves x minutes_per_window; +50 % margin, capped at sim_task_hours."""
    needed = waves(package) * SL["sim_minutes_per_window"]
    return hhmmss(min(SL["sim_task_hours"] * 60, needed * 1.5 + 10))


def analysis_time(per_task):
    return hhmmss(per_task * SL["analysis_minutes_per_data_point"] * 2 + 5)


def split_arrays(n_tasks, max_array):
    """[(offset, size), ...]: an array can hold at most max_array tasks."""
    out, offset = [], 0
    while offset < n_tasks:
        size = min(max_array, n_tasks - offset)
        out.append((offset, size))
        offset += size
    return out


def spread_package(n_windows):
    """Smallest package (a multiple of sim_cpus) that spreads the work over max_running tasks."""
    cpus = SL["sim_cpus"]
    per_task = math.ceil(n_windows / SL["max_running"])
    return max(cpus, math.ceil(per_task / cpus) * cpus)


def estimate_hours(n_windows, package):
    tasks = math.ceil(n_windows / package)
    rounds = math.ceil(tasks / SL["max_running"])
    return rounds * waves(min(package, n_windows)) * SL["sim_minutes_per_window"] / 60


# ---------------------------------------------------------------------------
# sbatch
# ---------------------------------------------------------------------------

def sbatch(name, limit, cpus, mem_flag, stage_args, array=None, dependency=None):
    """array = (size, offset) or None; stage_args = what job_array.slurm gets."""
    cmd = ["sbatch", "--parsable", "--ntasks=1", f"--job-name=apr_{name}", f"--time={limit}",
           f"--cpus-per-task={cpus}", mem_flag, f"--chdir={script_dir}"]
    if SL.get("partition"):
        cmd.append(f"--partition={SL['partition']}")
    if SL.get("account"):
        cmd.append(f"--account={SL['account']}")
    if array:
        size, _ = array
        cmd += [f"--array=0-{size - 1}%{min(size, SL['max_running'])}",
                f"--output={LOGS}/{name}_%A_%a.out"]
    else:
        cmd.append(f"--output={LOGS}/{name}_%j.out")
    if dependency:
        cmd.append("--dependency=afterany:" + ":".join(dependency))   # afterany: also if a task failed
    cmd.append(str(Path(script_dir) / "job_array.slurm"))
    return cmd + [str(a) for a in stage_args]


def run_sbatch(cmd, do_submit, placeholder):
    print("    " + shlex.join(cmd))
    if not do_submit:
        return placeholder
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode:
        if result.stderr:
            print(result.stderr.rstrip(), file=sys.stderr)
        if result.stdout:
            print(result.stdout.rstrip(), file=sys.stderr)
        raise SystemExit(result.returncode)
    job_id = result.stdout.strip().split(";")[0]
    print(f"    -> submitted as job {job_id}")
    return job_id


# ---------------------------------------------------------------------------
# The stages
# ---------------------------------------------------------------------------

def submit_sim(window_list, package, do_submit, stamp, show_hint=True):
    n = len(window_list)
    tasks = math.ceil(n / package)
    arrays = split_arrays(tasks, SL["max_array_size"])
    limit = sim_time(package)

    print(f"\n[sim] {n} windows to simulate")
    print(f"  work package: {package} windows per task on {SL['sim_cpus']} CPUs "
          f"({waves(package)} windows per CPU), time limit {limit}")
    print(f"  -> {tasks} tasks in {len(arrays)} array(s); at most {SL['max_running']} run at once; "
          f"about {estimate_hours(n, package):.1f} h if nodes are free")

    spread = spread_package(n)
    if show_hint and spread < package and tasks < SL["max_running"]:
        print(f"  hint: all {tasks} tasks would start at once, so smaller packages finish sooner:\n"
              f"        --windows-per-task {spread} -> {math.ceil(n / spread)} tasks, "
              f"about {estimate_hours(n, spread) * 60:.0f} min instead of {estimate_hours(n, package):.1f} h")

    worklist = LOGS / f"worklist_sim_{stamp}.txt"        # new file every time: running tasks keep theirs
    if do_submit:
        LOGS.mkdir(parents=True, exist_ok=True)
        worklist.write_text("\n".join(window_list) + "\n")
    ids = []
    for k, (offset, size) in enumerate(arrays):
        print(f"  array {k + 1}/{len(arrays)}: tasks {offset}-{offset + size - 1}")
        cmd = sbatch("sim", limit, SL["sim_cpus"], f"--mem-per-cpu={SL['sim_mem_per_cpu']}",
                     ["sim", worklist, package, offset], array=(size, offset))
        ids.append(run_sbatch(cmd, do_submit, f"<sim_job_{k + 1}>"))
    return ids


def submit_analysis(names, per_task, dependency, do_submit, stamp):
    per_task = max(1, per_task)
    if math.ceil(len(names) / per_task) > SL["max_array_size"]:
        per_task = math.ceil(len(names) / SL["max_array_size"])
    tasks = math.ceil(len(names) / per_task)
    limit = analysis_time(per_task)

    print(f"\n[analysis] {len(names)} data points -> {tasks} tasks of {per_task}, 1 CPU, time limit {limit}")
    worklist = LOGS / f"worklist_analysis_{stamp}.txt"
    if do_submit:
        LOGS.mkdir(parents=True, exist_ok=True)
        worklist.write_text("\n".join(names) + "\n")
    cmd = sbatch("analysis", limit, 1, f"--mem-per-cpu={SL['analysis_mem']}",
                 ["analysis", worklist, per_task, 0], array=(tasks, 0), dependency=dependency)
    return run_sbatch(cmd, do_submit, "<analysis_job>")


def submit_collect(dependency, do_submit):
    print("\n[collect] results.csv is built when the analysis has ended")
    cmd = sbatch("collect", "00:10:00", 1, "--mem-per-cpu=2G", ["collect"], dependency=dependency)
    run_sbatch(cmd, do_submit, "<collect_job>")


def show_status(dps):
    by_state, windows_left = {}, 0
    for dp in dps:
        st, left = scan(dp)
        by_state.setdefault(st, []).append(dp.name)
        windows_left += len(left)
    labels = [("analyzed", "analysis done"),
              ("needs_analysis", "simulated, analysis still to do"),
              ("needs_simulation", "windows ready, simulation not finished"),
              ("setup_incomplete", "setup incomplete (some window folders missing)"),
              ("no_windows", "no windows folder (setup not run or failed)")]
    print(f"{len(dps)} data points in {DATA}")
    for key, label in labels:
        names = by_state.get(key, [])
        if names:
            more = ", ..." if len(names) > 3 else ""
            print(f"  {len(names):6d}  {label}   e.g. {', '.join(names[:3])}{more}")
    print(f"  windows still to simulate: {windows_left}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["status", "sim", "analysis", "all"])
    ap.add_argument("--submit", action="store_true", help="really submit (default: only show)")
    ap.add_argument("--windows-per-task", type=int, help="size of a simulation work package")
    ap.add_argument("--points-per-task", type=int, help="data points per analysis task")
    ap.add_argument("--only", nargs="+", metavar="NAME")
    args = ap.parse_args()

    dps = data_points(args.only)
    if args.command == "status":
        show_status(dps)
        return

    window_list, analysis_names, skipped = [], [], []
    for dp in dps:
        st, left = scan(dp)
        if st == "needs_simulation":
            window_list += [f"{dp.name}/{w}" for w in left]
            if args.command == "all":          # these will need the analysis too
                analysis_names.append(dp.name)
        elif st == "needs_analysis":
            analysis_names.append(dp.name)
        elif st in ("setup_incomplete", "no_windows"):
            skipped.append(dp.name)
    if skipped:
        print(f"note: {len(skipped)} data point(s) without complete windows are left out "
              f"(e.g. {', '.join(skipped[:3])}); see 'status'")

    stamp = time.strftime("%Y%m%d_%H%M%S")
    sim_ids = None
    if args.command in ("sim", "all"):
        if window_list:
            package = args.windows_per_task or default_package()
            sim_ids = submit_sim(window_list, package, args.submit, stamp,
                                 show_hint=args.windows_per_task is None)
        else:
            print("\n[sim] nothing to simulate")

    if args.command in ("analysis", "all"):
        if analysis_names:
            job = submit_analysis(sorted(set(analysis_names)),
                                  args.points_per_task or SL["analysis_data_points_per_task"],
                                  sim_ids, args.submit, stamp)
            submit_collect([job], args.submit)
        else:
            print("\n[analysis] nothing to analyze")

    if not args.submit:
        print("\nNothing was submitted. Add --submit to submit.")


if __name__ == "__main__":
    main()
