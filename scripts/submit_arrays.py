#!/usr/bin/env python3
"""Submit setup, simulations and analysis as slurm arrays.

    python submit_arrays.py status                 # what is finished and what is left
    python submit_arrays.py setup --submit         # main.py for the data points that need it
    python submit_arrays.py all                    # show what would be submitted (nothing is submitted)
    python submit_arrays.py all --submit           # setup -> simulations -> analysis -> results.csv
    python submit_arrays.py sim --submit           # only the simulations
    python submit_arrays.py analysis --submit      # only the analysis

Setup
    If no data point folders exist yet, they are created from the XML automatically. Data points
    without complete windows are cut into tasks of a few data points; each task runs main.py for its
    data points. A failed data point gets a setup_failed.txt and is left out next time (--retry-failed).
    The 'all' command chains setup, simulations, analysis and collection. Setting up serially is
    also possible: python scripts/main.py.
Simulation
    Every window that still has to be simulated goes into ONE list (lines like iso_00001/a000).
    The list is cut into work packages of  sim_cpus x sim_windows_per_cpu  windows. One array task
    runs one package on sim_cpus CPUs, one window per CPU at a time.
    More tasks than max_array_size are split over several arrays.
Analysis
    The data points that are (or will be) fully simulated, a few per task. The analysis
    starts when all simulation arrays have ended; data points with unfinished windows are skipped.

Options:
    --windows-per-task N    size of a simulation work package (default: sim_cpus x sim_windows_per_cpu)
    --points-per-task N     data points per setup or analysis task (default: from config.toml)
    --retry-failed          setup: include data points that failed before
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

from apr_config import CFG
from datapoints import DATA, FAILED_FILE, WINDOWS, data_points, failure_reason, scan

LOGS = Path(CFG["directories"]["logs"])
SL = CFG["slurm"]


def ensure_data_directories():
    """Create data point folders from the XML when none exist yet."""
    if DATA.is_dir() and data_points():
        return

    pipeline_dir = Path(CFG["directories"]["pipeline"])
    creator = pipeline_dir / "create_data_directories.py"
    print(f"\n[directories] creating data point folders from {CFG['directories']['input_xml']}")
    subprocess.run(
        [sys.executable, str(creator), CFG["directories"]["input_xml"], str(DATA)],
        cwd=pipeline_dir,
        check=True,
    )


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
    cmd = ["sbatch", "--parsable","--ntasks=1", f"--job-name=apr_{name}", f"--time={limit}",
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

def submit_sim(window_list, package, do_submit, stamp, show_hint=True, dependency=None):
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
                     ["sim", worklist, package, offset], array=(size, offset), dependency=dependency)
        ids.append(run_sbatch(cmd, do_submit, f"<sim_job_{k + 1}>"))
    return ids


def points_resources(stage):
    """(cpus, memory flag) of the stages that work on whole data points."""
    if stage == "setup":
        return SL["setup_cpus"], f"--mem-per-cpu={SL['setup_mem']}"
    return 1, f"--mem-per-cpu={SL['analysis_mem']}"


def points_time(stage, per_task):
    if stage == "setup":
        return hhmmss(per_task * SL["setup_minutes_per_data_point"] * 1.5 + 10)
    return hhmmss(per_task * SL["analysis_minutes_per_data_point"] * 2 + 5)


def submit_points(stage, names, per_task, dependency, do_submit, stamp):
    """Array with one worklist line per data point: setup and analysis."""
    per_task = max(1, per_task)
    if math.ceil(len(names) / per_task) > SL["max_array_size"]:
        per_task = math.ceil(len(names) / SL["max_array_size"])
    tasks = math.ceil(len(names) / per_task)
    cpus, mem = points_resources(stage)
    limit = points_time(stage, per_task)

    print(f"\n[{stage}] {len(names)} data points -> {tasks} tasks of {per_task}, {cpus} CPU(s), time limit {limit}")
    rounds = math.ceil(tasks / SL["max_running"])
    print(f"  at most {SL['max_running']} tasks run at once -> about {rounds} round(s) of "
          f"~{per_task * SL[stage + '_minutes_per_data_point']} min if nodes are free")
    worklist = LOGS / f"worklist_{stage}_{stamp}.txt"
    if do_submit:
        LOGS.mkdir(parents=True, exist_ok=True)
        worklist.write_text("\n".join(names) + "\n")
    cmd = sbatch(stage, limit, cpus, mem, [stage, worklist, per_task, 0],
                 array=(tasks, 0), dependency=dependency)
    return run_sbatch(cmd, do_submit, f"<{stage}_job>")


def submit_collect(dependency, do_submit):
    print("\n[collect] results.csv is built when the analysis has ended")
    cmd = sbatch("collect", "00:10:00", 1, "--mem-per-cpu=1G", ["collect"], dependency=dependency)
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
              ("setup_failed", "setup failed (setup_failed.txt says why; 'setup --retry-failed')"),
              ("no_windows", "no windows folder yet (setup not run)")]
    print(f"{len(dps)} data points in {DATA}")
    for key, label in labels:
        names = by_state.get(key, [])
        if names:
            more = ", ..." if len(names) > 3 else ""
            print(f"  {len(names):6d}  {label}   e.g. {', '.join(names[:3])}{more}")
    for name in by_state.get("setup_failed", [])[:3]:
        print(f"           {name}: {failure_reason(DATA / name)}")
    print(f"  windows still to simulate: {windows_left}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["status", "setup", "sim", "analysis", "all"])
    ap.add_argument("--submit", action="store_true", help="really submit (default: only show)")
    ap.add_argument("--windows-per-task", type=int, help="size of a simulation work package")
    ap.add_argument("--points-per-task", type=int, help="data points per setup or analysis task")
    ap.add_argument("--retry-failed", action="store_true", help="setup: include data points that failed before")
    ap.add_argument("--only", nargs="+", metavar="NAME")
    args = ap.parse_args()

    if args.command != "status":
        ensure_data_directories()

    dps = data_points(args.only)
    if args.command == "status":
        show_status(dps)
        return

    window_list, analysis_names, setup_names, failed_names = [], [], [], []
    for dp in dps:
        st, left = scan(dp)
        has_failed_marker = (dp / FAILED_FILE).exists()
        if has_failed_marker:
            failed_names.append(dp.name)
            if args.command == "all" and args.retry_failed:
                setup_names.append(dp.name)
                window_list.extend(f"{dp.name}/{window}" for window in WINDOWS)
                analysis_names.append(dp.name)
        elif st == "needs_simulation":
            window_list += [f"{dp.name}/{w}" for w in left]
            if args.command == "all":
                analysis_names.append(dp.name)
        elif st == "needs_analysis":
            analysis_names.append(dp.name)
        elif st in ("no_windows", "setup_incomplete"):
            setup_names.append(dp.name)
            if args.command == "all":
                window_list.extend(f"{dp.name}/{window}" for window in WINDOWS)
                analysis_names.append(dp.name)
        elif st == "setup_failed":
            failed_names.append(dp.name)
    if args.command in ("sim", "analysis"):
        left_out = setup_names + failed_names
        if left_out:
            print(f"note: {len(left_out)} data point(s) without complete windows are left out "
                  f"(e.g. {', '.join(left_out[:3])}); run 'setup' first or see 'status'")
    elif args.command == "all" and failed_names and not args.retry_failed:
        print(f"note: {len(failed_names)} setup-failed data point(s) are left out "
              f"(e.g. {', '.join(failed_names[:3])}); --retry-failed includes them")

    stamp = time.strftime("%Y%m%d_%H%M%S")
    setup_ids = None
    if args.command == "setup":
        names = setup_names + (failed_names if args.retry_failed else [])
        if names:
            submit_points("setup", sorted(names), args.points_per_task or SL["setup_data_points_per_task"],
                          None, args.submit, stamp)
        elif not dps:
            print(f"\nno data point folders were created in {DATA}; check the configured input XML")
        else:
            print("\n[setup] nothing to set up")
        if failed_names and not args.retry_failed:
            print(f"note: {len(failed_names)} data point(s) failed earlier and are left out "
                  f"(e.g. {', '.join(failed_names[:3])}); --retry-failed includes them")
    elif args.command == "all":
        if setup_names:
            setup_ids = [submit_points(
                "setup", sorted(set(setup_names)),
                args.points_per_task or SL["setup_data_points_per_task"],
                None, args.submit, stamp,
            )]
        else:
            print("\n[setup] nothing to set up")

    sim_ids = None
    if args.command in ("sim", "all"):
        if window_list:
            package = args.windows_per_task or default_package()
            sim_ids = submit_sim(window_list, package, args.submit, stamp,
                                 show_hint=args.windows_per_task is None,
                                 dependency=setup_ids)
        else:
            print("\n[sim] nothing to simulate")

    if args.command in ("analysis", "all"):
        if analysis_names:
            dependencies = sim_ids if sim_ids is not None else setup_ids
            job = submit_points("analysis", sorted(set(analysis_names)),
                                  args.points_per_task or SL["analysis_data_points_per_task"],
                                  dependencies, args.submit, stamp)
            submit_collect([job], args.submit)
        else:
            print("\n[analysis] nothing to analyze")

    if not args.submit:
        print("\nNothing was submitted. Add --submit to submit.")


if __name__ == "__main__":
    main()
