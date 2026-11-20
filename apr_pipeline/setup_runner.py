"""The loop over data points shared by scripts/main.py and the four stage scripts.

    python apr_pipeline/<stage>.py              every data point whose next stage is <stage>
    python apr_pipeline/<stage>.py NAME ...     <stage> on these data points, whatever exists

Each stage module (structure, placement, parametrize, apr_setup) has run(dp) for one data point;
the module name is the stage name. A data point with simulation output is always skipped (and
each run(dp) refuses it too). A failure is written to the data point's setup_failed.txt as
"<stage>: <reason>" and the loop goes on; success removes the file.
"""
import argparse
import importlib
import subprocess
from pathlib import Path

from apr_config import CFG
from datapoints import (
    DONE,
    STAGES,
    TRAJECTORY,
    clear_failed,
    has_simulation_output,
    mark_failed,
    next_stage,
)
from guest_input import ensure_guest_input


def describe(error):
    """One-line reason, including the end of a failed command's output."""
    if isinstance(error, subprocess.CalledProcessError):
        output = (error.stderr or error.stdout or "").strip().splitlines()
        return f"{error}" + (f" -- {output[-1]}" if output else "")
    return str(error) or type(error).__name__


def run_data_points(names, plan, label="setup", data_dir=None):
    """Set up each named data point; returns (n_done, n_failed, n_skipped).

    plan(dp) returns the stages to run for dp, in order, or a string: why dp is skipped.
    """
    data_dir = Path(data_dir or CFG["directories"]["data"])
    n_done = n_failed = n_skipped = 0
    for i, name in enumerate(names, 1):
        dp = data_dir / name
        print(f"[{i}/{len(names)}] {name}", flush=True)

        if has_simulation_output(dp):
            todo = f"has simulation output (windows/*/{TRAJECTORY}); setup never touches it"
        elif not dp.is_dir():
            todo = f"no such data point folder: {dp}"
        else:
            todo = plan(dp)
        if isinstance(todo, str):
            print(name, f"SKIPPED: {todo}", flush=True)
            n_skipped += 1
            continue

        stage = "guest_input"
        try:
            ensure_guest_input(dp)
            for stage in todo:
                print(f"  {stage}", flush=True)
                importlib.import_module(stage).run(dp)
        except Exception as e:
            reason = f"{stage}: {describe(e)}"
            print(name, "FAILED in", reason, flush=True)
            mark_failed(dp, reason)
            n_failed += 1
            continue

        clear_failed(dp)
        n_done += 1

    print(f"{label} finished: {n_done} done, {n_failed} failed, {n_skipped} skipped", flush=True)
    return n_done, n_failed, n_skipped


def stage_main(stage, argv=None, data_dir=None):
    """Command line of a stage script. Exit code 1 when a data point failed."""
    later = STAGES[STAGES.index(stage) + 1:]
    parser = argparse.ArgumentParser(
        prog=f"{stage}.py",
        description=f"Run the {stage} stage on data points. Without names: every data point whose "
                    f"next stage is {stage} (its earlier outputs exist, its own do not). With names: "
                    f"{stage} on those, even if its output exists"
                    + (f" (the later stages, {', '.join(later)}, are not redone)." if later else ".")
                    + " Data points with simulation output are always skipped.")
    parser.add_argument("names", nargs="*", metavar="NAME", help="data points (folder names)")
    args = parser.parse_args(argv)
    data_dir = Path(data_dir or CFG["directories"]["data"])

    if args.names:
        names = args.names
    else:
        names = sorted(dp.name for dp in data_dir.iterdir()
                       if dp.is_dir() and next_stage(dp) == stage) if data_dir.is_dir() else []
        if not names:
            print(f"no data point in {data_dir} needs {stage}", flush=True)
            return 0

    def plan(dp):
        before = next_stage(dp)
        redo = before == DONE or (before in STAGES and STAGES.index(before) > STAGES.index(stage))
        if redo and later:
            print(f"  note: {dp.name} already has output of later stages, made from the old {stage} "
                  f"output; to redo those too: python scripts/main.py {dp.name} --from {stage}",
                  flush=True)
        return (stage,)

    n_failed = run_data_points(names, plan, label=stage, data_dir=data_dir)[1]
    return 1 if n_failed else 0
