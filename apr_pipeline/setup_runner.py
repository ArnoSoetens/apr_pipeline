"""The loop over data points shared by scripts/main.py and the four stage scripts.

    python apr_pipeline/<stage>.py              every data point whose next stage is <stage>
    python apr_pipeline/<stage>.py NAME ...     <stage> on these data points, whatever exists
                                                (also ranges and --list FILE, see selection.py)

Each stage module (structure, placement, parametrize, apr_setup) has run(dp) for one data point;
the module name is the stage name. A data point with simulation output is always skipped (and
each run(dp) refuses it too). A failure is written to the data point's setup_failed.txt as
"<stage>: <reason>" and the loop goes on; success removes the file.
"""
import argparse
import importlib
import shutil
import subprocess
import sys
from pathlib import Path

from apr_config import CFG
from datapoints import (
    DONE,
    STAGES,
    TRAJECTORY,
    clear_failed,
    has_simulation_output,
    is_candidate,
    mark_failed,
    next_stage,
)
from guest_input import ensure_guest_input
from selection import add_selection_args, selected


def describe(error):
    """One-line reason, including the end of a failed command's output."""
    if isinstance(error, subprocess.CalledProcessError):
        output = (error.stderr or error.stdout or "").strip().splitlines()
        return f"{error}" + (f" -- {output[-1]}" if output else "")
    return str(error) or type(error).__name__


def run_data_points(names, plan, label="setup", data_dir=None, prepare=None):
    """Set up each named data point; returns (n_done, n_failed, n_skipped).

    plan(dp) returns the stages to run for dp, in order, or a string: why dp is skipped.
    prepare(dp), if given, runs first for each existing data point (main.py --force deletes there).
    """
    data_dir = Path(data_dir or CFG["directories"]["data"])
    n_done = n_failed = n_skipped = 0
    for i, name in enumerate(names, 1):
        dp = data_dir / name
        print(f"[{i}/{len(names)}] {name}", flush=True)
        if prepare is not None and dp.is_dir():
            prepare(dp)

        if has_simulation_output(dp):
            todo = (f"has simulation output (windows/*/{TRAJECTORY}); setup never touches it "
                    "(main.py --force redoes it and deletes that output)")
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


def remove(paths):
    """Delete files and folders (a symlink itself, never what it points to)."""
    for path in paths:
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        else:
            path.unlink()


def disk_usage(paths):
    total = 0
    for path in paths:
        if path.is_dir() and not path.is_symlink():
            total += sum(f.lstat().st_size for f in path.rglob("*") if f.is_file() or f.is_symlink())
        else:
            total += path.lstat().st_size
    return total


def confirm(question, assume_yes=False, stdin=None, ask=input):
    """True when the user agrees. --yes agrees; without a terminal and without --yes: SystemExit."""
    if assume_yes:
        return True
    stdin = stdin or sys.stdin
    if not stdin.isatty():
        raise SystemExit("main.py: --force needs a confirmation: run it in a terminal, or add --yes")
    return ask(f"{question} [y/N] ").strip().lower() in ("y", "yes")


def stage_main(stage, argv=None, data_dir=None):
    """Command line of a stage script. Exit code 1 when a data point failed."""
    later = STAGES[STAGES.index(stage) + 1:]
    parser = argparse.ArgumentParser(
        prog=f"{stage}.py",
        description=f"Run the {stage} stage on data points. Without names: every data point whose "
                    f"next stage is {stage} (its earlier outputs exist, its own do not). With names: "
                    f"{stage} on those, even if its output exists"
                    + (f" (the later stages, {', '.join(later)}, are not redone)." if later else ".")
                    + " Data points with simulation output are always skipped. Names can be "
                    "ranges (iso_00001-iso_00010); --list FILE reads more from a file.")
    add_selection_args(parser)
    args = parser.parse_args(argv)
    data_dir = Path(data_dir or CFG["directories"]["data"])

    def needing_stage():
        names = sorted(dp.name for dp in data_dir.iterdir()
                       if is_candidate(dp) and next_stage(dp) == stage) if data_dir.is_dir() else []
        if not names:
            print(f"no data point in {data_dir} needs {stage}", flush=True)
        return names

    try:
        names = selected(args, todo=needing_stage, data_dir=data_dir)
    except (FileNotFoundError, ValueError) as e:
        parser.error(str(e))
    if not names:
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
