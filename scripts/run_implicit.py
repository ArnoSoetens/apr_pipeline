"""Guest-only GB simulations: the free guest and/or the guest alone in the p000 APR window.

    python run_implicit.py --implicit-only               every parametrized point without a finished run
    python run_implicit.py --p000 iso_00001 ...          these data points
    python run_implicit.py --implicit-only --p000 ...    both runs
    python run_implicit.py ... iso_00001-iso_00010       a range;  --list todo.txt  points from a file
    python run_implicit.py ... --cpus 4                  4 data points at a time (default 1)
    python run_implicit.py ... --clean                   after a successful run keep only .prmtop,
                                                         production.nc and the done marker
    python run_implicit.py ... --dry-run                 only list what would run

Output: <dp>/windows_implicit/implicit_only/ and <dp>/windows_implicit/implicit_p000/. A finished run
(done marker + production.nc) is never redone; an unfinished one starts again. A point needs its
<guest>.pdb, .mol2 and .frcmod (parametrize). A failed run gets run_failed.txt in its folder; the
others go on. The .in files are written into <data>/inputs/ from config.toml first
(implicit_only/*.in for the free guest; the APR window files for p000). Exit code 1 if a run failed.
"""
import argparse
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir.parent / "apr_pipeline"))   # apr_pipeline/ sits next to scripts/

from apr_config import CFG
from datapoints import IMPLICIT_ONLY, IMPLICIT_P000, all_data_points
from selection import add_selection_args, selected


def run_one(dp, run_name, clean):
    """(point, run, error or None). Runs in a worker process when --cpus > 1."""
    import implicit
    try:
        implicit.run(dp, run_name, clean=clean)
    except Exception as error:
        return Path(dp).name, run_name, str(error) or type(error).__name__
    return Path(dp).name, run_name, None


def main(argv=None):
    import implicit

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_selection_args(parser)
    parser.add_argument("--implicit-only", action="store_true", help="the free guest (windows_implicit/implicit_only)")
    parser.add_argument("--p000", action="store_true", help="the guest alone in the p000 window "
                                                            "(windows_implicit/implicit_p000)")
    parser.add_argument("--cpus", type=int, default=1, metavar="N", help="data points at a time (default 1)")
    parser.add_argument("--clean", action="store_true", help="after a successful run keep only .prmtop, "
                                                             "production.nc and the done marker")
    parser.add_argument("--dry-run", action="store_true", help="only list what would run")
    args = parser.parse_args(argv)
    runs = [name for name, wanted in ((IMPLICIT_ONLY, args.implicit_only), (IMPLICIT_P000, args.p000)) if wanted]
    if not runs:
        parser.error("choose --implicit-only, --p000 or both")
    if args.cpus < 1:
        parser.error("--cpus must be at least 1")

    def needing_a_run():
        return [dp.name for dp in all_data_points()
                if not implicit.missing_inputs(dp) and not all(implicit.finished(dp, r) for r in runs)]

    try:
        names = selected(args, todo=needing_a_run)
    except (FileNotFoundError, ValueError) as e:
        sys.exit(f"run_implicit.py: {e}")

    data_dir = Path(CFG["directories"]["data"])
    jobs, n_skipped = [], 0
    for name in names:
        dp = data_dir / name
        if not dp.is_dir():
            print(f"{name} SKIPPED: no such data point folder", flush=True)
            n_skipped += 1
            continue
        missing = implicit.missing_inputs(dp)
        if missing:
            print(f"{name} SKIPPED: no {', '.join(missing)} (parametrize first)", flush=True)
            n_skipped += 1
            continue
        for run_name in runs:
            if implicit.finished(dp, run_name):
                print(f"{name}/{run_name} SKIPPED: already done", flush=True)
                n_skipped += 1
            else:
                jobs.append((dp, run_name))

    if not jobs:
        print(f"nothing to run ({n_skipped} skipped)", flush=True)
        return 0
    print(f"{len(jobs)} guest-only run(s), {min(args.cpus, len(jobs))} at a time", flush=True)
    if args.dry_run:
        for dp, run_name in jobs:
            print(f"  {dp.name}/{run_name}")
        return 0

    from run_inputs import write_inputs
    write_inputs()

    n_done = n_failed = 0

    def report(name, run_name, error):
        nonlocal n_done, n_failed
        if error:
            n_failed += 1
            print(f"{name}/{run_name} FAILED: {error}", flush=True)
        else:
            n_done += 1
            print(f"{name}/{run_name} done", flush=True)

    if args.cpus == 1:
        for dp, run_name in jobs:
            print(f"{dp.name}/{run_name} started", flush=True)
            report(*run_one(dp, run_name, args.clean))
    else:
        with ProcessPoolExecutor(max_workers=args.cpus) as pool:
            futures = [pool.submit(run_one, dp, run_name, args.clean) for dp, run_name in jobs]
            for future in as_completed(futures):
                report(*future.result())

    print(f"implicit finished: {n_done} done, {n_failed} failed, {n_skipped} skipped", flush=True)
    return 1 if n_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
