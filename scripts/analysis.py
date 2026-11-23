"""Analyze simulated data points: binding free energy -> <dp>/results.json, plots in <dp>/plots/.

    python analysis.py                         every data point that is fully simulated and not analyzed yet
    python analysis.py iso_00001 ...           these data points, also when they were analyzed before
    python analysis.py iso_00001-iso_00010     a range;  --list todo.txt  the data points in a file
    python analysis.py --enhanced ...          also the convergence plot, the attach-window grid, the
                                               pull-window histograms and extra metadata in results.json
    python analysis.py --no-enhanced ...       without those (the default is analysis.enhanced in config.toml)

A data point that fails gets analysis_failed.txt (date and reason) and the loop goes on; a
successful analysis removes it. At the end results.csv is rebuilt from every results.json
(analysis.results_csv in config.toml); --no-collect leaves it alone (the slurm array tasks use
this, a separate job collects afterwards). The exit code is 0 even when a data point failed:
see the summary line and the analysis_failed.txt files.
"""
import argparse
import sys
from pathlib import Path

script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir.parent / "apr_pipeline"))   # apr_pipeline/ sits next to scripts/
sys.path.insert(0, str(script_dir))

from apr_config import CFG
from datapoints import (
    ANALYSIS_FAILED_FILE,
    DONE_MARKER,
    TOPOLOGY,
    TRAJECTORY,
    WINDOWS,
    all_data_points,
    clear_failed,
    mark_failed,
    scan,
    window_finished,
)
from selection import add_selection_args, selected

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
add_selection_args(parser)
parser.add_argument("--enhanced", action=argparse.BooleanOptionalAction, default=None,
                    help="also write the convergence plot, attach-window grid, pull-window histograms "
                         "and extra metadata (default: analysis.enhanced in config.toml)")
parser.add_argument("--clean-after-analysis", action="store_true",
                    help="remove transient files after each successful data-point analysis")
parser.add_argument("--no-collect", action="store_true",
                    help="do not rebuild results.csv at the end")
args = parser.parse_args()
enhanced = CFG["analysis"]["enhanced"] if args.enhanced is None else args.enhanced


def needing_analysis():
    names = [dp.name for dp in all_data_points() if scan(dp)[0] == "needs_analysis"]
    if not names:
        print("no data point needs analysis", flush=True)
    return names


def why_not(dp):
    """Why dp cannot be analyzed, or None when it can."""
    windows = dp / "windows"
    if not dp.is_dir():
        return f"no such data point folder: {dp}"
    if (dp / "results.json").is_file() and not all((windows / w / TOPOLOGY).is_file() for w in WINDOWS):
        return "cleaned up after analysis: its trajectories were removed"
    missing = [w for w in WINDOWS if not window_finished(windows / w)]
    if missing:
        return (f"not all windows are finished ({TRAJECTORY} + {DONE_MARKER}; "
                f"{len(missing)} missing, e.g. {', '.join(missing[:3])})")
    return None


try:
    names = selected(args, todo=needing_analysis)
except (FileNotFoundError, ValueError) as e:
    sys.exit(f"analysis.py: {e}")

if names:
    import apr_analysis                       # pAPRika/JAX: only imported when there is work
    if args.clean_after_analysis:
        from cleanup import cleanup_data_point

data_dir = Path(CFG["directories"]["data"])
n_done = n_failed = n_skipped = 0
for i, name in enumerate(names, 1):
    dp = data_dir / name
    print(f"[{i}/{len(names)}] {name}", flush=True)
    reason = why_not(dp)
    if reason:
        print(name, f"SKIPPED: {reason}", flush=True)
        n_skipped += 1
        continue

    try:
        out = apr_analysis.run(dp, enhanced=enhanced)
    except Exception as e:
        reason = str(e) or type(e).__name__
        print(name, f"FAILED: {reason}", flush=True)
        mark_failed(dp, reason, file=ANALYSIS_FAILED_FILE)
        n_failed += 1
        continue
    clear_failed(dp, file=ANALYSIS_FAILED_FILE)
    n_done += 1
    print(f"{name} done: {out['dG_bind']:.2f} +/- {out['dG_sem']:.2f} kcal/mol", flush=True)

    if args.clean_after_analysis:
        try:
            cleanup_data_point(dp)
        except Exception as e:
            print(f"cleanup {name} FAILED: {e}", flush=True)

print(f"analysis finished: {n_done} done, {n_failed} failed, {n_skipped} skipped", flush=True)

if not args.no_collect:
    import collect_results
    collect_results.main()
