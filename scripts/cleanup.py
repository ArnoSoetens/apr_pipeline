"""Remove transient files from analyzed data points, to save disk space.

    python cleanup.py --dry-run                 what would be removed from every analyzed point (nothing is deleted)
    python cleanup.py --dry-run -v iso_00001    the same, listing every file
    python cleanup.py                           clean every analyzed point (one with results.json)
    python cleanup.py iso_00001-iso_00010       these;  --list todo.txt  the data points in a file

What is removed (only from a point with results.json):
    in the data point folder   setup leftovers (datapoints.SETUP_TRANSIENT): complex/, sqm.*, ...
                               everything else is kept: guest_input.json, input PDBs, plots/, your files
    in windows/                everything except restraints.json, the trajectories of the first five
                               windows and the last attach window, and the first window's topology
                               (datapoints.scan() recognizes a cleaned point by that layout)
Window folders themselves are kept, symlinks are never followed, and only folders inside the data
folder are touched.
"""
import argparse
import os
import sys
from pathlib import Path

script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir.parent / "apr_pipeline"))

from apr_config import CFG
from datapoints import RESULTS, is_candidate, SETUP_TRANSIENT, TOPOLOGY, TRAJECTORY, WINDOWS


def kept_in_windows():
    """Paths relative to windows/ that cleanup keeps."""
    last_attach = next(window for window in reversed(WINDOWS) if window.startswith("a"))
    keep = {"restraints.json", f"{WINDOWS[0]}/{TOPOLOGY}"}
    keep |= {f"{window}/{TRAJECTORY}" for window in set(WINDOWS[:5]) | {last_attach}}
    return keep


def files_under(path):
    """Every file below path (path itself if it is a file); symlinks are listed, never followed."""
    if path.is_symlink() or path.is_file():
        return [path]
    found = []
    for root, _, filenames in os.walk(path, followlinks=False):
        found += [Path(root) / name for name in filenames]
    return found


def cleanup_plan(data_point, data_dir=None):
    """The files cleanup would remove from data_point, sorted. Raises for a folder outside data_dir."""
    data_dir = Path(CFG["directories"]["data"] if data_dir is None else data_dir).expanduser().resolve()
    data_point = Path(data_point).resolve()
    if data_point.parent != data_dir:
        raise ValueError(f"refusing to clean a folder outside {data_dir}: {data_point}")
    if not (data_point / RESULTS).is_file():
        return []

    remove = []
    for pattern in SETUP_TRANSIENT:
        for path in data_point.glob(pattern):
            remove += files_under(path)

    windows = data_point / "windows"
    if windows.is_dir() and not windows.is_symlink():
        keep = kept_in_windows()
        remove += [path for path in files_under(windows)
                   if path.relative_to(windows).as_posix() not in keep]
    return sorted(set(remove))


def remove_empty_dirs(data_point):
    """Remove folders left empty, except the windows/ folder and the window folders."""
    windows = data_point / "windows"
    removed = 0
    for root, directories, _ in os.walk(data_point, topdown=False, followlinks=False):
        for directory in directories:
            path = Path(root) / directory
            if path == windows or path.parent == windows or path.is_symlink():
                continue
            if not any(path.iterdir()):
                path.rmdir()
                removed += 1
    return removed


def size(paths):
    return sum(path.lstat().st_size for path in paths)


def human(n_bytes):
    for unit in ("B", "KB", "MB", "GB"):
        if n_bytes < 1024 or unit == "GB":
            return f"{n_bytes:.1f} {unit}" if unit != "B" else f"{n_bytes} B"
        n_bytes /= 1024


def cleanup_data_point(data_point, data_dir=None, dry_run=False, verbose=False):
    """Clean one data point (see the module docstring). Returns (files removed, bytes freed)."""
    data_point = Path(data_point)
    plan = cleanup_plan(data_point, data_dir)
    freed = size(plan)
    if verbose:
        for path in plan:
            print(f"    {path.relative_to(data_point.resolve())}")
    if not dry_run:
        for path in plan:
            path.unlink()
        remove_empty_dirs(data_point.resolve())
    print(f"cleanup {data_point.name}: {'would remove' if dry_run else 'removed'} "
          f"{len(plan)} file(s) ({human(freed)})", flush=True)
    return len(plan), freed


def main(argv=None):
    from selection import add_selection_args, selected

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_selection_args(parser)
    parser.add_argument("--dry-run", action="store_true", help="only show what would be removed")
    parser.add_argument("-v", "--verbose", action="store_true", help="list every file")
    parser.add_argument("--data-dir", type=Path, default=Path(CFG["directories"]["data"]),
                        help="data directory (default: the configured one)")
    args = parser.parse_args(argv)
    data_dir = args.data_dir.expanduser().resolve()
    if not data_dir.is_dir():
        sys.exit(f"cleanup.py: data directory does not exist: {data_dir}")

    def analyzed():
        return sorted(p.name for p in data_dir.iterdir()
                      if is_candidate(p) and not p.is_symlink() and (p / RESULTS).is_file())

    try:
        names = selected(args, todo=analyzed, data_dir=data_dir)
    except (FileNotFoundError, ValueError) as e:
        sys.exit(f"cleanup.py: {e}")

    total_files = total_bytes = 0
    for name in names:
        data_point = data_dir / name
        if not (data_point / RESULTS).is_file():
            print(f"cleanup {name}: SKIPPED: not analyzed (no {RESULTS}); nothing removed", flush=True)
            continue
        n_files, n_bytes = cleanup_data_point(data_point, data_dir, args.dry_run, args.verbose)
        total_files += n_files
        total_bytes += n_bytes
    print(f"cleanup finished: {'would remove' if args.dry_run else 'removed'} {total_files} file(s) "
          f"({human(total_bytes)}) from {len(names)} data point(s)", flush=True)


if __name__ == "__main__":
    main()
