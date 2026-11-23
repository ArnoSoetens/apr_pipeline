"""Which data points (and windows) a script works on. Shared by every script's command line.

    script.py                          every data point that needs this script (the script decides)
    script.py iso_00001 iso_00002      these data points
    script.py iso_00001-iso_00010      a range: same prefix, numbers counted with the same zero-padding;
                                       only numbers that have a folder (a note lists the others)
    script.py iso_00001/a000           one window          (only scripts that take windows)
    script.py iso_00001/a000-a005      a range of windows  (only scripts that take windows)
    script.py --list todo.txt          the same entries, one or more per line; blank lines and
                                       '# comments' are ignored

Entries may be mixed. Duplicates are dropped, the order is kept. A plain name without a folder
is kept, so the script reports it ("no such data point folder").
"""
import re
from pathlib import Path

from apr_config import CFG

RANGE = re.compile(r"^(?P<prefix>.*?)(?P<first>\d+)-(?P=prefix)(?P<last>\d+)$")


def add_selection_args(parser, windows=False):
    """Add NAME ... and --list FILE to an argparse parser."""
    what = "data points (folder names), ranges like iso_00001-iso_00010"
    if windows:
        what += ", or windows like iso_00001/a000 and iso_00001/a000-a005"
    parser.add_argument("names", nargs="*", metavar="NAME", help=f"{what} (default: all that need it)")
    parser.add_argument("--list", dest="list_file", metavar="FILE",
                        help="a text file with more of these entries; blank lines and # comments are ignored")


def expand_range(entry):
    """iso_00001-iso_00003 -> [iso_00001, iso_00002, iso_00003]; None when entry is not a range."""
    match = RANGE.match(entry)
    if not match:
        return None
    prefix, first, last = match["prefix"], match["first"], match["last"]
    if int(last) < int(first):
        raise ValueError(f"range {entry!r} counts down; write the lowest number first")
    if len(first) == len(last):
        width = len(first)
    elif not first.startswith("0") and not last.startswith("0"):
        width = 0                                   # p1-p10: numbers without zero-padding
    else:
        raise ValueError(f"range {entry!r}: both ends need the same number of digits")
    return [f"{prefix}{number:0{width}d}" for number in range(int(first), int(last) + 1)]


def check_name(name):
    if not name or Path(name).name != name or name in (".", ".."):
        raise ValueError(f"not a data point folder name: {name!r}")
    return name


def entries_from(args):
    """The raw entries: the NAME arguments, then the --list file's entries."""
    entries = list(args.names)
    if getattr(args, "list_file", None):
        path = Path(args.list_file)
        if not path.is_file():
            raise FileNotFoundError(f"--list file not found: {path}")
        for line in path.read_text().splitlines():
            entries.extend(line.split("#", 1)[0].split())
    return entries


def resolve(entries, data_dir=None, windows=False):
    """Expand entries into data point names, or (name, window or None) pairs when windows=True."""
    data_dir = Path(data_dir or CFG["directories"]["data"])
    out = []
    for entry in entries:
        name, _, window = entry.partition("/")
        if window and not windows:
            raise ValueError(f"{entry!r}: this script works on whole data points, not windows")

        if (data_dir / name).is_dir():              # a folder name that only looks like a range
            names = [name]
        else:
            names = expand_range(name)
            if names is None:
                names = [check_name(name)]
            else:
                present = [n for n in names if (data_dir / n).is_dir()]
                missing = [n for n in names if n not in present]
                if missing:
                    shown = ", ".join(missing[:5]) + (", ..." if len(missing) > 5 else "")
                    print(f"note: {name}: no folder for {len(missing)} of {len(names)} "
                          f"({shown})", flush=True)
                names = present
        window_names = [None]
        if window:
            window_names = expand_range(window) or [check_name(window)]
        for n in names:
            check_name(n)
            for w in window_names:
                out.append((n, w) if windows else n)
    return list(dict.fromkeys(out))                 # drop duplicates, keep the order


def selected(args, todo, data_dir=None, windows=False):
    """The script's selection: the entries given (NAME ... and --list), or todo() when none were.

    Returns data point names, or (name, window or None) pairs when windows=True. todo() returns
    the script's default selection in the same form.
    """
    if not args.names and not getattr(args, "list_file", None):
        return todo()
    return resolve(entries_from(args), data_dir=data_dir, windows=windows)
