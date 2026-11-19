"""What state is each data point in? Shared by scripts/main.py and scripts/submit_arrays.py.

Setup stages (next_stage() checks the files in the data point folder, latest stage first):
    structure    <name>.smi                          -> raw guest PDB
    placement    any guest PDB (see placement_inputs) -> <guest>.pdb + <complex>.pdb
    parametrize  <guest>.pdb + <complex>.pdb          -> <guest>.mol2 + <guest>.frcmod
    apr_setup    all four guest/complex files         -> windows/*/<complex>-dum.prmtop
    done         windows are built (scan() then tells simulation/analysis state)

States (checked on the files in the data point folder):
    no_windows        no windows/ folder yet (setup not run)
    setup_incomplete  windows/ exists but some window folder lacks its topology
    setup_failed      like the two above, and setup_failed.txt says why setup failed
    needs_simulation  all windows built; some window lacks its production trajectory
    needs_analysis    all window production trajectories exist; no results.json yet
    analyzed          results.json exists
"""
import time
from pathlib import Path

from apr_config import CFG, get_attach_fractions, get_pull_distances

DATA = Path(CFG["directories"]["data"])
TOPOLOGY = f'{CFG["naming"]["complex_prefix"]}-dum.prmtop'
TRAJECTORY = CFG["analysis"]["trajectory"]
FAILED_FILE = "setup_failed.txt"
SETUP_STATES = ("no_windows", "setup_incomplete", "setup_failed")     # setup has to be (re)done

# window folder names: a000 ... (attach) and p000 ... (pull)
WINDOWS = ([f"a{i:03d}" for i in range(len(get_attach_fractions())-1)] # the last attach window is the first pull window
           + [f"p{i:03d}" for i in range(len(get_pull_distances()))])

# --- setup stages -----------------------------------------------------------
GUEST_PDB = f'{CFG["naming"]["guest_prefix"]}.pdb'
GUEST_MOL2 = f'{CFG["naming"]["guest_prefix"]}.mol2'
GUEST_FRCMOD = f'{CFG["naming"]["guest_prefix"]}.frcmod'
COMPLEX_PDB = f'{CFG["naming"]["complex_prefix"]}.pdb'
NOT_GUEST_INPUT = {COMPLEX_PDB, "sqm.pdb"}   # PDBs in a data point folder that are never a guest input

STAGES = ("structure", "placement", "parametrize", "apr_setup")
DONE = "done"


def placement_inputs(dp):
    """Guest PDBs placement could read in this folder, sorted.

    Any .pdb counts except the complex PDB and antechamber's sqm.pdb. The guest PDB
    (<guest>.pdb) only counts when it is the only candidate: next to another PDB it is
    the output of an earlier (e.g. failed) placement, not the input.
    """
    candidates = sorted(p for p in Path(dp).glob("*.pdb")
                        if p.is_file() and p.name not in NOT_GUEST_INPUT)
    if len(candidates) > 1:
        candidates = [p for p in candidates if p.name != GUEST_PDB]
    return candidates


def next_stage(dp):
    """The first setup stage this data point still needs: one of STAGES, DONE, or None.

    Decided only by which files exist (existing outputs are trusted), latest stage first.
    None means the folder has no input any stage recognizes.
    """
    dp = Path(dp)
    has = lambda name: (dp / name).is_file()
    if has(GUEST_PDB) and has(COMPLEX_PDB):
        if not (has(GUEST_MOL2) and has(GUEST_FRCMOD)):
            return "parametrize"
        return "apr_setup" if scan(dp)[0] in SETUP_STATES else DONE
    if placement_inputs(dp):
        return "placement"
    if has(f"{dp.name}.smi"):
        return "structure"
    return None


def all_data_points(only=None, data_dir=DATA):
    """Folders in the data directory that any setup stage recognizes (including done ones), sorted."""
    data_dir = Path(data_dir)
    if not data_dir.is_dir():
        return []
    return sorted(p for p in data_dir.iterdir()
                  if p.is_dir() and (not only or p.name in only) and next_stage(p) is not None)


def scan(dp):
    """(state, windows still to simulate)"""
    dp = Path(dp)
    windows = dp / "windows"
    # Cleanup retains only the first window's topology, so recognize that analyzed state
    # before treating the missing per-window topologies as incomplete setup.
    if ((dp / "results.json").exists() and windows.is_dir()
            and (windows / WINDOWS[0] / TOPOLOGY).exists()
            and not all((windows / w / TOPOLOGY).exists() for w in WINDOWS)):
        return "analyzed", []
    if not (windows.is_dir() and all((windows / w / TOPOLOGY).exists() for w in WINDOWS)):
        if (dp / FAILED_FILE).exists():
            return "setup_failed", []
        return ("setup_incomplete" if windows.is_dir() else "no_windows"), []
    left = [w for w in WINDOWS if not (windows / w / TRAJECTORY).is_file()]
    if left:
        return "needs_simulation", left
    return "needs_analysis", []


def data_points(only=None):
    """Folders in the data directory that have a .smi file, sorted."""
    dps = sorted(p for p in DATA.iterdir() if p.is_dir() and (p / f"{p.name}.smi").exists())
    if only:
        dps = [p for p in dps if p.name in only]
    return dps


def prepared_data_points(only=None, data_dir=DATA):
    """Folders with all guest, parameter, and complex files needed for APR setup."""
    data_dir = Path(data_dir)
    required_files = (
        f'{CFG["naming"]["guest_prefix"]}.pdb',
        f'{CFG["naming"]["guest_prefix"]}.mol2',
        f'{CFG["naming"]["guest_prefix"]}.frcmod',
        f'{CFG["naming"]["complex_prefix"]}.pdb',
    )
    dps = sorted(
        path for path in data_dir.iterdir()
        if path.is_dir() and all((path / filename).is_file() for filename in required_files)
    )
    if only:
        dps = [path for path in dps if path.name in only]
    return dps


def setup_data_points(only=None, data_dir=DATA):
    """Data points with SMILES input or complete prepared APR setup inputs."""
    data_dir = Path(data_dir)
    if not data_dir.is_dir():
        return []
    dps = {
        path.name for path in data_dir.iterdir()
        if path.is_dir() and (path / f"{path.name}.smi").is_file()
        and (only is None or path.name in only)
    }
    dps.update(path.name for path in prepared_data_points(only, data_dir))
    return [data_dir / name for name in sorted(dps)]


def has_complete_run(dp):
    """Whether every configured window has its production trajectory."""
    dp = Path(dp)
    windows = dp / "windows"
    return windows.is_dir() and all(
        (windows / window / TRAJECTORY).is_file() for window in WINDOWS
    )


def mark_failed(dp, reason):
    (Path(dp) / FAILED_FILE).write_text(f"{time.strftime('%F %T')}\n{reason}\n")


def clear_failed(dp):
    (Path(dp) / FAILED_FILE).unlink(missing_ok=True)


def failure_reason(dp):
    lines = (Path(dp) / FAILED_FILE).read_text().splitlines()
    return lines[1][:100] if len(lines) > 1 else ""
