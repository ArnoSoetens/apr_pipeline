"""What state is each data point in? Shared by scripts/main.py and scripts/submit_arrays.py.

Setup stages (next_stage() checks the files in the data point folder, latest stage first):
    structure    <name>.smi                          -> raw guest PDB
    placement    any guest PDB (see placement_inputs) -> <guest>.pdb + <complex>.pdb
    parametrize  <guest>.pdb + <complex>.pdb          -> <guest>.mol2 + <guest>.frcmod
    apr_setup    all four guest/complex files         -> windows/*/<complex>-dum.prmtop, complex/
    done         windows are built (scan() then tells simulation/analysis state)

States (checked on the files in the data point folder):
    no_windows        no windows/ folder yet (setup not run)
    setup_incomplete  windows/ exists but some window folder lacks its topology
    setup_failed      like the two above, and setup_failed.txt says why setup failed
    needs_simulation  all windows built; some window is not finished (window_finished: done marker
                      + production trajectory; an interrupted window has only a partial trajectory)
    needs_analysis    every window finished; no results.json yet
    analyzed          results.json exists
"""
import subprocess
import sys
import time
from pathlib import Path

from apr_config import CFG, RESERVED, get_attach_fractions, get_pull_distances

DATA = Path(CFG["directories"]["data"])
TOPOLOGY = f'{CFG["naming"]["complex_prefix"]}-dum.prmtop'
TRAJECTORY = CFG["analysis"]["trajectory"]
FAILED_FILE = "setup_failed.txt"
ANALYSIS_FAILED_FILE = "analysis_failed.txt"
PLOTS = "plots"                                # analysis plots: <data point>/plots/
RESULTS = "results.json"
DONE_MARKER = CFG["run"]["done_marker"]       # run.sh writes it in a window folder when all steps finished
RUN_FAILED_FILE = "run_failed.txt"            # run.sh writes it in a window folder when a step failed
IMPLICIT_DIR = "windows_implicit"             # guest-only runs (scripts/run_implicit.py):
IMPLICIT_ONLY = "implicit_only"               #   <dp>/windows_implicit/implicit_only  free guest in GB
IMPLICIT_P000 = "implicit_p000"               #   <dp>/windows_implicit/implicit_p000  guest alone in the p000 window
# analysis plots in the data point folder itself, from before plots/ existed
LEGACY_PLOT_FILES = ("binding_free_energy.jpg", "convergence_plot.jpg",
                     "attach_windows_grid.jpg", "pull_window_histograms.jpg")
# setup leftovers in the data point folder that cleanup removes after a successful analysis:
# tleap's complex/ build folder and antechamber's scratch files
SETUP_TRANSIENT = ("complex", "sqm.in", "sqm.out", "sqm.pdb", "ATOMTYPE.INF", "ANTECHAMBER*", "leap.log")
SETUP_STATES = ("no_windows", "setup_incomplete", "setup_failed")     # setup has to be (re)done

# window folder names: a000 ... (attach) and p000 ... (pull)
WINDOWS = ([f"a{i:03d}" for i in range(len(get_attach_fractions())-1)] # the last attach window is the first pull window
           + [f"p{i:03d}" for i in range(len(get_pull_distances()))])

# --- setup stages -----------------------------------------------------------
GUEST_PDB = f'{CFG["naming"]["guest_prefix"]}.pdb'
GUEST_MOL2 = f'{CFG["naming"]["guest_prefix"]}.mol2'
GUEST_FRCMOD = f'{CFG["naming"]["guest_prefix"]}.frcmod'
COMPLEX_PDB = f'{CFG["naming"]["complex_prefix"]}.pdb'
CLASH_PDB = f'{CFG["naming"]["guest_prefix"]}_clash.pdb'   # guest pose a placement rejected (clash)
NOT_GUEST_INPUT = {COMPLEX_PDB, CLASH_PDB, "sqm.pdb"}   # PDBs in a data point folder that are never a guest input

STAGES = ("structure", "placement", "parametrize", "apr_setup")
DONE = "done"


def placement_inputs(dp):
    """Guest PDBs placement could read in this folder, sorted.

    Any .pdb counts except the complex PDB, the rejected pose <guest>_clash.pdb and
    antechamber's sqm.pdb. The guest PDB
    (<guest>.pdb) only counts when it is the only candidate: next to another PDB it is
    the output of an earlier (e.g. failed) placement, not the input.
    """
    candidates = sorted(p for p in Path(dp).glob("*.pdb")
                        if p.is_file() and p.name not in NOT_GUEST_INPUT)
    if len(candidates) > 1:
        candidates = [p for p in candidates if p.name != GUEST_PDB]
    return candidates


def start_stage(dp):
    """Where setup starts "from the start": structure with a <name>.smi, otherwise placement."""
    return "structure" if (Path(dp) / f"{Path(dp).name}.smi").is_file() else "placement"


def force_plan(dp, from_stage):
    """What main.py --force deletes before redoing setup from from_stage (existing paths only).

    Always the simulation and analysis output: windows/, complex/, results.json, plots/,
    analysis_failed.txt and old plots in the folder itself. From parametrize or earlier also the
    guest-only runs (windows_implicit/), which were built from the old guest parameters.
    Never inputs (.smi, input PDBs, guest_input.json) or anything else; the setup stages
    overwrite their own outputs (<guest>.pdb, <complex>.pdb, .mol2, .frcmod) themselves.
    """
    dp = Path(dp)
    names = ["windows", "complex", RESULTS, PLOTS, ANALYSIS_FAILED_FILE, *LEGACY_PLOT_FILES]
    if STAGES.index(from_stage) <= STAGES.index("parametrize"):
        names.append(IMPLICIT_DIR)
    return [dp / name for name in names if (dp / name).exists() or (dp / name).is_symlink()]


def has_simulation_output(dp):
    """Whether any window has a production trajectory (windows/*/production.nc).

    Such a data point has simulated data, so setup must never touch it again: not with
    --force, --from or by name. Cleaned-up points keep only a few trajectories (e.g.
    a000-a004), so one trajectory in one window is enough.
    """
    return any((Path(dp) / "windows").glob(f"*/{TRAJECTORY}"))


class HasSimulationOutput(RuntimeError):
    """Raised by a setup stage asked to work on a data point with simulation output."""


def refuse_simulated(dp):
    """First line of every setup stage's run(dp): never touch a data point with simulation output."""
    if has_simulation_output(dp):
        raise HasSimulationOutput(
            f"{Path(dp).name} has simulation output (windows/*/{TRAJECTORY}); setup never touches it")


def next_stage(dp):
    """The first setup stage this data point still needs: one of STAGES, DONE, or None.

    Decided only by which files exist (existing outputs are trusted), latest stage first.
    A data point with simulation output is always DONE. None means the folder has no
    input any stage recognizes.
    """
    dp = Path(dp)
    if has_simulation_output(dp):
        return DONE
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


def is_candidate(path):
    """A folder in the data folder that can be a data point (not inputs/ or logs/)."""
    return path.is_dir() and path.name not in RESERVED


def all_data_points(only=None, data_dir=DATA):
    """Folders in the data directory that any setup stage recognizes (including done ones), sorted."""
    data_dir = Path(data_dir)
    if not data_dir.is_dir():
        return []
    return sorted(p for p in data_dir.iterdir()
                  if is_candidate(p) and (not only or p.name in only) and next_stage(p) is not None)


def ensure_data_directories(data_dir=DATA, smiles_list=None):
    """Return the data points; create them from the SMILES list only when there are none yet.

    The SMILES list is only needed (and only read) when the data directory has no folder
    that any setup stage recognizes.
    """
    dps = all_data_points(data_dir=data_dir)
    if dps:
        return dps
    smiles_list = Path(smiles_list or CFG["directories"]["smiles_list"])
    if not smiles_list.is_file():
        raise FileNotFoundError(
            f"no data point folders in {data_dir}, and no SMILES list to create them from "
            f"({smiles_list}: put it in inputs/ and set input.smiles_list in inputs/config.toml)")
    print(f"[directories] creating data point folders from {smiles_list}", flush=True)
    subprocess.run(
        [sys.executable, str(Path(__file__).with_name("create_data_directories.py")),
         str(smiles_list), str(data_dir)],
        check=True,
    )
    return all_data_points(data_dir=data_dir)


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
    left = [w for w in WINDOWS if not window_finished(windows / w)]
    if left:
        return "needs_simulation", left
    return ("analyzed" if (dp / "results.json").exists() else "needs_analysis"), []


def window_finished(window_dir):
    """A window is finished when it has the done marker and its trajectory (run.sh skips it)."""
    window_dir = Path(window_dir)
    return (window_dir / DONE_MARKER).is_file() and (window_dir / TRAJECTORY).is_file()


def has_complete_run(dp):
    """Whether every configured window is finished (done marker + production trajectory)."""
    windows = Path(dp) / "windows"
    return windows.is_dir() and all(window_finished(windows / window) for window in WINDOWS)


def mark_failed(dp, reason, file=FAILED_FILE):
    (Path(dp) / file).write_text(f"{time.strftime('%F %T')}\n{reason}\n")


def clear_failed(dp, file=FAILED_FILE):
    (Path(dp) / file).unlink(missing_ok=True)


def failure_reason(dp, file=FAILED_FILE):
    lines = (Path(dp) / file).read_text().splitlines()
    return lines[1][:100] if len(lines) > 1 else ""
