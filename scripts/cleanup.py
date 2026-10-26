"""Remove transient files from a successfully analyzed data point."""
import os
import sys
from pathlib import Path

script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir.parent / "apr_pipeline"))

from apr_config import CFG
from datapoints import WINDOWS


PLOT_FILES = {
    "pull_window_histograms.jpg",
    "attach_windows_grid.jpg",
    "binding_free_energy.jpg",
    "convergence_plot.jpg",
}


def cleanup_data_point(data_point):
    data_dir = Path(CFG["directories"]["data"]).resolve()
    data_point = Path(data_point).resolve()
    if data_point.parent != data_dir or data_point == data_dir:
        raise ValueError(f"Refusing to clean data point outside {data_dir}: {data_point}")

    naming = CFG["naming"]
    retained_root_files = {
        f"{naming['guest_prefix']}.{suffix}"
        for suffix in ("pdb", "mol2", "frcmod")
    }
    retained_root_files.update({
        f"{naming['complex_prefix']}.pdb",
        f"{data_point.name}.smi",
        "results.json",
        *PLOT_FILES,
    })
    retained_windows = set(WINDOWS[:5])
    retained_topology = f"{naming['complex_prefix']}-dum.prmtop"

    removed = 0
    for root, _, filenames in os.walk(data_point, followlinks=False):
        for filename in filenames:
            path = Path(root) / filename
            relative = path.relative_to(data_point)
            keep = relative.as_posix() in retained_root_files

            if relative == Path("windows/restraints.json"):
                keep = True
            elif relative.parts and relative.parts[0] == "windows":
                if len(relative.parts) >= 3:
                    window_name = relative.parts[1]
                    keep = (
                        window_name in retained_windows and path.suffix.lower() == ".nc"
                    ) or (
                        window_name == "a000"
                        and relative.as_posix()
                        == f"windows/a000/{retained_topology}"
                    )

            if not keep:
                path.unlink()
                removed += 1

    print(f"cleanup {data_point.name}: removed {removed} file(s)")
