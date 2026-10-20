"""Analysis: APR binding free energy (TI, block bootstrap) -> results.json + plots."""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")                 # must come before importing pyplot
import matplotlib.pyplot as plt

from paprika import analysis
from paprika.io import load_restraints

from APR.apr_refactor.apr_pipeline.apr_config import CFG, get_attach_fractions, get_pull_distances, get_method

# --- values from config.toml ---
COMPLEX = CFG["naming"]["complex_prefix"]
INITIAL_DISTANCE = CFG["apr"]["initial_distance"]
METHOD = get_method()                     # "ti-block" or "mbar-autoc"
ENHANCED = CFG["analysis"]["enhanced"]
BOOT_CYCLES = CFG["analysis"]["bootstrap_cycles"]
TRAJECTORY = CFG["analysis"]["trajectory"]


def pull_window_histograms(simulation_data, out_path, restraint_index=0):
    for values in (w[restraint_index] for w in simulation_data["pull"]):
        plt.hist(values.magnitude, bins=50, alpha=0.5)
    plt.xlabel("Distance (Å)")
    plt.savefig(out_path, dpi=300)
    plt.close()


def attach_window_grid(fe, out_path, restraint_index=0, phase="attach", ncols=5):
    """One histogram per attach window: sampled restraint values, with the target as a red line.

    restraint_index 0 is the distance restraint (Å). For the angle/torsion restraints
    (1, 2) the values may need converting from radians to degrees.
    """
    data = fe.simulation_data[phase]                 # list per window, per restraint
    num_win = len(data)

    # ordered targets for this restraint, in the same window order as simulation_data
    _, _, _, _, _, ordered_targets, _ = fe.prepare_data(phase)
    targets = ordered_targets[restraint_index].magnitude

    nrows = int(np.ceil(num_win / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 2.5 * nrows),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes).flatten()

    for k in range(num_win):
        ax = axes[k]
        values = data[k][restraint_index].magnitude
        ax.hist(values, bins=40, color="steelblue", alpha=0.8)
        ax.axvline(targets[k], color="red", linestyle="--", linewidth=1)
        ax.set_title(f"a{k:03d}", fontsize=9)
        ax.tick_params(labelsize=7)

    for k in range(num_win, len(axes)):              # switch off unused slots
        axes[k].axis("off")

    fig.suptitle("Attach phase: restraint value distribution per window", y=1.0)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def _mag(x):
    """Strip pint units and give plain numbers / lists (so json can store them)."""
    return np.asarray(getattr(x, "magnitude", x)).tolist()


def _phase_details(res, phase):
    """Per-phase numbers for the enhanced results.json."""
    r = res[phase][METHOD]
    details = {"fe": _mag(r["fe"]), "sem": _mag(r["sem"]),
               "fe_per_window": _mag(r["fe_matrix"][0, :])}     # cumulative, window 0 -> window k
    if "sem_matrix" in r:
        details["sem_per_window"] = _mag(r["sem_matrix"][0, :])
    if "window_order" in res[phase]:
        details["window_order"] = [str(w) for w in res[phase]["window_order"]]
    return details


def calculate_free_energy(windows_dir, plot_path=None, grid_path=None, title_name="guest",
                          complex_prefix=COMPLEX, boot_cycles=BOOT_CYCLES, enhanced=ENHANCED):
    """Returns (results dict, simulation_data). dG values in kcal/mol."""

    attach_fractions = get_attach_fractions()
    initial_distance = INITIAL_DISTANCE
    pull_distances = get_pull_distances()
    guest_restraints = load_restraints(f"{windows_dir}/restraints.json")

    fe = analysis.fe_calc()
    fe.topology = f"{complex_prefix}-dum.prmtop"
    fe.trajectory = TRAJECTORY
    fe.path = str(windows_dir)
    fe.restraint_list = guest_restraints
    fe.collect_data()
    fe.methods = [METHOD]
    fe.ti_matrix = "full"
    fe.boot_cycles = boot_cycles
    fe.compute_free_energy()
    sim_data = fe.simulation_data

    fe.compute_ref_state_work([guest_restraints[0], guest_restraints[1], None, None, guest_restraints[2], None,])

    res = fe.results
    binding_affinity = -1 * (res["attach"][METHOD]["fe"] + res["pull"][METHOD]["fe"]
                             + res["ref_state_work"])
    sem = np.sqrt(res["attach"][METHOD]["sem"] ** 2 + res["pull"][METHOD]["sem"] ** 2)

    attach_fe = res["attach"][METHOD]["fe_matrix"][0, :].magnitude
    pull_fe = res["pull"][METHOD]["fe_matrix"][0, :].magnitude
    ref_work = res["ref_state_work"].magnitude
    dg, dg_sem = float(binding_affinity.magnitude), float(sem.magnitude)

    if plot_path:
        x_pull = pull_distances - initial_distance + attach_fractions[-1]
        plt.plot(attach_fractions, attach_fe, label="attach phase")
        plt.plot(x_pull, attach_fe[-1] + pull_fe, label="pull phase")
        plt.plot([x_pull[-1]] * 2, [attach_fe[-1] + pull_fe[-1], attach_fe[-1] + pull_fe[-1] + ref_work],
                 label="release phase")
        plt.ylabel("Work (kcal/mol)")
        plt.title(f"Binding affinity of {title_name} = {dg:0.2f} +/- {dg_sem:0.2f} kcal/mol")
        plt.legend()
        plt.savefig(plot_path, dpi=400)
        plt.close()

    if grid_path:
        attach_window_grid(fe, grid_path)

    out = {"dG_bind": dg, "dG_sem": dg_sem, "units": "kcal/mol",
           "ref_state_work": float(ref_work),
           "attach_fe": [float(x) for x in attach_fe], "pull_fe": [float(x) for x in pull_fe]}
    if enhanced:
        out["method"] = METHOD
        out["attach_fractions"] = attach_fractions
        out["pull_distances"] = [float(x) for x in pull_distances]
        out["phases"] = {phase: _phase_details(res, phase) for phase in ("attach", "pull")}
    return out, sim_data


def run(dp, complex_prefix=COMPLEX):
    dp = Path(dp)
    windows = dp / "windows"
    print(windows)

    # every window folder must have finished its production run
    missing = [w.name for w in sorted(windows.iterdir())
               if w.is_dir() and not (w / TRAJECTORY).exists()]
    if missing:
        raise FileNotFoundError(f"{len(missing)} window(s) without {TRAJECTORY}, e.g. {missing[:3]}")

    results, sim_data = calculate_free_energy(windows, plot_path=dp / "binding_free_energy.jpg",
                                              grid_path=(dp / "attach_windows_grid.jpg") if ENHANCED else None,
                                              title_name=dp.name, complex_prefix=complex_prefix,
                                              enhanced=ENHANCED)
    if ENHANCED:
        pull_window_histograms(sim_data, dp / "pull_window_histograms.jpg")
    (dp / "results.json").write_text(json.dumps(results, indent=2))
    return {"dG_bind": results["dG_bind"], "dG_sem": results["dG_sem"]}