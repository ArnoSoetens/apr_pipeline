"""APR binding free energy and sampling diagnostics -> results.json + plots."""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")                 # must come before importing pyplot
import matplotlib.pyplot as plt

from paprika import analysis
from paprika.io import load_restraints

from apr_config import CFG, get_attach_fractions, get_pull_distances, get_method

# --- values from config.toml ---
COMPLEX = CFG["naming"]["complex_prefix"]
INITIAL_DISTANCE = CFG["apr"]["initial_distance"]
METHOD = get_method()                     # "ti-block" or "mbar-autoc"
METHODS = ("ti-block", "mbar-autoc")
ENHANCED = CFG["analysis"]["enhanced"]
CONVERGENCE_PLOT = CFG["analysis"]["convergence_plot"]
BOOT_CYCLES = CFG["analysis"]["bootstrap_cycles"]
TRAJECTORY = CFG["analysis"]["trajectory"]
DISTANCE_THRESHOLD = 7.0
CONVERGENCE_FRACTIONS = np.linspace(0.1, 1.0, 10)


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


def _window_values(simulation_data, phase, restraint_index):
    values = [
        np.asarray(_mag(window[restraint_index]), dtype=float).reshape(-1)
        for window in simulation_data[phase]
    ]
    if not values or any(array.size == 0 for array in values):
        raise ValueError(f"No samples available for the {phase} phase")
    if any(not np.all(np.isfinite(array)) for array in values):
        raise ValueError(f"Non-finite samples found in the {phase} phase")
    return values


def _attachment_samples(simulation_data, attach_window_count):
    values = _window_values(simulation_data, "attach", restraint_index=0)
    if not 0 < attach_window_count <= len(values):
        raise ValueError(
            f"Expected {attach_window_count} attachment windows, found {len(values)}"
        )
    return np.concatenate(values[:attach_window_count])


def _fraction_outside_distance(samples, initial_distance, threshold):
    samples = np.asarray(samples, dtype=float).reshape(-1)
    if samples.size == 0 or not np.all(np.isfinite(samples)):
        raise ValueError("Attachment distance samples must be non-empty and finite")
    return float(np.count_nonzero(np.abs(samples - initial_distance) > threshold) / samples.size)


def _phase_details(res, phase):
    """One total free-energy estimate and SEM per method for the phase."""
    return {
        method: {
            "fe": float(_mag(res[phase][method]["fe"])),
            "sem": float(_mag(res[phase][method]["sem"])),
            "units": "kcal/mol",
        }
        for method in METHODS
    }


def convergence_plot(res, out_path):
    """Plot phase free-energy estimates against the fraction of trajectory frames used."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    for ax, phase in zip(axes, ("attach", "pull")):
        for method in METHODS:
            values = res[phase][method]
            estimates = [float(_mag(values["fraction_fe"][fraction]))
                         for fraction in CONVERGENCE_FRACTIONS]
            errors = [float(_mag(values["fraction_sem"][fraction]))
                      for fraction in CONVERGENCE_FRACTIONS]
            ax.errorbar(
                CONVERGENCE_FRACTIONS * 100,
                estimates,
                yerr=errors,
                marker="o",
                capsize=2,
                label=method,
            )
        ax.set_title(f"{phase.capitalize()} phase")
        ax.set_xlabel("Trajectory frames used (%)")
        ax.set_ylabel("Free energy (kcal/mol)")
        ax.set_xticks(CONVERGENCE_FRACTIONS * 100)
        ax.grid(alpha=0.3)
        ax.legend()
    fig.suptitle("Free-energy convergence with trajectory sampling")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def calculate_free_energy(windows_dir, plot_path=None, grid_path=None, title_name="guest",
                          complex_prefix=COMPLEX, boot_cycles=BOOT_CYCLES, enhanced=ENHANCED,
                          convergence_path=None):
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
    fe.methods = list(METHODS)
    fe.fractions = CONVERGENCE_FRACTIONS.tolist()
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

    if convergence_path:
        convergence_plot(res, convergence_path)

    x_pull = pull_distances - initial_distance + attach_fractions[-1]
    attach_window_count = len(attach_fractions) - 1
    attach_samples = _attachment_samples(sim_data, attach_window_count)
    attach_outside_fraction = _fraction_outside_distance(
        attach_samples, initial_distance, DISTANCE_THRESHOLD
    )

    if plot_path:
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
           "phases": {phase: _phase_details(res, phase) for phase in ("attach", "pull")},
           "attach_distance_outside_threshold": {
                   "metric": "fraction_of_samples",
                   "initial_distance_angstrom": float(initial_distance),
                   "threshold_angstrom": DISTANCE_THRESHOLD,
                   "condition": "abs(sample - initial_distance) > threshold",
                   "fraction": attach_outside_fraction,
                   "percent": 100.0 * attach_outside_fraction,
                   "sample_count": int(attach_samples.size),
                   "windows_included": [f"a{index:03d}" for index in range(attach_window_count)],
           }}
    if enhanced:
        out["method"] = METHOD
        out["attach_fractions"] = attach_fractions
        out["pull_distances"] = [float(x) for x in pull_distances]
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
                                              enhanced=ENHANCED,
                                              convergence_path=(
                                                  dp / "convergence_plot.jpg"
                                                  if CONVERGENCE_PLOT
                                                  else None
                                              ))
    if ENHANCED:
        pull_window_histograms(sim_data, dp / "pull_window_histograms.jpg")
    (dp / "results.json").write_text(json.dumps(results, indent=2))
    return {"dG_bind": results["dG_bind"], "dG_sem": results["dG_sem"]}