"""Analysis: APR binding free energy (TI, block bootstrap) -> results.json + plots."""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")                 # must come before importing pyplot
import matplotlib.pyplot as plt

from paprika import analysis
from paprika.io import load_restraints


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


def calculate_free_energy(windows_dir, plot_path=None, grid_path=None, title_name="guest",
                          host_prefix="BCD", guest_prefix="UNL", boot_cycles=1000):
    """Returns (results dict, simulation_data). dG values in kcal/mol."""

    attach_string = "0.00 0.01 0.0375 0.40 0.80 1.60 2.40 4.00 5.50 8.65 11.80 18.10 24.40 37.00 49.60 74.80 100.00"
    attach_fractions = [float(i) / 100 for i in attach_string.split()]
    initial_distance = 6.0
    pull_distances = np.arange(initial_distance, 14.5 + initial_distance, 0.5)
    guest_restraints = load_restraints(f"{windows_dir}/restraints.json")

    fe = analysis.fe_calc()
    fe.topology = f"{host_prefix}-{guest_prefix}-dum.prmtop"
    fe.trajectory = "production.nc"
    fe.path = str(windows_dir)
    fe.restraint_list = guest_restraints
    fe.collect_data()
    fe.methods = ["ti-block"]
    fe.ti_matrix = "full"
    fe.boot_cycles = boot_cycles
    fe.compute_free_energy()
    sim_data = fe.simulation_data

    fe.compute_ref_state_work([guest_restraints[0], guest_restraints[1], None, None, guest_restraints[2], None,])

    res = fe.results
    binding_affinity = -1 * (res["attach"]["ti-block"]["fe"] + res["pull"]["ti-block"]["fe"]
                             + res["ref_state_work"])
    sem = np.sqrt(res["attach"]["ti-block"]["sem"] ** 2 + res["pull"]["ti-block"]["sem"] ** 2)

    attach_fe = res["attach"]["ti-block"]["fe_matrix"][0, :].magnitude
    pull_fe = res["pull"]["ti-block"]["fe_matrix"][0, :].magnitude
    ref_work = res["ref_state_work"].magnitude
    dg, dg_sem = float(binding_affinity.magnitude), float(sem.magnitude)

    if plot_path:
        x_pull = pull_distances - 6 + attach_fractions[-1]
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
    return out, sim_data


def run(dp, host_prefix="BCD", guest_prefix="UNL"):
    dp = Path(dp)
    windows = dp / "windows"
    print(windows)

    # every window folder must have finished its production run
    missing = [w.name for w in sorted(windows.iterdir())
               if w.is_dir() and not (w / "production.nc").exists()]
    if missing:
        raise FileNotFoundError(f"{len(missing)} window(s) without production.nc, e.g. {missing[:3]}")

    results, sim_data = calculate_free_energy(windows, plot_path=dp / "binding_free_energy.jpg",
                                              grid_path=dp / "attach_windows_grid.jpg",
                                              title_name=dp.name, host_prefix=host_prefix,
                                              guest_prefix=guest_prefix)
    pull_window_histograms(sim_data, dp / "pull_window_histograms.jpg")
    (dp / "results.json").write_text(json.dumps(results, indent=2))
    return {"dG_bind": results["dG_bind"], "dG_sem": results["dG_sem"]}