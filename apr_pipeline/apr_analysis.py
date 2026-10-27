"""APR binding free energy and sampling diagnostics -> results.json + plots."""
import json
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")                 # must come before importing pyplot
import matplotlib.pyplot as plt

from paprika import analysis
from paprika.io import load_restraints
from paprika.restraints import create_window_list

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


def _read_distance_wall(disang_path):
    """Read the unique two-atom, four-bound distance wall from an Amber restraint file."""
    text = Path(disang_path).read_text()
    candidates = []
    for block in re.findall(r"&rst\b(.*?)&end", text, flags=re.IGNORECASE | re.DOTALL):
        atom_match = re.search(
            r"\biat\s*=\s*(.*?)(?=,\s*(?:igr[1-4]|r[1-4]|rk[23])\s*=)",
            block,
            flags=re.IGNORECASE | re.DOTALL,
        )
        if atom_match is None:
            continue
        atom_indices = [int(value) for value in re.findall(r"-?\d+", atom_match.group(1))]
        if len(atom_indices) != 2:
            continue

        parameters = {}
        for key in ("r1", "r2", "r3", "r4", "rk2", "rk3"):
            match = re.search(
                rf"\b{key}\s*=\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)",
                block,
                flags=re.IGNORECASE,
            )
            if match is None:
                break
            parameters[key] = float(match.group(1))
        if len(parameters) != 6:
            continue

        if (
            parameters["r1"] < parameters["r2"] < parameters["r3"] < parameters["r4"]
            and parameters["rk2"] > 0
            and parameters["rk3"] > 0
        ):
            if any(index <= 0 for index in atom_indices):
                raise ValueError(
                    f"Group-based distance walls are not supported in {disang_path}"
                )
            candidates.append({"atom_indices": atom_indices, **parameters})

    if len(candidates) > 1:
        raise ValueError(f"Multiple two-atom distance walls found in {disang_path}")
    return candidates[0] if candidates else None


def _amber_distance_wall_energy(distances, wall):
    """Return Amber's four-bound flat-bottom distance-restraint energy in kcal/mol."""
    distances = np.asarray(distances, dtype=float).reshape(-1)
    if distances.size == 0 or not np.all(np.isfinite(distances)):
        raise ValueError("Distance-wall samples must be non-empty and finite")

    r1, r2, r3, r4 = (wall[key] for key in ("r1", "r2", "r3", "r4"))
    rk2, rk3 = wall["rk2"], wall["rk3"]
    energy = np.zeros_like(distances)

    below_r1 = distances < r1
    between_r1_r2 = (distances >= r1) & (distances < r2)
    between_r3_r4 = (distances > r3) & (distances <= r4)
    above_r4 = distances > r4

    energy[below_r1] = (
        rk2 * (r2 - r1) ** 2
        + 2.0 * rk2 * (r2 - r1) * (r1 - distances[below_r1])
    )
    energy[between_r1_r2] = rk2 * (distances[between_r1_r2] - r2) ** 2
    energy[between_r3_r4] = rk3 * (distances[between_r3_r4] - r3) ** 2
    energy[above_r4] = (
        rk3 * (r4 - r3) ** 2
        + 2.0 * rk3 * (r4 - r3) * (distances[above_r4] - r4)
    )
    return energy


def _distance_wall_summary(windows_dir, simulation_data, window_order):
    """Summarize restraint energy and signed endpoint work for independent windows."""
    if len(simulation_data["attach"]) != len(window_order):
        raise ValueError(
            "Attachment window labels do not match the available simulation data"
        )

    attachment_windows = [
        (index, str(window))
        for index, window in enumerate(window_order)
        if str(window).startswith("a")
    ]
    if not attachment_windows:
        return None
    walls = [
        _read_distance_wall(Path(windows_dir) / window / "disang.rest")
        for _, window in attachment_windows
    ]
    if not any(wall is not None for wall in walls):
        return None
    if any(wall is None for wall in walls):
        raise ValueError("Distance wall is present in only some attachment windows")

    reference = walls[0]
    for (_, window), wall in zip(attachment_windows[1:], walls[1:]):
        if wall["atom_indices"] != reference["atom_indices"] or any(
            not np.isclose(wall[key], reference[key], rtol=0.0, atol=1e-5)
            for key in ("r1", "r2", "r3", "r4", "rk2", "rk3")
        ):
            raise ValueError(f"Distance-wall parameters differ in window {window}")

    per_window = {}
    all_energies = []
    all_active = []
    total_frames = 0
    for index, window in attachment_windows:
        distances = np.asarray(_mag(simulation_data["attach"][index][0]), dtype=float).reshape(-1)
        if distances.size == 0 or not np.all(np.isfinite(distances)):
            raise ValueError(f"Invalid distance-wall samples in window {window}")
        energies = _amber_distance_wall_energy(distances, reference)
        active = (distances < reference["r2"]) | (distances > reference["r3"])
        per_window[str(window)] = {
            "frame_count": int(distances.size),
            "mean_potential_energy_kcal_per_mol": float(np.mean(energies)),
            "maximum_potential_energy_kcal_per_mol": float(np.max(energies)),
            "active_frame_percent": float(100.0 * np.mean(active)),
            "net_work_by_wall_start_to_end_kcal_per_mol": float(
                energies[0] - energies[-1]
            ),
        }
        all_energies.append(energies)
        all_active.append(active)
        total_frames += distances.size

    pooled_energies = np.concatenate(all_energies)
    pooled_active = np.concatenate(all_active)
    return {
        "parameters": {
            "amber_atom_indices_1based": reference["atom_indices"],
            "r1_angstrom": reference["r1"],
            "r2_angstrom": reference["r2"],
            "r3_angstrom": reference["r3"],
            "r4_angstrom": reference["r4"],
            "rk2_kcal_per_mol_angstrom2": reference["rk2"],
            "rk3_kcal_per_mol_angstrom2": reference["rk3"],
        },
        "frame_count": int(total_frames),
        "mean_potential_energy_kcal_per_mol": float(np.mean(pooled_energies)),
        "maximum_potential_energy_kcal_per_mol": float(np.max(pooled_energies)),
        "active_frame_percent": float(100.0 * np.mean(pooled_active)),
        "per_window": per_window,
        "work_definition": (
            "For each independent window, net work done by the static wall from "
            "the first to last saved frame is U(first) - U(last). Per-window "
            "values are not summed because windows are independent trajectories."
        ),
    }


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

    fe.compute_ref_state_work(
        [
            guest_restraints[0],
            guest_restraints[1] if len(guest_restraints) > 1 else None,
            None,
            None,
            guest_restraints[2] if len(guest_restraints) > 2 else None,
            None,
        ]
    )

    res = fe.results
    all_window_names = create_window_list(guest_restraints)
    attach_window_names = [
        all_window_names[int(index)] for index in res["attach"]["window_order"]
    ]
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
    distance_wall = _distance_wall_summary(
        windows_dir,
        sim_data,
        attach_window_names,
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
    if distance_wall is not None:
        out["distance_wall"] = distance_wall
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