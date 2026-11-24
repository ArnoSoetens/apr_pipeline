"""APR binding free energy and sampling diagnostics -> results.json + plots."""
import json
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import matplotlib
matplotlib.use("Agg")                 # must come before importing pyplot
import matplotlib.pyplot as plt

from openff.units import unit as openff_unit
from paprika import analysis
from paprika.io import load_restraints, load_trajectory, read_restraint_data
from paprika.restraints import create_window_list

import parmed as pmd

from apr_config import CFG, get_attach_fractions, get_pull_distances, get_method
from apr_setup import _rim_plane_distance_from_dm1
from datapoints import PLOTS

# --- values from config.toml ---
COMPLEX = CFG["naming"]["complex_prefix"]
HOST = CFG["naming"]["host_prefix"]
GUEST = CFG["naming"]["guest_prefix"]
RIM_ATOM_GROUPS = CFG["apr"]["restraints"]["walls"]["distance"]["rim_atom_groups"]
INITIAL_DISTANCE = CFG["apr"]["initial_distance"]
METHOD = get_method()                     # "ti-block" or "mbar-autoc"
METHODS = ("ti-block", "mbar-autoc")
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


HALF_COLORS = ("tab:blue", "tab:orange")
HALF_LABELS = ("1st half", "2nd half")


def _halves(values):
    """First and second half of a window's frames (an odd extra frame goes to the second half)."""
    values = np.asarray(values)
    middle = len(values) // 2
    return values[:middle], values[middle:]


def _rim_heights_from_structure(structure, rim_atom_groups=RIM_ATOM_GROUPS, host_prefix=HOST):
    """(near, far) rim-plane heights in Å along +z from DM1, as apr_setup measures them."""
    dm1 = structure[":DM1"].coordinates[0]
    return tuple(sorted(
        _rim_plane_distance_from_dm1(structure, group, dm1, host_prefix)
        for group in rim_atom_groups
    ))


def rim_heights(windows_dir, complex_prefix=COMPLEX, rim_atom_groups=RIM_ATOM_GROUPS,
                host_prefix=HOST):
    """(near, far) rim heights from DM1, from the starting structure of windows/a000.

    The rim atoms are apr.restraints.walls.distance.rim_atom_groups in config.toml.
    """
    window = Path(windows_dir) / "a000"
    structure = pmd.load_file(str(window / f"{complex_prefix}-dum.prmtop"),
                              str(window / f"{complex_prefix}-dum.rst7"))
    return _rim_heights_from_structure(structure, rim_atom_groups, host_prefix)


def _grid_title(label):
    return f"{label} (λ=1)" if label == "p000" else label


def _window_grid(values, targets, labels, out_path, marks=(), mark_label="", title="",
                 xlabel="", ncols=5):
    """One histogram per window, 1st half of the frames blue and 2nd half orange, the target
    as a red dashed line and each x in marks as a dotted line. Returns the (closed) figure."""
    num_win = len(values)
    if len(labels) != num_win:
        raise ValueError(f"{len(labels)} labels for {num_win} attach windows")

    nrows = int(np.ceil(num_win / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 2.5 * nrows),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes).flatten()

    for k in range(num_win):
        ax = axes[k]
        window_values = np.asarray(values[k], dtype=float).reshape(-1)
        edges = np.histogram_bin_edges(window_values, bins=40)
        for part, color, name in zip(_halves(window_values), HALF_COLORS, HALF_LABELS):
            ax.hist(part, bins=edges, histtype="stepfilled", color=color, alpha=0.2)
            ax.hist(part, bins=edges, histtype="step", color=color, linewidth=1.2, label=name)
        ax.axvline(targets[k], color="red", linestyle="--", linewidth=1, label="target")
        for i, mark in enumerate(marks):
            ax.axvline(mark, color="0.3", linestyle=":", linewidth=1.2,
                       label=mark_label if i == 0 else None)
        ax.set_title(_grid_title(labels[k]), fontsize=9)
        ax.tick_params(labelsize=7)
        if k + ncols >= num_win:                     # no panel below: ticks and axis label
            ax.tick_params(labelbottom=True)
            ax.set_xlabel(xlabel, fontsize=8)

    for k in range(num_win, len(axes)):              # switch off unused slots
        axes[k].axis("off")

    handles, names = axes[0].get_legend_handles_labels()
    fig.suptitle(title, y=1.0)
    if num_win < len(axes):                          # legend in the first unused slot
        axes[num_win].legend(handles, names, loc="center", fontsize=9, frameon=False)
        fig.tight_layout()
    else:                                            # no free slot: one row above the panels
        fig.legend(handles, names, loc="upper center", ncol=len(names), fontsize=8,
                   frameon=False, bbox_to_anchor=(0.5, 0.985))
        fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return fig


def _attach_labels(labels, num_win):
    return labels if labels is not None else [f"a{k:03d}" for k in range(num_win)]


def attach_window_grid(fe, out_path, labels=None, rims=None, restraint_index=0, phase="attach",
                       ncols=5):
    """Distance restraint (DM1-G1, Å) per attach window: halves, target, rims as dotted lines.

    labels: window folder per panel (the λ=1 window is p000 with continuous_apr);
    a000, a001, ... when not given. Returns the (closed) figure.
    """
    data = fe.simulation_data[phase]                 # list per window, per restraint
    # ordered targets for this restraint, in the same window order as simulation_data
    targets = fe.prepare_data(phase)[5][restraint_index].magnitude
    values = [_mag(window[restraint_index]) for window in data]
    return _window_grid(values, targets, _attach_labels(labels, len(data)), out_path,
                        marks=tuple(rims or ()), mark_label="rim",
                        title="Attach phase: DM1-G1 distance per window",
                        xlabel="DM1-G1 distance (Å)", ncols=ncols)


def _degrees(angles):
    """Angles in degrees: pint quantities are converted, plain numbers are taken as radians
    (pAPRika's fe_calc reads angles in radians)."""
    if hasattr(angles, "to"):
        return np.asarray(angles.to("degree").magnitude, dtype=float)
    return np.degrees(np.asarray(angles, dtype=float))


def _first_angle_index(active_restraints):
    """Position of the first angle restraint (three masks) among the active restraints."""
    for index, restraint in enumerate(active_restraints):
        if restraint.mask3 is not None and restraint.mask4 is None:
            return index
    raise ValueError("no angle restraint among the attach restraints")


def attach_angle_grid(fe, out_path, labels=None, cone=None, phase="attach", ncols=5):
    """First guest angle restraint (DM2-DM1-G1, degrees) per attach window: halves, target,
    and the cone wall (r2, where its force starts) as a dotted line. Returns the figure."""
    data = fe.simulation_data[phase]
    prepared = fe.prepare_data(phase)
    index = _first_angle_index(prepared[3])
    targets = _degrees(prepared[5][index])
    values = [_degrees(window[index]) for window in data]
    return _window_grid(values, targets, _attach_labels(labels, len(data)), out_path,
                        marks=() if cone is None else (cone,), mark_label="cone wall",
                        title="Attach phase: DM2-DM1-G1 angle per window",
                        xlabel="DM2-DM1-G1 angle (degrees)", ncols=ncols)


def _read_cone_wall(disang_path):
    """r2 (degrees) of the cone wall in an Amber restraint file, or None when there is none.

    The cone is the three-atom flat-bottom block with r1 < r2 < r3 and rk2 > 0; the guest
    angle restraints always have r2 == r3.
    """
    text = Path(disang_path).read_text()
    for block in re.findall(r"&rst\b(.*?)&end", text, flags=re.IGNORECASE | re.DOTALL):
        atoms = re.search(r"\biat\s*=\s*(.*?)(?=,\s*r1\s*=)", block, flags=re.IGNORECASE | re.DOTALL)
        if atoms is None or len(re.findall(r"-?\d+", atoms.group(1))) != 3:
            continue
        values = {}
        for key in ("r1", "r2", "r3", "rk2"):
            match = re.search(rf"\b{key}\s*=\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)",
                              block, flags=re.IGNORECASE)
            if match is None:
                break
            values[key] = float(match.group(1))
        if len(values) == 4 and values["r1"] < values["r2"] < values["r3"] and values["rk2"] > 0:
            return values["r2"]
    return None


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


def _read_distance_walls(disang_path):
    """Read all two-atom, four-bound distance walls from an Amber restraint file."""
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
            parameters["r1"] <= parameters["r2"] < parameters["r3"] < parameters["r4"]
            and parameters["rk2"] > 0
            and parameters["rk3"] > 0
        ):
            if any(index <= 0 for index in atom_indices):
                raise ValueError(
                    f"Group-based distance walls are not supported in {disang_path}"
                )
            candidates.append({"atom_indices": atom_indices, **parameters})

    return candidates


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


def _distance_wall_summary(windows_dir, simulation_data, window_order, topology=None):
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
    walls_by_window = [
        _read_distance_walls(Path(windows_dir) / window / "disang.rest")
        for _, window in attachment_windows
    ]
    if not any(walls_by_window):
        return None
    if any(len(walls) != len(walls_by_window[0]) for walls in walls_by_window):
        raise ValueError("Distance-wall count differs between attachment windows")

    reference_walls = walls_by_window[0]
    for (_, window), walls in zip(attachment_windows[1:], walls_by_window[1:]):
        for wall_index, (reference, wall) in enumerate(zip(reference_walls, walls)):
            if wall["atom_indices"] != reference["atom_indices"] or any(
                not np.isclose(wall[key], reference[key], rtol=0.0, atol=1e-5)
                for key in ("r1", "r2", "r3", "r4", "rk2", "rk3")
            ):
                raise ValueError(
                    f"Distance-wall parameters differ in window {window} "
                    f"for wall {wall_index + 1}"
                )

    summaries = []
    for wall_index, reference in enumerate(reference_walls):
        per_window = {}
        all_energies = []
        all_active = []
        total_frames = 0
        for simulation_index, window in attachment_windows:
            if wall_index == 0:
                # The first generated wall is DM1-G1, the APR distance restraint.
                distances = _mag(simulation_data["attach"][simulation_index][0])
            else:
                if topology is None:
                    raise ValueError(
                        "A topology is required to analyze additional distance walls"
                    )
                trajectory = load_trajectory(
                    str(Path(windows_dir) / window), TRAJECTORY, topology
                )
                atom1, atom2 = reference["atom_indices"]
                distance_restraint = SimpleNamespace(
                    mask1=f"@{atom1}",
                    mask2=f"@{atom2}",
                    mask3=None,
                    mask4=None,
                )
                distances = _mag(read_restraint_data(trajectory, distance_restraint))

            distances = np.asarray(distances, dtype=float).reshape(-1)
            if distances.size == 0 or not np.all(np.isfinite(distances)):
                raise ValueError(f"Invalid distance-wall samples in window {window}")
            energies = _amber_distance_wall_energy(distances, reference)
            active = energies > 0.0
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
        summaries.append({
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
        })

    if len(summaries) == 1:
        result = summaries[0]
    else:
        result = {"walls": summaries}
    result["work_definition"] = (
        "For each independent window, net work done by the static wall from "
        "the first to last saved frame is U(first) - U(last). Per-window "
        "values are not summed because windows are independent trajectories."
    )
    return result


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


def _ti_per_frame(prepared, phase):
    """dU per frame and the integration variable per window, as pAPRika's run_ti builds them.

    attach/release: dU = sum_r k_max,r (x_r - x0_r)^2 (kcal/mol), x-axis lambda = k / k_max of
    the first restraint; pull: dU = sum_r 2 k_max,r (x_r - x0_r) (kcal/mol/Å), x-axis the
    target of the first restraint. Units are those prepare_data converted to (Å, radians).
    Dihedrals are wrapped around their target like pAPRika does. Nothing is modified in place.
    Returns (x per window, list of dU arrays).
    """
    num_win, _, _, active, force_constants, targets, values = prepared
    fcs = [np.asarray(f.magnitude, dtype=float) for f in force_constants]
    tgs = [np.asarray(t.magnitude, dtype=float) for t in targets]
    max_fcs = [float(np.max(f)) for f in fcs]
    period = float((360.0 * openff_unit.degree).to(openff_unit.radian).magnitude)

    dl_vals, dU = np.zeros(num_win), []
    for k in range(num_win):
        total = 0.0
        for r, restraint in enumerate(active):
            x = values[k][r]
            if hasattr(x, "to"):
                x = x.to(targets[r].units)
            x = np.array(getattr(x, "magnitude", x), dtype=float)    # a copy
            target = tgs[r][k]
            if restraint.mask3 is not None and restraint.mask4 is not None:
                x[x < target - period / 2] += period
                x[x > target + period / 2] -= period
            if phase == "pull":
                total = total + 2.0 * max_fcs[r] * (x - target)
            else:
                total = total + max_fcs[r] * (x - target) ** 2
        dU.append(np.asarray(total, dtype=float))
        dl_vals[k] = tgs[0][k] if phase == "pull" else fcs[0][k] / max_fcs[0]
    return dl_vals, dU


def _ti_integral(dl_vals, means, sems, rng, cycles=BOOT_CYCLES):
    """pAPRika's TI-block integral and SEM from window means and SEMs: bootstrap normal draws,
    an Akima spline with 100 points per interval, and the [0, -1] element of its "full" matrix
    (paprika.analysis.bootstrap.integrate_bootstraps)."""
    from scipy.interpolate import Akima1DInterpolator

    x_intp = np.concatenate(
        [np.linspace(dl_vals[k - 1], dl_vals[k], num=100, endpoint=False)
         for k in range(1, len(dl_vals))] + [[dl_vals[-1]]])
    samples = rng.normal(means, sems, size=(cycles, len(means)))
    end = len(x_intp) - 1                       # integrate_bootstraps stops before the last point
    integrals = np.array([np.trapz(Akima1DInterpolator(dl_vals, y)(x_intp)[:end], x_intp[:end])
                          for y in samples])
    return float(np.mean(integrals)), float(np.std(integrals))


def _block_sem(values):
    """pAPRika's blocking SEM: max over block sizes, on the get_nearest_max(N) frames."""
    from paprika.analysis.utils import get_block_sem, get_nearest_max

    return float(get_block_sem(values[:get_nearest_max(len(values))]))


def _convergence_curves(dl_vals, dU, phase, fractions=CONVERGENCE_FRACTIONS,
                        cycles=BOOT_CYCLES, seed=0):
    """TI estimate, SEM and total frames for the first ("forward") and last ("reverse")
    fraction f of every window. The SEM is a new blocking analysis of each piece (pAPRika only
    rescales the full-data SEM). The pull work gets pAPRika's minus sign."""
    sign = -1.0 if phase == "pull" else 1.0
    curves = {"fractions": np.asarray(fractions, dtype=float)}
    for direction in ("forward", "reverse"):
        fe, sem, frames = [], [], []
        for i, fraction in enumerate(fractions):
            parts = []
            for window in dU:
                m = int(fraction * len(window))
                # reverse: the last m frames, newest first, so blocking drops the oldest ones
                parts.append(window[:m] if direction == "forward" else window[::-1][:m])
            means = np.array([np.mean(part) for part in parts])
            sems = np.array([_block_sem(part) for part in parts])
            value, error = _ti_integral(dl_vals, means, sems,
                                        np.random.default_rng([seed, i]), cycles)
            fe.append(sign * value)
            sem.append(error)
            frames.append(sum(len(part) for part in parts))
        curves[direction] = {"fe": np.array(fe), "sem": np.array(sem), "frames": np.array(frames)}
    return curves


def convergence_plot(fe, out_path, phases=("attach", "pull"), fractions=CONVERGENCE_FRACTIONS,
                     cycles=BOOT_CYCLES):
    """TI-block free energy from the first (forward) and last (reverse) fraction of every window,
    and the SEM against the number of frames on log-log axes with an N^-1/2 line.

    One row per phase. Returns (figure, {phase: curves}).
    """
    styles = {"forward": dict(marker="o", color="tab:blue", label="forward (first f)"),
              "reverse": dict(marker="s", color="tab:orange", label="reverse (last f)")}
    fig, axes = plt.subplots(len(phases), 2, figsize=(12, 4.5 * len(phases)), squeeze=False)
    all_curves = {}
    for row, phase in zip(axes, phases):
        dl_vals, dU = _ti_per_frame(fe.prepare_data(phase), phase)
        curves = all_curves[phase] = _convergence_curves(dl_vals, dU, phase, fractions, cycles)
        percent = curves["fractions"] * 100
        for direction, offset in (("forward", 0.0), ("reverse", 0.8)):
            c = curves[direction]
            row[0].errorbar(percent + offset, c["fe"], yerr=c["sem"], capsize=2, **styles[direction])
            row[1].loglog(c["frames"], c["sem"], **styles[direction])
        full = curves["forward"]
        reference = full["sem"][-1] * np.sqrt(full["frames"][-1] / full["frames"])
        row[1].loglog(full["frames"], reference, "k--", linewidth=1, label=r"$\propto N^{-1/2}$")

        row[0].set_title(f"{phase.capitalize()} phase (TI-block)")
        row[0].set_xlabel("Trajectory frames used per window (%)")
        row[0].set_ylabel("Free energy ± 1 SEM (kcal/mol)")
        row[0].set_xticks(percent)
        row[1].set_title(f"{phase.capitalize()} phase: SEM vs frames")
        row[1].set_xlabel("Frames used (all windows)")
        row[1].set_ylabel("SEM (kcal/mol)")
        for ax in row:
            ax.grid(alpha=0.3, which="both")
            ax.legend(fontsize=8)
    fig.suptitle("Free-energy convergence: forward and reverse")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    return fig, all_curves


def _reduced_potentials(prepared, beta):
    """u_kln[k, l, n] = beta * sum_r k_r,l (x_r,k,n - x0_r,l)^2: frame n of window k in the
    potential of window l, as pAPRika's run_mbar builds it (dihedrals wrapped around the target
    of window l; on copies). beta in mol/kcal. Returns (u_kln, N_k)."""
    num_win, data_points, max_points, active, force_constants, targets, values = prepared
    fcs = [np.asarray(f.magnitude, dtype=float) for f in force_constants]
    tgs = [np.asarray(t.magnitude, dtype=float) for t in targets]
    period = float((360.0 * openff_unit.degree).to(openff_unit.radian).magnitude)
    N_k = np.asarray(data_points, dtype=int)
    u_kln = np.zeros((num_win, num_win, max_points))
    for k in range(num_win):
        for r, restraint in enumerate(active):
            x = values[k][r]
            if hasattr(x, "to"):
                x = x.to(targets[r].units)
            x = np.asarray(getattr(x, "magnitude", x), dtype=float)
            deviation = x[None, :] - tgs[r][:, None]                 # (window l, frame n)
            if restraint.mask3 is not None and restraint.mask4 is not None:
                deviation = (deviation + period / 2) % period - period / 2
            u_kln[k, :, :N_k[k]] += beta * fcs[r][:, None] * deviation ** 2
    return u_kln, N_k


def _min_neighbour_overlap(u_kln, N_k):
    """Smallest MBAR overlap-matrix element between neighbouring windows (O[k, k+1] and
    O[k+1, k]); all frames are used. Rule of thumb: below about 0.03 the pair overlaps poorly."""
    import pymbar

    overlap = pymbar.MBAR(u_kln, N_k).compute_overlap()["matrix"]
    pairs = [overlap[k, k + 1] for k in range(len(N_k) - 1)]
    pairs += [overlap[k + 1, k] for k in range(len(N_k) - 1)]
    return float(min(pairs))


def phase_overlaps(fe, phases=("attach", "pull")):
    """{phase: smallest neighbour overlap} for results.json."""
    beta = float(fe.beta.to("mol/kcal").magnitude)
    return {phase: _min_neighbour_overlap(*_reduced_potentials(fe.prepare_data(phase), beta))
            for phase in phases}


def monotonicity_check(dl_vals, dU, alpha=0.05):
    """<R>_lambda must not increase with lambda (d<R>/dlambda = -beta Var(R) <= 0, so the attach
    work W(lambda) is concave). R is the attach dU per frame (full force constants).

    Windows are sorted by lambda; a rise between neighbours is a violation when
    z = (m_j - m_i) / sqrt(s_i^2 + s_j^2) > z_crit, with s the blocking SEM and z_crit one-sided,
    Bonferroni-corrected over the K-1 pairs: Phi^-1(1 - alpha / (K-1)).
    """
    from statistics import NormalDist

    order = np.argsort(dl_vals, kind="stable")
    lam = np.asarray(dl_vals, dtype=float)[order]
    means = np.array([np.mean(dU[i]) for i in order])
    sems = np.array([_block_sem(np.asarray(dU[i])) for i in order])
    pairs = len(lam) - 1
    z_crit = NormalDist().inv_cdf(1.0 - alpha / max(pairs, 1))
    rise = means[1:] - means[:-1]
    combined = np.sqrt(sems[1:] ** 2 + sems[:-1] ** 2)
    z = rise / np.maximum(combined, np.finfo(float).tiny)
    violations = [{"lambda_i": float(lam[k]), "lambda_j": float(lam[k + 1]),
                   "rise": float(rise[k]), "z": float(z[k])}
                  for k in range(pairs) if z[k] > z_crit]
    return {"monotonic": not violations, "z_crit": float(z_crit),
            "max_z": float(np.max(z)) if pairs else None, "violations": violations}


ROTAMER_CENTERS = (60.0, 180.0, 300.0)      # g+, t, g- (degrees)
ROTAMER_CORE = 30.0                         # a frame belongs to a state within +-30 deg of it
MAX_TRANSITIONS_COLOUR = 10                 # colour scale of the transition heatmap: 0-10


def _guest_torsions(topology, resname=GUEST):
    """Rotatable guest torsions: RDKit's non-ring torsion list (TorsionFingerprints) of the guest
    residue, built from the bonds of the simulated topology (MDTraj) with every bond single, so
    no atom types or bond orders are needed. Returns [(0-based atom indices, "A-B-C-D")]."""
    from rdkit import Chem
    from rdkit.Chem import TorsionFingerprints

    residue_atoms = [a for a in topology.atoms if a.residue.name == resname]
    if not residue_atoms:
        raise ValueError(f"no residue named {resname} in the topology")
    if len({a.residue.index for a in residue_atoms}) > 1:
        raise ValueError(f"more than one residue named {resname} in the topology")
    local = {atom.index: i for i, atom in enumerate(residue_atoms)}

    mol = Chem.RWMol()
    for atom in residue_atoms:
        rd_atom = Chem.Atom(int(atom.element.atomic_number))
        rd_atom.SetNoImplicit(True)
        mol.AddAtom(rd_atom)
    for a, b in topology.bonds:
        if a.index in local and b.index in local:
            mol.AddBond(local[a.index], local[b.index], Chem.BondType.SINGLE)
    mol = mol.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(mol)

    nonring, _ = TorsionFingerprints.CalculateTorsionLists(mol)
    torsions = []
    for atom_sets, _ in nonring:
        atoms = [residue_atoms[i] for i in atom_sets[0]]
        torsions.append(([a.index for a in atoms], "-".join(a.name for a in atoms)))
    return sorted(torsions, key=lambda t: sorted(t[0][1:3]))     # by the central bond


def count_transitions(theta_deg, centers=ROTAMER_CENTERS, core=ROTAMER_CORE):
    """Transitions between rotamer states, counted with core regions: a frame is assigned to a
    state only within +-core of its centre; between cores it keeps its last state, so jitter
    across a boundary is not counted."""
    theta = np.mod(np.asarray(theta_deg, dtype=float), 360.0)
    c = np.asarray(centers, dtype=float)
    distance = np.abs((theta[:, None] - c[None, :] + 180.0) % 360.0 - 180.0)
    state = np.where((distance <= core).any(axis=1), np.argmin(distance, axis=1), -1)
    assigned = state[state >= 0]
    return int(np.count_nonzero(assigned[1:] != assigned[:-1]))


def _js_divergence(a, b, bins):
    """Jensen-Shannon divergence in bits: 0 = identical histograms, 1 = no overlap."""
    p = np.histogram(a, bins=bins)[0].astype(float)
    q = np.histogram(b, bins=bins)[0].astype(float)
    p, q = p / max(p.sum(), 1.0), q / max(q.sum(), 1.0)
    m = 0.5 * (p + q)
    with np.errstate(divide="ignore", invalid="ignore"):
        kl_p = np.where(p > 0, p * np.log2(p / m), 0.0).sum()
        kl_q = np.where(q > 0, q * np.log2(q / m), 0.0).sum()
    return float(0.5 * kl_p + 0.5 * kl_q)


def _js_halves(theta_deg, bin_width=10.0):
    """JS divergence between the 1st and 2nd half of a dihedral time series (0-360 deg bins)."""
    theta = np.mod(np.asarray(theta_deg, dtype=float), 360.0)
    return _js_divergence(*_halves(theta), np.arange(0.0, 360.0 + bin_width, bin_width))


def plot_dihedral_overview(windows, names, transitions, js, out_path):
    """Two heatmaps, windows (rows) x dihedrals (columns): transitions (colour 0-10, >= 10 shown
    as 10, the count written in each cell) and the JS divergence of the 1st vs 2nd half.
    Returns the (closed) figure."""
    transitions = np.asarray(transitions)
    js = np.asarray(js, dtype=float)
    width = max(10.0, 2 * (0.5 * len(names) + 2.0))
    fig, axes = plt.subplots(1, 2, figsize=(width, 0.28 * len(windows) + 2.5), sharey=True)
    image = axes[0].imshow(transitions, aspect="auto", cmap="viridis", vmin=0,
                           vmax=MAX_TRANSITIONS_COLOUR, interpolation="nearest")
    for i in range(transitions.shape[0]):
        for j in range(transitions.shape[1]):
            colour = "black" if transitions[i, j] >= 0.6 * MAX_TRANSITIONS_COLOUR else "white"
            axes[0].text(j, i, str(int(transitions[i, j])), ha="center", va="center",
                         fontsize=6, color=colour)
    bar = fig.colorbar(image, ax=axes[0], extend="max", fraction=0.046)
    bar.set_label(f"transitions (≥{MAX_TRANSITIONS_COLOUR} shown as {MAX_TRANSITIONS_COLOUR})")
    axes[0].set_title("Rotamer transitions per window")

    js_image = axes[1].imshow(js, aspect="auto", cmap="magma_r", vmin=0.0,
                              vmax=max(0.2, float(js.max()) if js.size else 0.2),
                              interpolation="nearest")
    fig.colorbar(js_image, ax=axes[1], fraction=0.046).set_label("JS divergence (bits)")
    axes[1].set_title("1st vs 2nd half (JS divergence)")

    for ax in axes:
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=90, fontsize=7)
    axes[0].set_yticks(range(len(windows)))
    axes[0].set_yticklabels(windows, fontsize=7)
    fig.suptitle("Guest dihedrals: g+/t/g- states (cores ±30° around 60/180/300°)")
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return fig


def dihedral_overview(windows_dir, out_path, complex_prefix=COMPLEX, trajectory=TRAJECTORY):
    """Transitions and half-vs-half JS divergence of every rotatable guest torsion in every
    window folder (attach and pull, in folder order) -> plot. Returns the figure, or None when
    the guest has no rotatable torsions."""
    import mdtraj as md

    windows_dir = Path(windows_dir)
    windows = sorted(w.name for w in windows_dir.iterdir()
                     if w.is_dir() and (w / trajectory).is_file())
    prmtop = windows_dir / windows[0] / f"{complex_prefix}-dum.prmtop"
    topology = md.load_prmtop(str(prmtop))
    torsions = _guest_torsions(topology)
    if not torsions:
        print("no rotatable guest torsions: no dihedral plot", flush=True)
        return None
    needed = sorted({i for quad, _ in torsions for i in quad})
    local = np.array([[needed.index(i) for i in quad] for quad, _ in torsions])

    transitions = np.zeros((len(windows), len(torsions)), dtype=int)
    js = np.zeros((len(windows), len(torsions)))
    for row, window in enumerate(windows):
        traj = md.load(str(windows_dir / window / trajectory), top=str(prmtop),
                       atom_indices=needed)
        angles = np.degrees(md.compute_dihedrals(traj, local, periodic=False))
        for col in range(len(torsions)):
            transitions[row, col] = count_transitions(angles[:, col])
            js[row, col] = _js_halves(angles[:, col])
    return plot_dihedral_overview(windows, [name for _, name in torsions], transitions, js,
                                  out_path)


def calculate_free_energy(windows_dir, plot_path=None, grid_path=None, title_name="guest",
                          complex_prefix=COMPLEX, boot_cycles=BOOT_CYCLES, enhanced=False,
                          convergence_path=None, angle_grid_path=None, dihedral_path=None):
    """Returns (results dict, simulation_data). dG values in kcal/mol.

    The free energy at several fractions of the trajectory (CONVERGENCE_FRACTIONS) is only
    computed for the convergence plot (convergence_path); otherwise just the full data set.
    """

    attach_fractions = get_attach_fractions()
    initial_distance = INITIAL_DISTANCE
    pull_distances = get_pull_distances()
    guest_restraints = load_restraints(f"{windows_dir}/restraints.json")

    fe = analysis.fe_calc()
    topology = f"{complex_prefix}-dum.prmtop"
    fe.topology = topology
    fe.trajectory = TRAJECTORY
    fe.path = str(windows_dir)
    fe.restraint_list = guest_restraints
    fe.collect_data()
    fe.methods = list(METHODS)
    fe.fractions = CONVERGENCE_FRACTIONS.tolist() if convergence_path else [1.0]
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
        convergence_plot(fe, convergence_path)

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
        topology=topology,
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
        try:
            rims = rim_heights(windows_dir, complex_prefix)
        except Exception as e:                       # a plot detail must not block results.json
            rims = None
            print(f"warning: no rim lines in the attach-window grid: {e}", flush=True)
        attach_window_grid(fe, grid_path, labels=attach_window_names, rims=rims)

    if angle_grid_path:
        try:
            cone = _read_cone_wall(Path(windows_dir) / "a000" / "disang.rest")
        except OSError as e:                         # a plot detail must not block results.json
            cone = None
            print(f"warning: no cone line in the attach angle grid: {e}", flush=True)
        attach_angle_grid(fe, angle_grid_path, labels=attach_window_names, cone=cone)

    if dihedral_path:
        try:
            dihedral_overview(windows_dir, dihedral_path, complex_prefix)
        except Exception as e:                       # a plot detail must not block results.json
            print(f"warning: no dihedral plot: {e}", flush=True)

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
        out["overlap"] = phase_overlaps(fe)
        out["monotonicity_attach"] = monotonicity_check(
            *_ti_per_frame(fe.prepare_data("attach"), "attach"))
    return out, sim_data


def run(dp, complex_prefix=COMPLEX, enhanced=False):
    """Analyze one fully simulated data point: <dp>/results.json and plots in <dp>/plots/.

    Always binding_free_energy.jpg; enhanced=True also the convergence plot, the attach-window
    distance and angle grids, the guest-dihedral heatmaps, the pull-window histograms and extra
    metadata in results.json.
    """
    dp = Path(dp)
    windows = dp / "windows"

    # every window folder must have finished its production run
    missing = [w.name for w in sorted(windows.iterdir())
               if w.is_dir() and not (w / TRAJECTORY).exists()]
    if missing:
        raise FileNotFoundError(f"{len(missing)} window(s) without {TRAJECTORY}, e.g. {missing[:3]}")

    plots = dp / PLOTS
    plots.mkdir(exist_ok=True)
    results, sim_data = calculate_free_energy(
        windows,
        plot_path=plots / "binding_free_energy.jpg",
        grid_path=(plots / "attach_windows_grid.jpg") if enhanced else None,
        title_name=dp.name, complex_prefix=complex_prefix, enhanced=enhanced,
        convergence_path=(plots / "convergence_plot.jpg") if enhanced else None,
        angle_grid_path=(plots / "attach_angle_grid.jpg") if enhanced else None,
        dihedral_path=(plots / "dihedral_overview.jpg") if enhanced else None,
    )
    if enhanced:
        pull_window_histograms(sim_data, plots / "pull_window_histograms.jpg")
    (dp / "results.json").write_text(json.dumps(results, indent=2))
    return results
