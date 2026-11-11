#!/usr/bin/env python
"""Plot a custom atom-triplet angle distribution across selected APR windows."""
import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from paprika.io import load_restraints


DEFAULT_DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "PFBA"
DEFAULT_WINDOWS = ("a000", "a001", "a002", "a003", "a004", "a015")


def load_angle_restraint(data_dir):
    restraint_path = data_dir / "windows" / "restraints.json"
    restraints = load_restraints(str(restraint_path))
    angle_restraints = [
        restraint
        for restraint in restraints
        if restraint.restraint_type.value == "angle"
    ]
    if not angle_restraints:
        raise ValueError(f"No angle restraint found in {restraint_path}")
    return angle_restraints[0]


def load_window_angles(data_dir, window, restraint):
    window_dir = data_dir / "windows" / window
    trajectory_path = window_dir / "production.nc"
    topology_path = window_dir / "BCD-UNL-dum.prmtop"
    trajectory = md.load_netcdf(str(trajectory_path), top=str(topology_path))

    indices = [restraint.index1, restraint.index2, restraint.index3]
    if any(index is None or len(index) != 1 for index in indices):
        raise ValueError(
            "The first angle restraint must select exactly one atom for each "
            "of its three masks"
        )
    angle_indices = [
        int(index[0]) - int(restraint.amber_index) for index in indices
    ]
    if any(index < 0 or index >= trajectory.n_atoms for index in angle_indices):
        raise ValueError(
            f"Angle restraint atom indices are outside the topology in {window_dir}"
        )
    labels = [
        f"{trajectory.topology.atom(index).residue.name}:"
        f"{trajectory.topology.atom(index).name}"
        for index in angle_indices
    ]
    angles_deg = np.rad2deg(
        md.compute_angles(
            trajectory,
            np.asarray([angle_indices]),
            periodic=False,
        )[:, 0]
    )
    if not np.all(np.isfinite(angles_deg)):
        raise ValueError(f"Non-finite angle samples found in window {window}")

    return trajectory.time, angles_deg, labels


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Histogram the first Paprika angle restraint for selected "
            "attachment windows."
        )
    )
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--windows", nargs="+", default=DEFAULT_WINDOWS)
    parser.add_argument("--bins", type=int, default=36)
    args = parser.parse_args()

    if args.bins <= 0:
        parser.error("--bins must be a positive integer")
    if len(set(args.windows)) != len(args.windows):
        parser.error("--windows cannot contain duplicates")
    output_dir = args.output_dir or args.data_dir / "analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    restraint = load_angle_restraint(args.data_dir)
    labels = None
    targets = np.asarray(restraint.phase["attach"]["targets"].magnitude)
    samples = {}
    samples_path = None
    plot_path = None
    for window in args.windows:
        if not window.startswith("a"):
            parser.error(f"{window!r} is not an attachment window (expected aNNN)")
        window_index = int(window[1:])
        if window_index >= targets.size:
            parser.error(f"No attachment target found for window {window}")
        times, angles, window_labels = load_window_angles(
            args.data_dir, window, restraint
        )
        if labels is None:
            labels = window_labels
        elif labels != window_labels:
            raise ValueError(f"Angle restraint atoms differ in window {window}")
        samples[window] = (times, angles, float(targets[window_index]))

    stem = "restraint_angle_" + "_".join(label.replace(":", "-") for label in labels)
    samples_path = output_dir / f"{stem}_samples.csv"
    with samples_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["window", "frame", "time_ps", "angle_deg", "target_deg"])
        for window, (times, angles, target) in samples.items():
            writer.writerows(
                (window, frame, f"{time:.6f}", f"{angle:.6f}", f"{target:.6f}")
                for frame, (time, angle) in enumerate(zip(times, angles))
            )

    ncols = min(3, len(samples))
    nrows = (len(samples) + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4 * ncols, 3.5 * nrows), sharex=True, sharey=True
    )
    axes = np.atleast_1d(axes).ravel()
    for ax, (window, (_, angles, target)) in zip(axes, samples.items()):
        ax.hist(angles, bins=args.bins, range=(0, 180), color="steelblue")
        ax.set_title(f"{window} (n={len(angles)}, mean={np.mean(angles):.1f}°)")
        ax.axvline(target, color="red", linestyle="--", linewidth=1)
        ax.set_xlim(0, 185)
        ax.set_xlabel("Angle (degrees)")
        ax.set_ylabel("Frames")
    for ax in axes[len(samples) :]:
        ax.set_visible(False)
    fig.suptitle(
        "Paprika angle restraint: " + " - ".join(labels) + " by attachment window"
    )
    fig.tight_layout()
    plot_path = output_dir / f"{stem}_histogram.png"
    fig.savefig(plot_path, dpi=200)
    plt.close(fig)

    for window, (_, angles, target) in samples.items():
        print(
            f"{window}: {len(angles)} frames, "
            f"mean={np.mean(angles):.2f}°, std={np.std(angles):.2f}°, "
            f"target={target:.2f}°"
        )
    print(f"Histogram: {plot_path}")
    print(f"Samples: {samples_path}")


if __name__ == "__main__":
    main()
