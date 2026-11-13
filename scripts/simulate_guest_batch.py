#!/usr/bin/env python3
"""Run guest-only implicit-solvent simulations across data points in parallel."""

import argparse
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "apr_pipeline"))

from apr_config import CFG

SIMULATOR = SCRIPT_DIR / "simulate_guest.py"
APR_SIMULATOR = SCRIPT_DIR / "simulate_guest_apr.py"
GUEST_PREFIX = CFG["naming"]["guest_prefix"]
DATA_DIR = Path(CFG["directories"]["data"])
SIMULATIONS = {
    "guest": (SIMULATOR, "guest_implicit"),
    "apr-p000": (APR_SIMULATOR, "guest_apr_p000"),
}


def _positive_int(value):
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def _get_cpu_count(requested):
    if requested is not None:
        return requested
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus:
        return _positive_int(slurm_cpus)
    return os.cpu_count() or 1


def _data_points(data_dir):
    if not data_dir.is_dir():
        raise FileNotFoundError(f"Data directory does not exist: {data_dir}")
    return sorted(
        path
        for path in data_dir.iterdir()
        if path.is_dir()
    )


def _guest_files(data_point):
    return {
        "pdb": data_point / f"{GUEST_PREFIX}.pdb",
        "mol2": data_point / f"{GUEST_PREFIX}.mol2",
        "frcmod": data_point / f"{GUEST_PREFIX}.frcmod",
    }


def _clean_output_dir(output_dir):
    retained_topologies = [
        path for path in output_dir.iterdir()
        if path.is_file() and path.suffix == ".prmtop"
    ]
    production_trajectory = output_dir / "production.nc"
    if not retained_topologies or not production_trajectory.is_file():
        raise FileNotFoundError(
            f"Cannot clean {output_dir}: expected at least one .prmtop "
            "and production.nc"
        )

    for path in output_dir.iterdir():
        if path.is_file() and (
            path.name == "production.nc" or path.suffix == ".prmtop"
        ):
            continue
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        else:
            path.unlink()


def _run_one(data_point, files, simulator, output_name, pmemd_command, clean=False):
    output_dir = data_point / output_name
    log_path = output_dir / "batch.log"
    command = [
        sys.executable,
        str(simulator),
        str(files["pdb"]),
        str(files["mol2"]),
        "--frcmod",
        str(files["frcmod"]),
        "--output-dir",
        str(output_dir),
        "--pmemd-command",
        pmemd_command,
    ]
    env = os.environ.copy()
    env.update(
        {
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
        }
    )
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
            env=env,
        )
        log_content = result.stdout + result.stderr
        log_path.write_text(log_content, encoding="utf-8")
        if clean and result.returncode == 0:
            try:
                _clean_output_dir(output_dir)
            except OSError as exc:
                log_path.write_text(
                    f"{log_content}\nCleaning failed: {exc}\n", encoding="utf-8"
                )
                return data_point.name, 1, log_path, str(exc)
    except OSError as exc:
        return data_point.name, 1, log_path, str(exc)
    return data_point.name, result.returncode, log_path, ""


def run_batch(
    data_dir=DATA_DIR,
    cpus=None,
    simulator=None,
    simulation="guest",
    clean=False,
):
    data_dir = Path(data_dir).expanduser().resolve()
    worker_count = _get_cpu_count(cpus)
    points = _data_points(data_dir)
    if simulation not in SIMULATIONS:
        raise ValueError(f"Unknown simulation type: {simulation}")
    default_simulator, output_name = SIMULATIONS[simulation]
    simulator = simulator or default_simulator
    jobs = []
    skipped = 0

    for data_point in points:
        files = _guest_files(data_point)
        missing = [str(path) for path in files.values() if not path.is_file()]
        if missing:
            skipped += 1
            print(f"[skip] {data_point.name}: missing {', '.join(missing)}", flush=True)
        else:
            jobs.append((data_point, files))

    if not jobs:
        print(f"No data points have all required guest files under {data_dir}.")
        return 0

    print(
        f"Running {len(jobs)} guest simulation(s) on up to {worker_count} CPU(s).",
        flush=True,
    )
    failures = []
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = {
            executor.submit(
                _run_one, point, files, simulator, output_name, "pmemd", clean
            ): point
            for point, files in jobs
        }
        for future in as_completed(futures):
            point_name, return_code, log_path, error = future.result()
            if return_code:
                failures.append(point_name)
                detail = f"see {log_path}" if not error else error
                print(f"[failed] {point_name}: {detail}", flush=True)
            else:
                print(f"[done] {point_name}", flush=True)

    print(
        f"Batch complete: {len(jobs) - len(failures)} succeeded, "
        f"{len(failures)} failed, {skipped} skipped.",
        flush=True,
    )
    if failures:
        print(f"Failed data points: {', '.join(sorted(failures))}", flush=True)
        return 1
    return 0


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Run a guest-only GB simulation for every data point "
            "with guest PDB, MOL2, and FRCMOD files."
        )
    )
    parser.add_argument(
        "--simulation",
        choices=sorted(SIMULATIONS),
        default="guest",
        help="Simulation workflow to run (default: guest)",
    )
    parser.add_argument(
        "--cpus-per-task",
        type=_positive_int,
        help="Maximum concurrent data points (default: SLURM_CPUS_PER_TASK or CPU count)",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DATA_DIR,
        help=f"Data-point directory (default: {DATA_DIR})",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="after successful simulations, remove output files except .prmtop and production.nc",
    )
    args = parser.parse_args()
    return run_batch(
        args.data_dir,
        args.cpus_per_task,
        simulation=args.simulation,
        clean=args.clean,
    )


if __name__ == "__main__":
    raise SystemExit(main())
