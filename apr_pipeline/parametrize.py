import argparse
from pathlib import Path
import subprocess
import sys

from apr_config import CFG
from guest_input import ensure_guest_input, net_charge


def parametrize(
    data_dir,
    guest_pdb=f'{CFG["naming"]["guest_prefix"]}.pdb',
    guest_prefix=CFG["naming"]["guest_prefix"],
):
    """antechamber + parmchk2 on the guest PDB; the charge comes from guest_input.net_charge()."""
    charge = net_charge(data_dir)

    result = subprocess.run(
    [
        "antechamber",
        "-i", f"{guest_pdb}",
        "-fi", "pdb",
        "-o", f"{guest_prefix}.mol2",
        "-fo", "mol2",
        "-c", CFG["parametrize"]["charge_method"],
        "-pf", "y",
        "-nc", str(charge),
        "-rn", guest_prefix,
        "-at", CFG["parametrize"]["atom_type"]
    ],
    cwd=data_dir,
    check=True,
    capture_output=True,
    text=True,
    )

    result2 = subprocess.run(
        [
            "parmchk2",
            "-i", f"{guest_prefix}.mol2",
            "-f", "mol2",
            "-o", f"{guest_prefix}.frcmod",
            "-s", str(CFG["parametrize"]["parmchk_set"])

        ],
        cwd=data_dir,
        check=True,
        capture_output=True,
        text=True,
    )
    return


def run(dp):
    """Parametrize stage: <guest>.pdb -> <guest>.mol2 + <guest>.frcmod (charge: guest_input.net_charge())."""
    dp = Path(dp)
    ensure_guest_input(dp)
    guest_pdb = dp / f'{CFG["naming"]["guest_prefix"]}.pdb'
    if not guest_pdb.is_file():
        raise FileNotFoundError(f"guest PDB not found: {guest_pdb}")
    parametrize(dp)


def parametrize_data_points(data_dir=CFG["directories"]["data"], names=None):
    """Parametrize data-point folders containing the configured guest PDB."""
    data_dir = Path(data_dir)
    if not data_dir.is_dir():
        raise FileNotFoundError(f"data directory does not exist: {data_dir}")

    guest_pdb = f'{CFG["naming"]["guest_prefix"]}.pdb'
    if names:
        directories = []
        for name in names:
            if Path(name).name != name or name in (".", ".."):
                raise ValueError(f"data-point name must be a directory name: {name!r}")
            directory = data_dir / name
            if not directory.is_dir():
                raise FileNotFoundError(f"data-point directory does not exist: {directory}")
            directories.append(directory)
    else:
        directories = sorted(
            directory for directory in data_dir.iterdir()
            if directory.is_dir() and (directory / guest_pdb).is_file()
        )

    if not directories:
        print(f"No data-point folders containing {guest_pdb} found in {data_dir}.", file=sys.stderr)
        return 1

    failed = 0
    parametrized = 0
    for directory in directories:
        try:
            run(directory)
            parametrized += 1
            print(f"{directory.name}: parametrized {guest_pdb}", flush=True)
        except (OSError, ValueError, subprocess.CalledProcessError) as error:
            failed += 1
            if isinstance(error, subprocess.CalledProcessError):
                details = (error.stderr or error.stdout or "").strip() or str(error)
            else:
                details = str(error)
            print(f"{directory.name}: FAILED: {details}", flush=True)

    print(
        f"parametrization finished: {parametrized} parametrized, {failed} failed",
        flush=True,
    )
    return 1 if failed else 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Parametrize guest PDBs in data-point folders. The guest charge comes from "
                    "config.toml parametrize.net_charge, or with \"auto\" from each guest_input.json."
    )
    parser.add_argument(
        "names",
        nargs="*",
        metavar="DATA_POINT",
        help="data-point folders to process (default: folders containing the guest PDB)",
    )
    args = parser.parse_args(argv)
    return parametrize_data_points(names=args.names or None)


if __name__ == "__main__":
    raise SystemExit(main())
