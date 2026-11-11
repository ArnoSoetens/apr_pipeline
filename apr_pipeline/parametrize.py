import argparse
from pathlib import Path
import subprocess
import sys

from apr_config import CFG


def parametrize(
    data_dir,
    guest_pdb=f'{CFG["naming"]["guest_prefix"]}.pdb',
    guest_prefix=CFG["naming"]["guest_prefix"],
    net_charge=None,
):
    configured_charge = CFG["parametrize"]["net_charge"]
    if configured_charge == "auto":
        if net_charge is None:
            raise ValueError("net_charge is required when parametrize.net_charge = 'auto'")
        charge = net_charge
    else:
        charge = configured_charge

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


def parametrize_data_points(data_dir=CFG["directories"]["data"], names=None, net_charge=None):
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
        pdb_path = directory / guest_pdb
        if not pdb_path.is_file():
            failed += 1
            print(f"{directory.name}: FAILED: guest PDB not found: {pdb_path}", flush=True)
            continue

        try:
            parametrize(directory, net_charge=net_charge)
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
        description="Parametrize guest PDBs in data-point folders."
    )
    parser.add_argument(
        "names",
        nargs="*",
        metavar="DATA_POINT",
        help="data-point folders to process (default: folders containing the guest PDB)",
    )
    parser.add_argument(
        "--net-charge",
        type=int,
        help="overall guest charge (required when parametrize.net_charge = 'auto')",
    )
    args = parser.parse_args(argv)

    configured_charge = CFG["parametrize"]["net_charge"]
    if configured_charge == "auto" and args.net_charge is None:
        parser.error(
            "--net-charge is required because parametrize.net_charge is configured as 'auto'"
        )
    if configured_charge != "auto" and args.net_charge is not None:
        parser.error(
            "--net-charge cannot be used when parametrize.net_charge is a fixed integer"
        )

    return parametrize_data_points(names=args.names or None, net_charge=args.net_charge)


if __name__ == "__main__":
    raise SystemExit(main())
