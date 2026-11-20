from pathlib import Path
import subprocess

from apr_config import CFG
from datapoints import refuse_simulated
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
    refuse_simulated(dp)
    dp = Path(dp)
    ensure_guest_input(dp)
    guest_pdb = dp / f'{CFG["naming"]["guest_prefix"]}.pdb'
    if not guest_pdb.is_file():
        raise FileNotFoundError(f"guest PDB not found: {guest_pdb}")
    parametrize(dp)


def main(argv=None):
    """python parametrize.py [NAME ...]; the guest charge comes from config.toml
    parametrize.net_charge, or with "auto" from each guest_input.json."""
    from setup_runner import stage_main
    return stage_main("parametrize", argv)


if __name__ == "__main__":
    raise SystemExit(main())
