import subprocess

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
