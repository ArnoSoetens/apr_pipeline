import subprocess

from APR.apr_refactor.apr_pipeline.apr_config import CFG

def parametrize(data_dir, guest_pdb=f'{CFG["naming"]["guest_prefix"]}.pdb', guest_prefix=CFG["naming"]["guest_prefix"]):
    result = subprocess.run(
    [
        "antechamber",
        "-i", f"{guest_pdb}",
        "-fi", "pdb",
        "-o", f"{guest_prefix}.mol2",
        "-fo", "mol2",
        "-c", CFG["parametrize"]["charge_method"],
        "-pf", "y",
        "-nc", str(CFG["parametrize"]["net_charge"]),
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
