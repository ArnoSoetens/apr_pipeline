import subprocess

def parametrize(data_dir, guest_pdb="UNL.pdb", guest_prefix="UNL"):
    result = subprocess.run(
    [
        "antechamber",
        "-i", f"{guest_pdb}",
        "-fi", "pdb",
        "-o", f"{guest_prefix}.mol2",
        "-fo", "mol2",
        "-c", "bcc",
        "-pf", "y",
        "-nc", "-1",
        "-at", "gaff2"
    ],
    cwd=data_dir,
    check=True,
    capture_output=True,
    text=True,
    )

    result2 = subprocess.run(
        [
            "parmchk2",
            "-i", "UNL.mol2",
            "-f", "mol2",
            "-o", f"{guest_prefix}.frcmod",
            "-s", "2"

        ],
        cwd=data_dir,
        check=True,
        capture_output=True,
        text=True,
    )
    return
