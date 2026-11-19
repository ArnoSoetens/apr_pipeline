import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import simulate_guest


class GuestInputTests(unittest.TestCase):
    def write_guest_files(self, directory, pdb_names=("C1", "O1"), mol2_names=("C1", "O1")):
        pdb = directory / "guest.pdb"
        mol2 = directory / "guest.mol2"
        pdb.write_text(
            "".join(
                f"HETATM{index:5d} {name:>4s} UNL A{1:4d}    "
                "0.000   0.000   0.000  1.00  0.00           C\n"
                for index, name in enumerate(pdb_names, start=1)
            ),
            encoding="utf-8",
        )
        mol2.write_text(
            "@<TRIPOS>MOLECULE\nUNL\n"
            f"{len(mol2_names)} 0 1 0 0\nSMALL\nUSER_CHARGES\n"
            "@<TRIPOS>ATOM\n"
            + "".join(
                f"{index} {name} 0.0 0.0 0.0 C.3 1 UNL 0.0\n"
                for index, name in enumerate(mol2_names, start=1)
            )
            + "@<TRIPOS>BOND\n",
            encoding="utf-8",
        )
        return pdb, mol2

    def test_guest_atom_names_may_be_in_different_order(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            pdb, mol2 = self.write_guest_files(
                Path(temp_dir), pdb_names=("C1", "O1"), mol2_names=("O1", "C1")
            )
            simulate_guest.validate_guest_inputs(pdb, mol2)

    def test_mismatched_atom_names_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            pdb, mol2 = self.write_guest_files(
                Path(temp_dir), pdb_names=("C1", "O1"), mol2_names=("C1", "N1")
            )
            with self.assertRaisesRegex(ValueError, "same atom names"):
                simulate_guest.validate_guest_inputs(pdb, mol2)

    def test_rendered_inputs_match_config_without_apr_restraints(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            simulate_guest.render_inputs(output_dir)
            minimize = (output_dir / "minimize.in").read_text(encoding="utf-8")
            equilibration = (output_dir / "equilibrate.in").read_text(encoding="utf-8")
            production = (output_dir / "production.in").read_text(encoding="utf-8")
            run = simulate_guest.CFG["run"]

        self.assertIn(f"igb =           {run['igb']}", minimize)
        self.assertIn(f"dt =            {run['dt']}", equilibration)
        self.assertIn(f"temp0 =         {run['temperature']}", production)
        self.assertIn("ntr =           0", production)
        self.assertNotIn("DISANG", production)
        self.assertIn("ntx =           5", production)
        self.assertIn("irest =         1", production)


if __name__ == "__main__":
    unittest.main()
