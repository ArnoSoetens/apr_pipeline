import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Geometry import Point3D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "apr_pipeline"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import placement
from datapoints import CLASH_PDB, COMPLEX_PDB, GUEST_PDB, next_stage
from test_main_guard import temporary_project


def molecule_at(x):
    molecule = Chem.MolFromSmiles("C")
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    conformer.SetAtomPosition(0, Point3D(x, 0.0, 0.0))
    molecule.AddConformer(conformer)
    return molecule


class PlacementClashDiagnosticsTests(unittest.TestCase):
    def test_clash_saves_guest_pose_to_clash_pdb(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory)
            (output_dir / "UNL.pdb").write_text("input\n")

            with self.assertRaisesRegex(ValueError, "guest atom 0"):
                placement.generate_pdb(
                    molecule_at(0.0),
                    molecule_at(0.25),
                    out_dir=output_dir,
                    clash_cutoff=1.0,
                )

            guest = Chem.MolFromPDBFile(
                str(output_dir / "UNL_clash.pdb"), removeHs=False
            )
            self.assertIsNotNone(guest)
            self.assertEqual(guest.GetNumAtoms(), 1)
            self.assertEqual((output_dir / "UNL.pdb").read_text(), "input\n")
            self.assertFalse((output_dir / "BCD-UNL.pdb").exists())
            self.assertAlmostEqual(
                guest.GetConformer().GetAtomPosition(0).x, 0.25, places=2
            )
            self.assertFalse((output_dir / "clash.json").exists())

    def test_successful_retry_removes_previous_failed_guest_pose(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory)
            with self.assertRaisesRegex(ValueError, "guest atom 0"):
                placement.generate_pdb(
                    molecule_at(0.0),
                    molecule_at(0.25),
                    out_dir=output_dir,
                    clash_cutoff=1.0,
                )

            placement.generate_pdb(
                molecule_at(0.0),
                molecule_at(10.0),
                out_dir=output_dir,
                clash_cutoff=1.0,
            )

            guest = Chem.MolFromPDBFile(
                str(output_dir / "UNL.pdb"), removeHs=False
            )
            self.assertIsNotNone(guest)
            self.assertAlmostEqual(
                guest.GetConformer().GetAtomPosition(0).x, 10.0, places=2
            )
            self.assertFalse((output_dir / "UNL_clash.pdb").exists())



class PlacementScriptClashTests(unittest.TestCase):
    """placement.py for real on a temporary project (never data/), with a cutoff that always clashes."""

    def test_clash_keeps_the_input_and_a_later_success_cleans_up(self):
        from structure import smiles_to_pdb

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = temporary_project(root)
            config = root / "config.toml"
            normal = config.read_text()
            config.write_text(re.sub(r"(?m)^clash_cutoff\s*=\s*[0-9.]+", "clash_cutoff = 100", normal))
            env = {**os.environ, "APR_CONFIG": str(config)}

            dp = data / "iso_00001"
            dp.mkdir()
            smiles_to_pdb("[O-]C(=O)C(F)(F)C(F)(F)C(F)(F)F", str(dp / GUEST_PDB))  # as structure leaves it
            unplaced = (dp / GUEST_PDB).read_text()
            (dp / "guest_input.json").write_text('{"net_charge": -1}')

            def place():
                return subprocess.run([sys.executable, str(ROOT / "apr_pipeline" / "placement.py"),
                                       "iso_00001"], env=env, capture_output=True, text=True, timeout=300)

            result = place()
            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            self.assertIn("FAILED in placement: clash", result.stdout)
            self.assertTrue((dp / CLASH_PDB).is_file())
            self.assertEqual((dp / GUEST_PDB).read_text(), unplaced)
            self.assertFalse((dp / COMPLEX_PDB).exists())
            self.assertIn("guest atom", (dp / "setup_failed.txt").read_text())
            self.assertEqual(next_stage(dp), "placement")

            config.write_text(normal)
            result = place()
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertFalse((dp / CLASH_PDB).exists())
            self.assertFalse((dp / "setup_failed.txt").exists())
            self.assertEqual(next_stage(dp), "parametrize")


if __name__ == "__main__":
    unittest.main()
