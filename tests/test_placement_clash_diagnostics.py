import sys
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Geometry import Point3D

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))
import placement


def molecule_at(x):
    molecule = Chem.MolFromSmiles("C")
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    conformer.SetAtomPosition(0, Point3D(x, 0.0, 0.0))
    molecule.AddConformer(conformer)
    return molecule


class PlacementClashDiagnosticsTests(unittest.TestCase):
    def test_clash_saves_guest_pose_to_normal_guest_pdb(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory)

            with self.assertRaisesRegex(ValueError, "guest atom 0"):
                placement.generate_pdb(
                    molecule_at(0.0),
                    molecule_at(0.25),
                    out_dir=output_dir,
                    clash_cutoff=1.0,
                )

            guest = Chem.MolFromPDBFile(
                str(output_dir / "UNL.pdb"), removeHs=False
            )
            self.assertIsNotNone(guest)
            self.assertEqual(guest.GetNumAtoms(), 1)
            self.assertAlmostEqual(
                guest.GetConformer().GetAtomPosition(0).x, 0.25, places=2
            )
            self.assertFalse((output_dir / "clash.json").exists())

    def test_successful_retry_overwrites_previous_failed_guest_pose(self):
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


if __name__ == "__main__":
    unittest.main()
