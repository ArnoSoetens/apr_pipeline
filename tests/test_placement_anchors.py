"""G1/G2 anchors in placement: automatic (smallest turn, lower end is G1) or chosen in guest_input.json.

placement.py runs for real (RDKit only) on a temporary project, never on data/.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "apr_pipeline"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import placement
from datapoints import GUEST_PDB
from test_main_guard import temporary_project

PFBA_ANION = "[O-]C(=O)C(F)(F)C(F)(F)C(F)(F)F"


def molecule_with_coordinates(coordinates):
    molecule = Chem.MolFromSmiles("CCC")
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index, (x, y, z) in enumerate(coordinates):
        conformer.SetAtomPosition(index, Point3D(x, y, z))
    molecule.AddConformer(conformer)
    return molecule


class AnchorGeometryTests(unittest.TestCase):
    def test_without_flip_g1_to_g2_always_ends_up_along_plus_axis(self):
        # b1 -> b2 points almost straight down: the smallest turn would keep it down
        molecule = molecule_with_coordinates([(0.0, 0.0, 1.0), (0.1, 0.0, -1.0), (0.0, 2.0, 0.0)])
        for allow_flip, sign in ((True, -1), (False, 1)):
            placed = placement.place_guest_in_cavity(
                molecule, 0, 1, cavity_center=[0, 0, 0], cavity_axis=[0, 0, 1],
                allow_flip=allow_flip)[0].GetConformer().GetPositions()
            axis = placed[1] - placed[0]
            with self.subTest(allow_flip=allow_flip):
                np.testing.assert_allclose(axis / np.linalg.norm(axis), [0, 0, sign], atol=1e-9)

    def test_order_anchors_puts_the_lower_atom_first(self):
        coords = [(0, 0, 3.0), (0, 0, -1.0)]
        self.assertEqual(placement.order_anchors(coords, 0, 1, [0, 0, 1]), (1, 0))
        self.assertEqual(placement.order_anchors(coords, 1, 0, [0, 0, 1]), (1, 0))
        self.assertEqual(placement.order_anchors(coords, 0, 1, [0, 0, -1]), (0, 1))

    def test_check_serials(self):
        with tempfile.TemporaryDirectory() as tmp:
            pdb = Path(tmp) / "g.pdb"
            molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
            from rdkit.Chem import AllChem
            AllChem.EmbedMolecule(molecule, randomSeed=1)
            Chem.MolToPDBFile(molecule, str(pdb))
            placement.check_serials(pdb)            # RDKit writes 1..N

            shifted = [line[:6] + f"{int(line[6:11]) + 100:5d}" + line[11:]
                       if line.startswith(("ATOM", "HETATM")) else line
                       for line in pdb.read_text().splitlines()]
            pdb.write_text("\n".join(shifted) + "\n")
            with self.assertRaisesRegex(ValueError, "does not number its atoms 1 to 9"):
                placement.check_serials(pdb)


class PlacementAnchorTests(unittest.TestCase):
    """placement.py on PFBA, given upright and upside down."""

    @classmethod
    def setUpClass(cls):
        from structure import backbone_anchors, smiles_to_pdb

        cls._tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls._tmp.name)
        cls.data = temporary_project(cls.root)
        cls.env = {**os.environ, "APR_CONFIG": str(cls.root / "config.toml")}

        cls.source = cls.root / "pfba.pdb"
        mol = smiles_to_pdb(PFBA_ANION, str(cls.source))
        mol = Chem.MolFromPDBFile(str(cls.source), removeHs=False)
        cls.chain_ends = backbone_anchors(mol)
        # the same guest turned 180 degrees about x: upside down
        flipped = Chem.Mol(mol)
        conformer = flipped.GetConformer()
        for i, (x, y, z) in enumerate(conformer.GetPositions()):
            conformer.SetAtomPosition(i, Point3D(x, -y, -z))
        cls.flipped_source = cls.root / "pfba_flipped.pdb"
        Chem.MolToPDBFile(flipped, str(cls.flipped_source), flavor=2)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def data_point(self, name, source, guest_input=None):
        dp = self.data / name
        dp.mkdir()
        (dp / "guest.pdb").write_text(Path(source).read_text())
        if guest_input is not None:
            (dp / "guest_input.json").write_text(json.dumps(guest_input))
        return dp

    def place(self, *names):
        result = subprocess.run([sys.executable, str(ROOT / "apr_pipeline" / "placement.py"), *names],
                                env=self.env, capture_output=True, text=True, timeout=300)
        return result.returncode, result.stdout + result.stderr

    @staticmethod
    def placed(dp):
        return Chem.MolFromPDBFile(str(dp / GUEST_PDB), removeHs=False).GetConformer().GetPositions()

    @staticmethod
    def saved_anchors(dp):
        return json.loads((dp / "guest_input.json").read_text())["anchors"]

    def assert_g1_below_g2(self, dp, serials):
        coords = self.placed(dp)
        self.assertLess(coords[serials[0] - 1][2], coords[serials[1] - 1][2])

    def test_automatic_anchors_keep_the_input_orientation_and_name_the_lower_end_g1(self):
        upright = self.data_point("auto_upright", self.source)
        flipped = self.data_point("auto_flipped", self.flipped_source)
        code, out = self.place("auto_upright", "auto_flipped")
        self.assertEqual(code, 0, out)

        a, b = (i + 1 for i in self.chain_ends)
        self.assertEqual(sorted(self.saved_anchors(upright)), sorted([a, b]))
        self.assertEqual(self.saved_anchors(flipped), self.saved_anchors(upright)[::-1])
        for dp in (upright, flipped):
            self.assert_g1_below_g2(dp, self.saved_anchors(dp))
        # the two poses stay 180 degrees apart: chain end a is low in one, high in the other
        self.assertLess(self.placed(upright)[a - 1][2] * self.placed(flipped)[a - 1][2], 0)

        # placing again with the saved anchors gives the same pose
        before = self.placed(flipped)
        code, out = self.place("auto_flipped")
        self.assertEqual(code, 0, out)
        np.testing.assert_allclose(self.placed(flipped), before, atol=1e-3)

    def test_chosen_anchors_are_kept_and_g1_ends_up_below_g2(self):
        a, b = (i + 1 for i in self.chain_ends)
        for name, source in (("chosen_upright", self.source), ("chosen_flipped", self.flipped_source)):
            for order in ([a, b], [b, a]):
                with self.subTest(source=name, anchors=order):
                    dp = self.data / f"{name}_{order[0]}"
                    dp.mkdir()
                    (dp / "guest.pdb").write_text(source.read_text())
                    (dp / "guest_input.json").write_text(json.dumps({"anchors": order}))
                    code, out = self.place(dp.name)
                    self.assertEqual(code, 0, out)
                    self.assertEqual(self.saved_anchors(dp), order)
                    self.assert_g1_below_g2(dp, order)

    def test_chosen_anchors_need_atoms_numbered_from_one(self):
        dp = self.data_point("renumbered", self.source, {"anchors": [1, 5]})
        pdb = dp / "guest.pdb"
        pdb.write_text("\n".join(
            line[:6] + f"{int(line[6:11]) + 100:5d}" + line[11:]
            if line.startswith(("ATOM", "HETATM")) else line
            for line in pdb.read_text().splitlines()) + "\n")
        code, out = self.place("renumbered")
        self.assertNotEqual(code, 0)
        self.assertIn("does not number its atoms 1 to", out)
        self.assertFalse((dp / GUEST_PDB).exists())

        # with automatic anchors the numbering does not matter
        (dp / "guest_input.json").write_text(json.dumps({"anchors": None}))
        code, out = self.place("renumbered")
        self.assertEqual(code, 0, out)

    def test_chosen_anchors_beyond_the_last_atom_are_refused(self):
        self.data_point("too_far", self.source, {"anchors": [1, 999]})
        code, out = self.place("too_far")
        self.assertNotEqual(code, 0)
        self.assertIn("has only", out)


if __name__ == "__main__":
    unittest.main()
