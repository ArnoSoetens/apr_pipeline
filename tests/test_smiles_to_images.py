import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import smiles_to_images


def run_main(args):
    output = io.StringIO()
    with contextlib.redirect_stdout(output):
        smiles_to_images.main([str(arg) for arg in args])
    return output.getvalue()


class SmilesToImagesTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.input = self.root / "molecules.txt"
        self.out = self.root / "images"

    def tearDown(self):
        self._tmp.cleanup()

    def test_default_carbon_limit_filters_molecules_above_c9(self):
        self.input.write_text("CCCCCCCCC c9\nCCCCCCCCCC c10\n")
        output = run_main([self.input, self.out, "--size", "80"])
        self.assertTrue((self.out / "grid_001.png").is_file())
        self.assertIn("1 molecules kept, 1 over 9 carbons filtered", output)

    def test_carbon_limit_can_be_increased(self):
        self.input.write_text("CCCCCCCCC c9\nCCCCCCCCCC c10\n")
        output = run_main([self.input, self.out, "--size", "80", "--max-carbons", "10"])
        self.assertIn("2 molecules kept, 0 over 10 carbons filtered", output)

    def test_without_restraints_no_guides_are_drawn(self):
        self.input.write_text("CCCCCC alpha\n")
        run_main([self.input, self.out, "--size", "120", "--single"])
        (single,) = (self.out / "single").glob("*.png")
        from PIL import Image
        colors = {color for _, color in Image.open(single).convert("RGB").getcolors(120 * 120)}
        self.assertNotIn((190, 35, 35), colors)
        self.assertNotIn((225, 105, 105), colors)

    def test_restraints_write_guided_grid_and_single_images(self):
        self.input.write_text("CCCCCC alpha\nbeta CCCCO\n")
        output = run_main([self.input, self.out, "--restraints", "--size", "240", "--cols", "2",
                           "--single"])
        self.assertTrue((self.out / "grid_001.png").is_file())
        self.assertEqual(len(list((self.out / "single").glob("*.png"))), 2)
        self.assertIn("2 molecules kept", output)
        self.assertIn("0 without restraint anchors", output)

    def test_molecule_without_anchors_is_reported_not_drawn(self):
        self.input.write_text("CCCC butane\nC methane\n")
        output = run_main([self.input, self.out, "--restraints", "--size", "120"])
        self.assertIn("1 molecules kept", output)
        self.assertIn("1 without restraint anchors", output)
        self.assertIn("methane", output)


class RestraintDrawingTests(unittest.TestCase):
    def test_pipeline_restraint_atoms_are_on_the_vertical_x_zero_axis(self):
        molecule, anchors = smiles_to_images.prepare_molecule(Chem.MolFromSmiles("CCCCCC"))
        conformer = molecule.GetConformer()
        first = conformer.GetAtomPosition(anchors[0])
        second = conformer.GetAtomPosition(anchors[1])
        self.assertAlmostEqual(first.x, 0.0, places=7)
        self.assertAlmostEqual(second.x, 0.0, places=7)
        self.assertAlmostEqual((first.y + second.y) / 2, 0.0, places=7)
        self.assertFalse(conformer.Is3D())

    def test_anchors_are_the_pipelines(self):
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))
        from structure import backbone_anchors

        smiles = "CCOCCCN"                       # ether: the chain runs through the O
        _, anchors = smiles_to_images.prepare_molecule(Chem.MolFromSmiles(smiles))
        self.assertEqual(anchors, backbone_anchors(Chem.MolFromSmiles(smiles)))
        self.assertEqual(set(anchors), {0, 5})

    def test_rendered_image_contains_both_guide_line_styles(self):
        molecule, anchors = smiles_to_images.prepare_molecule(Chem.MolFromSmiles("CCCCCC"))
        image = smiles_to_images.draw_molecule(molecule, "test", 240, anchors)
        colors = {color: count for count, color in image.getcolors(maxcolors=240 * 268)}
        self.assertGreater(colors.get((190, 35, 35), 0), 10)
        self.assertGreater(colors.get((225, 105, 105), 0), 10)


if __name__ == "__main__":
    unittest.main()
