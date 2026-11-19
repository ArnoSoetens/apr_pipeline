import sys
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))
import smiles_to_2d


class SmilesTo2DTests(unittest.TestCase):
    def test_pipeline_restraint_atoms_are_on_the_vertical_x_zero_axis(self):
        molecule = Chem.MolFromSmiles("CCCCCC")

        molecule, anchors = smiles_to_2d.prepare_molecule(molecule)

        conformer = molecule.GetConformer()
        first = conformer.GetAtomPosition(anchors[0])
        second = conformer.GetAtomPosition(anchors[1])
        self.assertAlmostEqual(first.x, 0.0, places=7)
        self.assertAlmostEqual(second.x, 0.0, places=7)
        self.assertAlmostEqual((first.y + second.y) / 2, 0.0, places=7)
        self.assertFalse(conformer.Is3D())

    def test_rendered_image_contains_both_guide_line_styles(self):
        molecule, anchors = smiles_to_2d.prepare_molecule(Chem.MolFromSmiles("CCCCCC"))

        image = smiles_to_2d.draw_molecule(molecule, "test", 240, anchors)

        colors = {color: count for count, color in image.getcolors(maxcolors=240 * 268)}
        self.assertGreater(colors.get((190, 35, 35), 0), 10)
        self.assertGreater(colors.get((225, 105, 105), 0), 10)

    def test_main_writes_guided_grid_and_single_images(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            input_path = root / "molecules.txt"
            input_path.write_text("CCCCCC alpha\nbeta CCCCO\n")
            output_path = root / "images"

            smiles_to_2d.main(
                [
                    str(input_path),
                    str(output_path),
                    "--size",
                    "240",
                    "--cols",
                    "2",
                    "--single",
                ]
            )

            self.assertTrue((output_path / "grid_001.png").is_file())
            single_images = list((output_path / "single").glob("*.png"))
            self.assertEqual(len(single_images), 2)


if __name__ == "__main__":
    unittest.main()
