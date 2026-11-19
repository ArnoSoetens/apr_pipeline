import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))
import smiles_to_images


class SmilesToImagesTests(unittest.TestCase):
    def test_default_carbon_limit_filters_molecules_above_c9(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            input_path = root / "molecules.txt"
            input_path.write_text("CCCCCCCCC c9\nCCCCCCCCCC c10\n")
            output_path = root / "images"
            output = io.StringIO()

            with contextlib.redirect_stdout(output):
                smiles_to_images.main(
                    [str(input_path), str(output_path), "--size", "80"]
                )

            self.assertTrue((output_path / "grid_001.png").is_file())
            self.assertIn("1 molecules kept, 1 over 9 carbons filtered", output.getvalue())

    def test_carbon_limit_can_be_increased(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            input_path = root / "molecules.txt"
            input_path.write_text("CCCCCCCCC c9\nCCCCCCCCCC c10\n")
            output_path = root / "images"
            output = io.StringIO()

            with contextlib.redirect_stdout(output):
                smiles_to_images.main(
                    [
                        str(input_path),
                        str(output_path),
                        "--size",
                        "80",
                        "--max-carbons",
                        "10",
                    ]
                )

            self.assertIn("2 molecules kept, 0 over 10 carbons filtered", output.getvalue())


if __name__ == "__main__":
    unittest.main()
