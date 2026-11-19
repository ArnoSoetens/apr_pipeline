import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))
from create_data_directories import isomer_rows


class IsomerRowsTests(unittest.TestCase):
    def test_skips_directories_without_matching_smiles_files(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            data_dir = Path(temporary_directory)
            with_smiles = data_dir / "with_smiles"
            without_smiles = data_dir / "without_smiles"
            with_smiles.mkdir()
            without_smiles.mkdir()
            (with_smiles / "with_smiles.smi").write_text("CCO\n")
            (data_dir / "not_a_directory").write_text("ignored")

            self.assertEqual(isomer_rows(data_dir), [["with_smiles", "CCO"]])


if __name__ == "__main__":
    unittest.main()
