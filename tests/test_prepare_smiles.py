import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))

import prepare_smiles
import create_data_directories


class PrepareSmilesTests(unittest.TestCase):
    def test_functionalization_can_be_disabled(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            smiles_path = root / "smiles.txt"
            smiles_path.write_text("CCC\n\n")
            directories = {
                "smiles_list": str(smiles_path),
            }
            config = {
                "directories": directories,
                "functionalize": {
                    "enabled": False,
                    "groups": {"invalid group is ignored": "not-smiles"},
                    "site_strategy": "backbone_ends",
                },
            }

            with patch.object(prepare_smiles, "CFG", config):
                prepare_smiles.main()

            self.assertEqual(smiles_path.read_text(), "CCC\n")

    def test_functionalization_writes_products_to_the_same_smiles_list(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            smiles_path = root / "smiles.txt"
            smiles_path.write_text("CC\n")
            config = {
                "directories": {
                    "smiles_list": str(smiles_path),
                },
                "functionalize": {
                    "enabled": True,
                    "groups": {"hydroxyl": "O"},
                    "site_strategy": "backbone_ends",
                },
            }

            with patch.object(prepare_smiles, "CFG", config):
                prepare_smiles.main()

            self.assertEqual(smiles_path.read_text(), "CCO\n")

    def test_read_and_write_plain_smiles_skips_blank_lines(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            smiles_input = root / "input.txt"
            smiles_input.write_text("CC\n \n")
            smiles_path = root / "smiles.txt"

            records, skipped = prepare_smiles.read_smiles(smiles_input)
            prepare_smiles.write_smiles(
                [record["smiles"] for record in records], smiles_path
            )

            self.assertEqual(skipped, 1)
            self.assertEqual(smiles_path.read_text(), "CC\n")

    def test_hydrogen_is_replaced_and_products_are_deduplicated(self):
        products = prepare_smiles.make_products(
            [{"isomer_id": "mol_1", "smiles": "CC"}],
            {"hydroxyl": "O"},
            "backbone_ends",
        )

        self.assertEqual(products, ["CCO"])

    def test_fallback_replaces_f_only_at_terminal_cf3(self):
        products = prepare_smiles.make_products(
            [{"isomer_id": "mol_1", "smiles": "CC(F)(F)F"}],
            {"hydroxyl": "O"},
            "backbone_ends",
        )

        self.assertEqual(set(products), {"OCC(F)(F)F", "CC(O)(F)F"})

    def test_group_formal_charges_are_preserved(self):
        products = prepare_smiles.make_products(
            [{"isomer_id": "mol_1", "smiles": "CC"}],
            {
                "carboxylate": "C(=O)[O-]",
                "sulfonate": "S(=O)(=O)[O-]",
                "ammonium": "[NH3+]",
            },
            "backbone_ends",
        )

        self.assertTrue(products)
        for smiles in products:
            molecule = Chem.MolFromSmiles(smiles)
            self.assertIsNotNone(molecule)
            self.assertIn(Chem.GetFormalCharge(molecule), (-1, 1))

    def test_directory_loader_writes_one_folder_to_smiles_csv(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            smiles_path = root / "smiles.txt"
            smiles_path.write_text("CCO\nCCN\n")
            output_dir = root / "isomers"

            with patch.object(
                create_data_directories.sys,
                "argv",
                [
                    "create_data_directories.py",
                    str(smiles_path),
                    str(output_dir),
                ],
            ):
                create_data_directories.main()

            with open(output_dir / "isomers.csv", newline="") as mapping_file:
                rows = list(csv.DictReader(mapping_file))

            self.assertEqual(
                rows,
                [
                    {"folder": "iso_00001", "smiles": "CCO"},
                    {"folder": "iso_00002", "smiles": "CCN"},
                ],
            )
            self.assertEqual(
                (output_dir / "iso_00001" / "iso_00001.smi").read_text(),
                "CCO\n",
            )

    def test_directory_loader_rejects_non_txt_input(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            smiles_path = root / "smiles.xml"
            smiles_path.write_text("<smiles>CCO</smiles>\n")
            output_dir = root / "isomers"

            with (
                patch.object(
                    create_data_directories.sys,
                    "argv",
                    [
                        "create_data_directories.py",
                        str(smiles_path),
                        str(output_dir),
                    ],
                ),
                self.assertRaises(SystemExit) as error,
            ):
                create_data_directories.main()

            self.assertEqual(error.exception.code, 1)


if __name__ == "__main__":
    unittest.main()
