import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))
import parametrize as parametrize_module


class ParametrizeCliTests(unittest.TestCase):
    def test_scans_only_folders_with_guest_pdb_and_uses_charge(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            data_dir = Path(temporary_directory)
            with_guest = data_dir / "with_guest"
            without_guest = data_dir / "without_guest"
            with_guest.mkdir()
            without_guest.mkdir()
            (with_guest / "UNL.pdb").write_text("PDB\n")

            with patch.object(parametrize_module, "parametrize") as run_parametrize:
                result = parametrize_module.parametrize_data_points(
                    data_dir=data_dir, net_charge=-1
                )

        self.assertEqual(result, 0)
        run_parametrize.assert_called_once_with(with_guest, net_charge=-1)

    def test_standalone_auto_charge_requires_net_charge_argument(self):
        with patch.dict(parametrize_module.CFG["parametrize"], {"net_charge": "auto"}):
            with self.assertRaises(SystemExit) as error:
                parametrize_module.main([])

        self.assertEqual(error.exception.code, 2)

    def test_existing_main_call_style_still_supplies_guest_charge(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            with patch.dict(parametrize_module.CFG["parametrize"], {"net_charge": "auto"}):
                with patch.object(parametrize_module.subprocess, "run") as run:
                    parametrize_module.parametrize(
                        temporary_directory, net_charge=-1
                    )

        self.assertEqual(run.call_count, 2)
        antechamber_arguments = run.call_args_list[0].args[0]
        self.assertEqual(antechamber_arguments[antechamber_arguments.index("-nc") + 1], "-1")


if __name__ == "__main__":
    unittest.main()
