import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))
import parametrize as parametrize_module


class ParametrizeCliTests(unittest.TestCase):
    def test_scans_only_folders_with_guest_pdb(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            data_dir = Path(temporary_directory)
            with_guest = data_dir / "with_guest"
            without_guest = data_dir / "without_guest"
            with_guest.mkdir()
            without_guest.mkdir()
            (with_guest / "UNL.pdb").write_text("PDB\n")

            with patch.object(parametrize_module, "parametrize") as run_parametrize:
                result = parametrize_module.parametrize_data_points(data_dir=data_dir)

        self.assertEqual(result, 0)
        run_parametrize.assert_called_once_with(with_guest)

    def test_auto_charge_comes_from_guest_input_json(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            (Path(temporary_directory) / "guest_input.json").write_text(json.dumps({"net_charge": -1}))
            with patch.dict(parametrize_module.CFG["parametrize"], {"net_charge": "auto"}):
                with patch.object(parametrize_module.subprocess, "run") as run:
                    parametrize_module.parametrize(temporary_directory)

        self.assertEqual(run.call_count, 2)
        antechamber_arguments = run.call_args_list[0].args[0]
        self.assertEqual(antechamber_arguments[antechamber_arguments.index("-nc") + 1], "-1")

    def test_auto_charge_without_guest_input_fails_that_point_only(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            data_dir = Path(temporary_directory)
            for name in ("a", "b"):
                (data_dir / name).mkdir()
                (data_dir / name / "UNL.pdb").write_text("PDB\n")
            (data_dir / "b" / "guest_input.json").write_text(json.dumps({"net_charge": 0}))

            with patch.dict(parametrize_module.CFG["parametrize"], {"net_charge": "auto"}):
                with patch.object(parametrize_module.subprocess, "run") as run:
                    result = parametrize_module.parametrize_data_points(data_dir=data_dir)

        self.assertEqual(result, 1)               # a failed ...
        self.assertEqual(run.call_count, 2)       # ... b was still parametrized


if __name__ == "__main__":
    unittest.main()
