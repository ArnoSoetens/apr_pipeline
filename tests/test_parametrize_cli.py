import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))
import parametrize as parametrize_module
import setup_runner
from datapoints import COMPLEX_PDB, GUEST_PDB


def placed(dp):
    dp.mkdir()
    for name in (GUEST_PDB, COMPLEX_PDB):
        (dp / name).write_text("PDB\n")


class ParametrizeCliTests(unittest.TestCase):
    def test_without_names_only_folders_whose_next_stage_is_parametrize(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            data_dir = Path(temporary_directory)
            with_guest = data_dir / "with_guest"
            placed(with_guest)
            unplaced = data_dir / "unplaced"
            unplaced.mkdir()
            (unplaced / "my_guest.pdb").write_text("PDB\n")

            with patch.object(parametrize_module, "parametrize") as run_parametrize:
                result = setup_runner.stage_main("parametrize", [], data_dir=data_dir)

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
                placed(data_dir / name)
            (data_dir / "b" / "guest_input.json").write_text(json.dumps({"net_charge": 0}))

            with patch.dict(parametrize_module.CFG["parametrize"], {"net_charge": "auto"}):
                with patch.object(parametrize_module.subprocess, "run") as run:
                    result = setup_runner.stage_main("parametrize", ["a", "b"], data_dir=data_dir)
            failed = (data_dir / "a" / "setup_failed.txt").read_text()

        self.assertEqual(result, 1)
        self.assertIn("parametrize: set net_charge", failed)               # a failed ...
        self.assertEqual(run.call_count, 2)       # ... b was still parametrized


if __name__ == "__main__":
    unittest.main()
