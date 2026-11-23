import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

tests_dir = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(tests_dir / "apr_pipeline"))
import datapoints
from datapoints import WINDOWS, ensure_data_directories, has_complete_run


class PreparedDataPointTests(unittest.TestCase):
    def test_completed_run_requires_a_finished_trajectory_in_every_window(self):
        trajectory = datapoints.CFG["analysis"]["trajectory"]
        with tempfile.TemporaryDirectory() as temporary_directory:
            data_point = Path(temporary_directory) / "prepared"
            windows_dir = data_point / "windows"
            for window in WINDOWS:
                window_dir = windows_dir / window
                window_dir.mkdir(parents=True)
                (window_dir / trajectory).write_text("trajectory")
                (window_dir / datapoints.DONE_MARKER).write_text("")

            self.assertTrue(has_complete_run(data_point))
            (windows_dir / WINDOWS[-1] / datapoints.DONE_MARKER).unlink()   # interrupted: partial trajectory
            self.assertFalse(has_complete_run(data_point))
            (windows_dir / WINDOWS[-1] / trajectory).unlink()
            self.assertFalse(has_complete_run(data_point))


class EnsureDataDirectoriesTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.data_dir = self.root / "data"
        self.data_dir.mkdir()
        self.smiles = self.root / "smiles.txt"

    def tearDown(self):
        self._tmp.cleanup()

    def test_existing_data_points_need_no_smiles_list(self):
        # placed but unparametrized guests, and no SMILES list anywhere (the original bug)
        dp = self.data_dir / "PFBA_m_down"
        dp.mkdir()
        for name in ("pfba_m_down.pdb", "UNL.pdb", "BCD-UNL.pdb"):
            (dp / name).write_text("x\n")

        with patch.object(datapoints.subprocess, "run") as run:
            dps = ensure_data_directories(self.data_dir, self.smiles)

        run.assert_not_called()
        self.assertEqual(dps, [dp])

    def test_existing_data_points_win_over_smiles_list(self):
        dp = self.data_dir / "iso_00001"
        dp.mkdir()
        (dp / "iso_00001.smi").write_text("CC\n")
        self.smiles.write_text("CCC\n")

        with patch.object(datapoints.subprocess, "run") as run:
            ensure_data_directories(self.data_dir, self.smiles)

        run.assert_not_called()

    def test_empty_data_dir_is_filled_from_smiles_list(self):
        self.smiles.write_text("CC\n")

        with patch.object(datapoints.subprocess, "run") as run:
            ensure_data_directories(self.data_dir, self.smiles)

        run.assert_called_once_with(
            [sys.executable, str(tests_dir / "apr_pipeline" / "create_data_directories.py"),
             str(self.smiles), str(self.data_dir)],
            check=True,
        )

    def test_empty_data_dir_without_smiles_list_is_a_clear_error(self):
        with self.assertRaisesRegex(FileNotFoundError, "no data point folders"):
            ensure_data_directories(self.data_dir, self.smiles)

    def test_creates_real_folders_from_smiles_list(self):
        self.smiles.write_text("CC\nCCC\n")

        dps = ensure_data_directories(self.data_dir, self.smiles)

        self.assertEqual([dp.name for dp in dps], ["iso_00001", "iso_00002"])


if __name__ == "__main__":
    unittest.main()
