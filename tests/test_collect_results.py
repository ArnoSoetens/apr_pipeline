import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))

from collect_results import collect_result_rows
from datapoints import WINDOWS
from apr_config import CFG


class CollectResultsTests(unittest.TestCase):
    def test_collects_completed_prepared_results_without_smiles(self):
        trajectory = CFG["analysis"]["trajectory"]
        with tempfile.TemporaryDirectory() as temporary_directory:
            data_dir = Path(temporary_directory)
            data_point = data_dir / "prepared"
            data_point.mkdir()
            (data_point / "results.json").write_text(
                '{"dG_bind": -1.5, "dG_sem": 0.2}'
            )
            windows_dir = data_point / "windows"
            for window in WINDOWS:
                window_dir = windows_dir / window
                window_dir.mkdir(parents=True)
                (window_dir / trajectory).write_text("trajectory")

            rows, missing = collect_result_rows(data_dir)

        self.assertEqual(rows, [["prepared", "", -1.5, 0.2]])
        self.assertEqual(missing, 0)

    def test_ignores_incomplete_run_without_results(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            data_dir = Path(temporary_directory)
            data_point = data_dir / "incomplete"
            (data_point / "windows").mkdir(parents=True)

            rows, missing = collect_result_rows(data_dir)

        self.assertEqual(rows, [])
        self.assertEqual(missing, 0)


if __name__ == "__main__":
    unittest.main()
