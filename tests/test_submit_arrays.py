import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import submit_arrays


class SubmitArraysTests(unittest.TestCase):
    def test_default_package_reserves_25_percent_of_task_time(self):
        with patch.dict(
            submit_arrays.SL,
            {"sim_cpus": 4, "sim_task_hours": 4, "sim_minutes_per_window": 5},
        ):
            self.assertEqual(submit_arrays.default_package(), 4 * 36)

    def test_default_package_rejects_task_too_short_for_a_window(self):
        with patch.dict(
            submit_arrays.SL,
            {"sim_cpus": 1, "sim_task_hours": 1, "sim_minutes_per_window": 50},
        ):
            with self.assertRaisesRegex(ValueError, "reserving 25%"):
                submit_arrays.default_package()

    def test_simulation_time_limit_uses_25_percent_margin(self):
        with patch.dict(
            submit_arrays.SL,
            {"sim_cpus": 4, "sim_task_hours": 4, "sim_minutes_per_window": 5},
        ):
            package = submit_arrays.default_package()
            self.assertEqual(submit_arrays.sim_time(package), "03:55:00")

    def test_analysis_defaults_isolate_data_points_with_more_memory(self):
        self.assertEqual(submit_arrays.SL["analysis_data_points_per_task"], 1)
        self.assertEqual(submit_arrays.points_resources("analysis"), (1, "--mem-per-cpu=4G"))


if __name__ == "__main__":
    unittest.main()
