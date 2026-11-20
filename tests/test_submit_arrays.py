import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import submit_arrays
from datapoints import (COMPLEX_PDB, FAILED_FILE, GUEST_FRCMOD, GUEST_MOL2, GUEST_PDB, TOPOLOGY,
                        TRAJECTORY, WINDOWS)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_main_guard import temporary_project


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



class StatusBySetupStageTests(unittest.TestCase):
    """submit_arrays.py status/setup for real on a temporary project (never data/; nothing submitted)."""

    def test_status_and_setup_follow_next_stage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = temporary_project(root)

            def point(name, *files):
                (data / name).mkdir()
                for filename in files:
                    (data / name / filename).write_text("x\n")
                return data / name

            point("s1", "s1.smi")
            point("p1", "x.pdb")
            point("q1", GUEST_PDB, COMPLEX_PDB)
            point("r1", GUEST_PDB, COMPLEX_PDB, GUEST_MOL2, GUEST_FRCMOD)
            failed = point("f1", "f1.smi", GUEST_PDB)
            (failed / FAILED_FILE).write_text("2026-11-20 10:00:00\nplacement: clash: guest atom 3\n")
            built = point("w1", GUEST_PDB, COMPLEX_PDB, GUEST_MOL2, GUEST_FRCMOD)
            for window in WINDOWS:
                (built / "windows" / window).mkdir(parents=True)
                (built / "windows" / window / TOPOLOGY).write_text("x\n")
            (built / "windows" / WINDOWS[0] / TRAJECTORY).write_text("x\n")
            (built / FAILED_FILE).write_text("2026-11-20 10:00:00\nold\n")   # stale: setup is done
            point("junk", "notes.txt")

            env = {**os.environ, "APR_CONFIG": str(root / "config.toml")}

            def submit_arrays_py(*args):
                result = subprocess.run([sys.executable, str(Path(submit_arrays.__file__)), *args],
                                        env=env, capture_output=True, text=True, timeout=120)
                self.assertEqual(result.returncode, 0, result.stderr)
                return result.stdout

            status = submit_arrays_py("status")
            for expected in ("6 data points", "1  need structure", "1  need placement",
                             "1  need parametrize", "1  need apr_setup", "1  failed last time",
                             "f1: placement: clash", "1  windows ready, simulation not finished",
                             f"windows still to simulate: {len(WINDOWS) - 1}",
                             "e.g. junk"):
                self.assertIn(expected, status)

            setup = submit_arrays_py("setup")              # dry run: shows, submits nothing
            self.assertIn("[setup] 4 data points", setup)   # s1, p1, q1, r1; not f1, w1
            self.assertIn("1 data point(s) failed earlier", setup)
            self.assertIn("[setup] 5 data points", submit_arrays_py("setup", "--retry-failed"))


if __name__ == "__main__":
    unittest.main()
