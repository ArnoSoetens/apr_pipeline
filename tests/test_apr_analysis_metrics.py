import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))

from apr_analysis import (
    _amber_distance_wall_energy,
    _attachment_samples,
    _distance_wall_summary,
    _fraction_outside_distance,
    _phase_details,
    _read_distance_walls,
)


class AnalysisMetricTests(unittest.TestCase):
    def test_attachment_samples_exclude_shared_pull_start_window(self):
        attach_windows = [[np.array([6.0])] for _ in range(16)]
        attach_windows.append([np.array([100.0])])  # p000, shared with pulling

        samples = _attachment_samples({"attach": attach_windows}, attach_window_count=16)

        self.assertEqual(samples.size, 16)
        self.assertTrue(np.all(samples == 6.0))

    def test_attachment_samples_reject_an_invalid_window_count(self):
        with self.assertRaisesRegex(ValueError, "Expected 3 attachment windows, found 2"):
            _attachment_samples({"attach": [[np.array([1.0])], [np.array([2.0])]]}, 3)

    def test_attachment_distance_fraction_uses_absolute_deviation_from_initial(self):
        samples = np.array([6.0, -1.1, 13.0, 13.1, 6.5, 13.0])
        self.assertAlmostEqual(
            _fraction_outside_distance(samples, initial_distance=6.0, threshold=7.0),
            2 / 6,
        )

    def test_phase_results_have_one_general_estimate_per_method(self):
        results = {
            "attach": {
                "ti-block": {"fe": 2.0, "sem": 0.2},
                "mbar-autoc": {"fe": 1.9, "sem": 0.1},
            }
        }

        details = _phase_details(results, "attach")

        self.assertEqual(set(details), {"ti-block", "mbar-autoc"})
        self.assertEqual(details["ti-block"], {"fe": 2.0, "sem": 0.2, "units": "kcal/mol"})
        self.assertEqual(details["mbar-autoc"], {"fe": 1.9, "sem": 0.1, "units": "kcal/mol"})

    def test_distance_wall_parser_selects_distance_and_ignores_angle(self):
        contents = (
            "&rst iat= 149,148, r1= 0, r2= 6, r3= 6, r4= 999, "
            "rk2= 10, rk3= 10, &end\n"
            "&rst iat= 150,149,148, r1= 90, r2= 120, r3= 180, r4= 180, "
            "rk2= 100, rk3= 100, &end\n"
            "&rst iat= 149,148, r1= 0, r2= 1, r3= 11, r4= 12, "
            "rk2= 100, rk3= 100, &end\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            restraint_file = Path(directory) / "disang.rest"
            restraint_file.write_text(contents)

            wall = _read_distance_walls(restraint_file)[0]

        self.assertEqual(wall["atom_indices"], [149, 148])
        self.assertEqual((wall["r1"], wall["r2"], wall["r3"], wall["r4"]), (0, 1, 11, 12))

    def test_amber_wall_energy_has_flat_harmonic_and_linear_regions(self):
        wall = {"r1": 0.0, "r2": 1.0, "r3": 11.0, "r4": 12.0, "rk2": 100.0, "rk3": 100.0}
        energies = _amber_distance_wall_energy(
            np.array([-0.5, 0.5, 1.0, 6.0, 11.0, 11.5, 12.0, 12.5]),
            wall,
        )

        np.testing.assert_allclose(energies, [200.0, 25.0, 0.0, 0.0, 0.0, 25.0, 100.0, 200.0])

    def test_distance_wall_summary_reports_potential_and_per_window_net_work(self):
        wall_line = (
            "&rst iat= 149,148, r1= 0, r2= 1, r3= 11, r4= 12, "
            "rk2= 100, rk3= 100, &end\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            windows_dir = Path(directory)
            for name in ("a000", "a001", "p000"):
                window = windows_dir / name
                window.mkdir()
                (window / "disang.rest").write_text(wall_line if name.startswith("a") else "")

            result = _distance_wall_summary(
                windows_dir,
                {
                    "attach": [
                        [np.array([0.5, 1.0, 11.5, 12.5])],
                        [np.array([1.5, 10.0, 12.0])],
                        [np.array([100.0])],
                    ]
                },
                ["a000", "a001", "p000"],
            )

        self.assertEqual(result["frame_count"], 7)
        self.assertEqual(set(result["per_window"]), {"a000", "a001"})
        self.assertAlmostEqual(result["per_window"]["a000"]["net_work_by_wall_start_to_end_kcal_per_mol"], -175.0)
        self.assertAlmostEqual(result["per_window"]["a001"]["net_work_by_wall_start_to_end_kcal_per_mol"], -100.0)
        self.assertNotIn("p000", result["per_window"])

    def test_distance_wall_summary_supports_both_rim_walls(self):
        wall_lines = (
            "&rst iat= 149,148, r1= 5, r2= 6, r3= 999, r4= 1000, "
            "rk2= 100, rk3= 100, &end\n"
            "&rst iat= 150,148, r1= 0, r2= 0, r3= 9, r4= 10, "
            "rk2= 100, rk3= 100, &end\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            windows_dir = Path(directory)
            for name in ("a000", "a001"):
                window = windows_dir / name
                window.mkdir()
                (window / "disang.rest").write_text(wall_lines)

            with (
                patch(
                    "apr_analysis.load_trajectory",
                    side_effect=lambda *args: Path(args[0]).name,
                ),
                patch(
                    "apr_analysis.read_restraint_data",
                    side_effect=lambda *args: (
                        np.array([8.5, 9.5]) if args[1].mask1 == "@150" else None
                    ),
                ),
            ):
                result = _distance_wall_summary(
                    windows_dir,
                    {
                        "attach": [
                            [np.array([5.5, 6.0])],
                            [np.array([5.0, 6.5])],
                        ]
                    },
                    ["a000", "a001"],
                    topology="complex-dum.prmtop",
                )

        self.assertEqual(len(result["walls"]), 2)
        self.assertEqual(
            [wall["parameters"]["amber_atom_indices_1based"] for wall in result["walls"]],
            [[149, 148], [150, 148]],
        )
        self.assertEqual(result["walls"][1]["frame_count"], 4)
        self.assertAlmostEqual(
            result["walls"][1]["per_window"]["a000"][
                "net_work_by_wall_start_to_end_kcal_per_mol"
            ],
            -25.0,
        )

    def test_distance_wall_summary_returns_none_when_wall_is_absent(self):
        with tempfile.TemporaryDirectory() as directory:
            window = Path(directory) / "a000"
            window.mkdir()
            (window / "disang.rest").write_text(
                "&rst iat= 149,148, r1= 0, r2= 6, r3= 6, r4= 999, "
                "rk2= 10, rk3= 10, &end\n"
            )
            result = _distance_wall_summary(
                directory,
                {"attach": [[np.array([6.0])]]},
                ["a000"],
            )

        self.assertIsNone(result)


if __name__ == "__main__":
    unittest.main()
