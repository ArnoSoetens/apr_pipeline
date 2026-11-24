"""apr_analysis.run(): where the plots go, and which extras --enhanced adds. The pAPRika
calculation is mocked; the numbers themselves are checked against real data separately."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))
import apr_analysis
from datapoints import TRAJECTORY, WINDOWS


class StopAfterSetup(Exception):
    pass


class FakeFeCalc:
    """Records what calculate_free_energy configures, then stops."""
    seen = {}

    def collect_data(self):
        pass

    def compute_free_energy(self):
        FakeFeCalc.seen["fractions"] = self.fractions
        raise StopAfterSetup


class ApraAnalysisRunTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.dp = Path(self._tmp.name) / "iso_00001"
        for window in WINDOWS:
            (self.dp / "windows" / window).mkdir(parents=True)
            (self.dp / "windows" / window / TRAJECTORY).write_text("x\n")

    def tearDown(self):
        self._tmp.cleanup()

    def run_with_mocked_calculation(self, enhanced):
        calls = {}

        def fake_calculate(windows, **kwargs):
            calls.update(kwargs)
            return {"dG_bind": -3.0, "dG_sem": 0.2}, {"pull": []}

        with patch.object(apr_analysis, "calculate_free_energy", side_effect=fake_calculate), \
             patch.object(apr_analysis, "pull_window_histograms") as histograms:
            apr_analysis.run(self.dp, enhanced=enhanced)
        return calls, histograms

    def test_default_only_the_binding_plot_in_plots(self):
        calls, histograms = self.run_with_mocked_calculation(enhanced=False)
        plots = self.dp / "plots"
        self.assertEqual(calls["plot_path"], plots / "binding_free_energy.jpg")
        self.assertIsNone(calls["grid_path"])
        self.assertIsNone(calls["angle_grid_path"])
        self.assertIsNone(calls["convergence_path"])
        self.assertFalse(calls["enhanced"])
        histograms.assert_not_called()
        self.assertTrue(plots.is_dir())
        self.assertEqual(json.loads((self.dp / "results.json").read_text())["dG_bind"], -3.0)

    def test_enhanced_adds_the_extra_plots(self):
        calls, histograms = self.run_with_mocked_calculation(enhanced=True)
        plots = self.dp / "plots"
        self.assertEqual(calls["grid_path"], plots / "attach_windows_grid.jpg")
        self.assertEqual(calls["angle_grid_path"], plots / "attach_angle_grid.jpg")
        self.assertEqual(calls["convergence_path"], plots / "convergence_plot.jpg")
        self.assertTrue(calls["enhanced"])
        self.assertEqual(histograms.call_args.args[1], plots / "pull_window_histograms.jpg")

    def test_convergence_fractions_only_for_the_convergence_plot(self):
        for convergence_path, expected in ((None, [1.0]),
                                           ("c.jpg", apr_analysis.CONVERGENCE_FRACTIONS.tolist())):
            with self.subTest(convergence_path=convergence_path), \
                 patch.object(apr_analysis, "load_restraints", return_value=[]), \
                 patch.object(apr_analysis.analysis, "fe_calc", FakeFeCalc), \
                 self.assertRaises(StopAfterSetup):
                apr_analysis.calculate_free_energy(self.dp / "windows", convergence_path=convergence_path)
            self.assertEqual(FakeFeCalc.seen["fractions"], expected)


if __name__ == "__main__":
    unittest.main()
