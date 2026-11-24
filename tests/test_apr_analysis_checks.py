"""The --enhanced diagnostics in apr_analysis: attach-window grid with halves and rim lines."""
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import parmed as pmd
from parmed.topologyobjects import Atom

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))
import apr_analysis
from apr_analysis import (
    _convergence_curves,
    _halves,
    _ti_integral,
    _ti_per_frame,
    convergence_plot,
    _read_cone_wall,
    attach_angle_grid,
    _rim_heights_from_structure,
    attach_window_grid,
)


def ring_structure(rim_heights=(11.0, 3.0), dm1_z=-2.0):
    """DM1 plus two horizontal rings of 7 oxygens (BCD atoms 1-7 and 8-14) at the given
    heights above DM1; the first group is the far rim, to check that the order doesn't matter."""
    structure = pmd.Structure()
    for group, height in enumerate(rim_heights):
        for i in range(7):
            angle = 2.0 * np.pi * i / 7
            atom = Atom(name=f"O{group * 7 + i + 1}", atomic_number=8, mass=16.0)
            atom.xx, atom.xy, atom.xz = 5.0 * np.cos(angle), 5.0 * np.sin(angle), dm1_z + height
            structure.add_atom(atom, "BCD", 1)
    dummy = Atom(name="DUM", atomic_number=6, mass=12.0)
    dummy.xx, dummy.xy, dummy.xz = 0.0, 0.0, dm1_z
    structure.add_atom(dummy, "DM1", 2)
    return structure


class FakeFe:
    """simulation_data and prepare_data, the two things attach_window_grid reads."""

    def __init__(self, windows, targets):
        self.simulation_data = {"attach": [[w] for w in windows]}
        self._targets = targets

    def prepare_data(self, phase):
        return None, None, None, None, None, [SimpleNamespace(magnitude=np.asarray(self._targets))], None


def vertical_lines(ax):
    """(x, linestyle, colour) of every axvline on ax."""
    out = []
    for line in ax.get_lines():
        x = line.get_xdata()
        if len(x) == 2 and x[0] == x[1]:
            out.append((float(x[0]), line.get_linestyle(), line.get_color()))
    return out


class HalvesTests(unittest.TestCase):
    def test_odd_frame_count_gives_the_extra_frame_to_the_second_half(self):
        first, second = _halves(np.arange(401))
        self.assertEqual((len(first), len(second)), (200, 201))
        self.assertEqual(first[-1] + 1, second[0])


class RimHeightTests(unittest.TestCase):
    def test_near_and_far_rim_from_dm1_whatever_the_group_order(self):
        structure = ring_structure(rim_heights=(11.0, 3.0))
        near, far = _rim_heights_from_structure(
            structure, rim_atom_groups=[list(range(1, 8)), list(range(8, 15))], host_prefix="BCD")
        self.assertAlmostEqual(near, 3.0, places=6)
        self.assertAlmostEqual(far, 11.0, places=6)

    def test_rim_heights_reads_the_a000_structure(self):
        loaded = []

        def fake_load(prmtop, rst7):
            loaded.append((prmtop, rst7))
            return ring_structure(rim_heights=(4.0, 9.5))

        with patch.object(apr_analysis.pmd, "load_file", side_effect=fake_load):
            near, far = apr_analysis.rim_heights(
                "/w", complex_prefix="X",
                rim_atom_groups=[list(range(1, 8)), list(range(8, 15))], host_prefix="BCD")
        self.assertEqual(loaded, [("/w/a000/X-dum.prmtop", "/w/a000/X-dum.rst7")])
        self.assertAlmostEqual(near, 4.0, places=6)
        self.assertAlmostEqual(far, 9.5, places=6)


class AttachGridTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        self.windows = [rng.normal(6.0, 1.0, 400) for _ in range(3)]
        self.targets = [6.0, 6.0, 6.0]
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name) / "grid.jpg"

    def tearDown(self):
        self.tmp.cleanup()

    def test_halves_target_and_rim_lines_in_every_panel(self):
        fig = attach_window_grid(FakeFe(self.windows, self.targets), self.out,
                                 labels=["a000", "a001", "p000"], rims=(3.0, 11.0))
        self.assertTrue(self.out.is_file())
        panels = fig.axes[:3]
        for ax in panels:
            lines = vertical_lines(ax)
            self.assertIn((6.0, "--", "red"), lines)
            dotted = sorted(x for x, style, _ in lines if style == ":")
            self.assertEqual(dotted, [3.0, 11.0])
            outline_colours = {p.get_edgecolor()[:3] for p in ax.patches if not p.get_fill()}
            self.assertEqual(len(outline_colours), 2)          # 1st and 2nd half
        self.assertEqual([ax.get_title() for ax in panels], ["a000", "a001", "p000 (λ=1)"])
        self.assertFalse(fig.axes[3].axison)                    # unused slot switched off

    def test_halves_are_plotted_separately(self):
        windows = [np.concatenate([np.full(200, 5.0), np.full(200, 7.0)])]
        fig = attach_window_grid(FakeFe(windows, [6.0]), self.out, labels=["a000"])
        ax = fig.axes[0]
        outlines = [p for p in ax.patches if not p.get_fill()]
        self.assertEqual(len(outlines), 2)
        means = []
        for patch in outlines:
            xy = patch.get_xy()
            tops = xy[xy[:, 1] > 0]
            means.append(tops[:, 0].mean())
        self.assertLess(means[0], 6.0)                          # 1st half sits at 5 Å
        self.assertGreater(means[1], 6.0)                       # 2nd half sits at 7 Å

    def test_no_rim_lines_without_rims(self):
        fig = attach_window_grid(FakeFe(self.windows, self.targets), self.out)
        for ax in fig.axes[:3]:
            self.assertFalse([x for x, style, _ in vertical_lines(ax) if style == ":"])
        self.assertEqual(fig.axes[2].get_title(), "a002")

    def test_legend_goes_in_the_first_unused_slot_or_above_the_panels(self):
        fig = attach_window_grid(FakeFe(self.windows, self.targets), self.out, rims=(3.0, 11.0))
        legend = fig.axes[3].get_legend()                       # 3 windows, 5 slots
        self.assertEqual([t.get_text() for t in legend.get_texts()],
                         ["1st half", "2nd half", "target", "rim"])
        self.assertFalse(fig.legends)
        fig = attach_window_grid(FakeFe(self.windows[:2], [6.0, 6.0]), self.out, ncols=2)
        self.assertEqual(len(fig.legends), 1)                   # every slot used

    def test_label_count_must_match_the_windows(self):
        with self.assertRaisesRegex(ValueError, "2 labels for 3 attach windows"):
            attach_window_grid(FakeFe(self.windows, self.targets), self.out, labels=["a000", "a001"])


# the a000 lines of wall_test/iso_00001 (distance, angle and torsion restraints at lambda = 0)
# plus a cone wall as apr_setup writes it
DISANG = """\
&rst iat= 156,148,          r1=    0.00000, r2=    6.00000, r3=    6.00000, r4=  999.00000, rk2=    0.00000, rk3=    0.00000,  &end
&rst iat= 157,156,148,      r1=    0.00000, r2=  180.00000, r3=  180.00000, r4=  180.00000, rk2=    0.00000, rk3=    0.00000,  &end
&rst iat= 156,148,149,      r1=    0.00000, r2=  180.00000, r3=  180.00000, r4=  180.00000, rk2=  100.00000, rk3=  100.00000,  &end
&rst iat= 157,156,75,       r1=    0.00000, r2=  141.06222, r3=  141.06222, r4=  180.00000, rk2=  100.00000, rk3=  100.00000,  &end
&rst iat= 156,148,          r1=    1.67671, r2=    2.67671, r3=   11.65314, r4=   12.65314, rk2=  100.00000, rk3=  100.00000,  &end
"""
CONE = "&rst iat= 157,156,148,      r1=  127.00000, r2=  132.00000, r3=  180.00000, r4=  180.00000, rk2=  100.00000, rk3=  100.00000,  &end\n"


class ConeWallTests(unittest.TestCase):
    def read(self, text):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "disang.rest"
            path.write_text(text)
            return _read_cone_wall(path)

    def test_cone_r2_is_found_between_the_angle_restraints(self):
        self.assertEqual(self.read(DISANG + CONE), 132.0)

    def test_no_cone_gives_none(self):
        self.assertIsNone(self.read(DISANG))      # guest angles (r2 == r3) and the distance wall


class FakeAngleFe:
    """Attach data as pAPRika has it: a distance restraint first, then an angle in radians."""

    def __init__(self, angles_deg, target_deg=180.0):
        from openff.units import unit
        self.unit = unit
        self.simulation_data = {"attach": [
            [np.full(len(a), 6.0) * unit.angstrom, np.radians(a) * unit.radian] for a in angles_deg]}
        self._targets = [np.full(len(angles_deg), 6.0) * unit.angstrom,
                         np.radians(np.full(len(angles_deg), target_deg)) * unit.radian]

    def prepare_data(self, phase):
        restraints = [SimpleNamespace(mask3=None, mask4=None),
                      SimpleNamespace(mask3=":UNL@C1", mask4=None)]
        return None, None, None, restraints, None, self._targets, None


class AttachAngleGridTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name) / "angles.jpg"

    def tearDown(self):
        self.tmp.cleanup()

    def test_angles_in_degrees_with_target_and_cone_line(self):
        angles = [np.linspace(150.0, 170.0, 100), np.linspace(160.0, 179.0, 100)]
        fig = attach_angle_grid(FakeAngleFe(angles), self.out, labels=["a000", "p000"], cone=132.0)
        self.assertTrue(self.out.is_file())
        for ax, window in zip(fig.axes[:2], angles):
            lines = vertical_lines(ax)
            self.assertIn((180.0, "--", "red"), [(round(x, 6), s, c) for x, s, c in lines])
            self.assertEqual([x for x, style, _ in lines if style == ":"], [132.0])
            edges = [p.get_xy()[:, 0] for p in ax.patches if not p.get_fill()]
            self.assertAlmostEqual(min(e.min() for e in edges), window.min(), places=6)
            self.assertAlmostEqual(max(e.max() for e in edges), window.max(), places=6)
        self.assertEqual(fig.axes[1].get_title(), "p000 (λ=1)")
        self.assertEqual(fig.axes[1].get_xlabel(), "DM2-DM1-G1 angle (degrees)")

    def test_without_cone_no_dotted_line(self):
        fig = attach_angle_grid(FakeAngleFe([np.linspace(150.0, 170.0, 100)]), self.out)
        self.assertFalse([x for x, style, _ in vertical_lines(fig.axes[0]) if style == ":"])

    def test_refuses_when_there_is_no_angle_restraint(self):
        fe = FakeAngleFe([np.linspace(150.0, 170.0, 100)])
        fe.prepare_data = lambda phase: (None, None, None, [SimpleNamespace(mask3=None, mask4=None)],
                                         None, None, None)
        with self.assertRaisesRegex(ValueError, "no angle restraint"):
            attach_angle_grid(fe, self.out)


class AxisLabelTests(unittest.TestCase):
    def test_x_label_and_ticks_only_on_panels_without_a_panel_below(self):
        with tempfile.TemporaryDirectory() as tmp:
            windows = [np.random.default_rng(i).normal(6.0, 1.0, 50) for i in range(7)]
            fig = attach_window_grid(FakeFe(windows, [6.0] * 7), Path(tmp) / "g.jpg")
        labelled = [k for k, ax in enumerate(fig.axes[:7]) if ax.get_xlabel()]
        self.assertEqual(labelled, [2, 3, 4, 5, 6])                # 7 windows in rows of 5


def prepared(values, force_constants, targets, units=("angstrom", "kcal/mol/angstrom**2"),
             masks=(None, None)):
    """prepare_data's tuple for one restraint: values per window, force constant and target per
    window, all as pint quantities like pAPRika has them."""
    from openff.units import unit
    restraint = SimpleNamespace(mask3=masks[0], mask4=masks[1])
    return (len(values), [len(v) for v in values], max(len(v) for v in values), [restraint],
            [np.asarray(force_constants, dtype=float) * unit(units[1])],
            [np.asarray(targets, dtype=float) * unit(units[0])],
            [[np.asarray(v, dtype=float) * unit(units[0])] for v in values])


class SequenceRng:
    """Stands in for numpy's Generator: hands out fixed samples."""

    def __init__(self, samples):
        self.samples = samples

    def normal(self, means, sems, size):
        return self.samples


class TiPerFrameTests(unittest.TestCase):
    def test_attach_du_is_full_k_times_squared_deviation(self):
        values = [[5.0, 7.0], [6.5], [8.0]]
        dl, dU = _ti_per_frame(prepared(values, [0.0, 5.0, 10.0], [6.0, 6.0, 6.0]), "attach")
        np.testing.assert_allclose(dl, [0.0, 0.5, 1.0])
        np.testing.assert_allclose(dU[0], [10.0, 10.0])
        np.testing.assert_allclose(dU[1], [2.5])
        np.testing.assert_allclose(dU[2], [40.0])

    def test_pull_du_is_two_k_times_deviation_and_x_is_the_target(self):
        values = [[6.5], [6.0], [9.0]]
        dl, dU = _ti_per_frame(prepared(values, [10.0] * 3, [6.0, 7.0, 8.0]), "pull")
        np.testing.assert_allclose(dl, [6.0, 7.0, 8.0])
        np.testing.assert_allclose(np.concatenate(dU), [10.0, -20.0, 20.0])

    def test_dihedrals_are_wrapped_around_the_target_without_touching_the_data(self):
        values = [np.radians([-175.0])]
        data = prepared(values, [100.0], np.radians([170.0]), units=("radian", "kcal/mol/radian**2"),
                        masks=(":A", ":B"))
        _, dU = _ti_per_frame(data, "attach")
        np.testing.assert_allclose(dU[0], [100.0 * np.radians(15.0) ** 2])
        np.testing.assert_allclose(data[6][0][0].magnitude, np.radians([-175.0]))


class TiIntegralTests(unittest.TestCase):
    def test_same_number_as_paprika_integrate_bootstraps_full_matrix(self):
        from paprika.analysis.bootstrap import integrate_bootstraps
        dl = np.array([0.0, 0.002, 0.01, 0.1, 0.4, 1.0])
        rng = np.random.default_rng(3)
        samples = rng.normal([5.0, 4.0, 3.0, 2.5, 2.0, 1.8], 0.2, size=(50, 6))
        x_intp = np.concatenate([np.linspace(dl[k - 1], dl[k], 100, endpoint=False)
                                 for k in range(1, 6)] + [[dl[-1]]])
        mean, sem = integrate_bootstraps(dl, samples, x_intp=x_intp, matrix="full")
        ours = _ti_integral(dl, np.zeros(6), np.zeros(6), SequenceRng(samples), cycles=50)
        self.assertAlmostEqual(ours[0], mean[0, -1], places=12)
        self.assertAlmostEqual(ours[1], sem[0, -1], places=12)


class ConvergenceCurveTests(unittest.TestCase):
    dl = np.array([0.0, 0.25, 0.5, 0.75, 1.0])

    def windows(self, n, seed=0, drift=0.0):
        rng = np.random.default_rng(seed)
        out = []
        for k in range(len(self.dl)):
            window = rng.normal(10.0 - 5.0 * self.dl[k], 1.0, n)
            window[n // 2:] += drift
            out.append(window)
        return out

    def test_forward_and_reverse_agree_with_all_frames(self):
        # 360 frames: get_nearest_max(360) = 360, so both directions block the same frames
        curves = _convergence_curves(self.dl, self.windows(360), "attach", cycles=200)
        self.assertAlmostEqual(curves["forward"]["fe"][-1], curves["reverse"]["fe"][-1], places=12)
        self.assertAlmostEqual(curves["forward"]["sem"][-1], curves["reverse"]["sem"][-1], places=12)
        self.assertEqual(curves["forward"]["frames"][-1], 5 * 360)
        self.assertEqual(curves["forward"]["frames"][0], 5 * 36)

    def test_independent_samples_give_sem_close_to_n_to_the_minus_half(self):
        curves = _convergence_curves(self.dl, self.windows(4000, seed=1), "attach", cycles=200)
        f = curves["forward"]
        slope = np.polyfit(np.log(f["frames"][2:]), np.log(f["sem"][2:]), 1)[0]
        self.assertGreater(slope, -0.65)
        self.assertLess(slope, -0.35)

    def test_drift_separates_forward_from_reverse(self):
        curves = _convergence_curves(self.dl, self.windows(2000, seed=2, drift=1.0), "attach",
                                     cycles=200)
        i = 4                                                       # f = 0.5
        gap = curves["reverse"]["fe"][i] - curves["forward"]["fe"][i]
        combined = np.hypot(curves["forward"]["sem"][i], curves["reverse"]["sem"][i])
        self.assertAlmostEqual(gap, 1.0, delta=0.15)                # the drift, integrated over 0-1
        self.assertGreater(gap / combined, 5)

    def test_pull_gets_the_minus_sign(self):
        windows = [np.full(100, 2.0) for _ in self.dl]
        curves = _convergence_curves(self.dl, windows, "pull", cycles=20)
        # like pAPRika, the integral stops 1/100 of the last interval (0.25) before the end
        self.assertAlmostEqual(curves["forward"]["fe"][-1], -2.0 * (1.0 - 0.0025), places=9)


class ConvergencePlotTests(unittest.TestCase):
    def test_four_panels_with_log_axes_on_the_right(self):
        rng = np.random.default_rng(4)
        attach = prepared([rng.normal(6.0, 1.0 / (1 + k), 200) for k in range(4)],
                          [0.0, 2.0, 5.0, 10.0], [6.0] * 4)
        pull = prepared([rng.normal(6.0 + k, 0.3, 200) for k in range(4)], [10.0] * 4,
                        [6.0, 7.0, 8.0, 9.0])
        fe = SimpleNamespace(prepare_data=lambda phase: attach if phase == "attach" else pull)
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "convergence.jpg"
            fig, curves = convergence_plot(fe, out, cycles=20)
            self.assertTrue(out.is_file())
        self.assertEqual(len(fig.axes), 4)
        self.assertEqual([ax.get_xscale() for ax in fig.axes], ["linear", "log", "linear", "log"])
        self.assertEqual([ax.get_yscale() for ax in fig.axes], ["linear", "log", "linear", "log"])
        self.assertEqual(set(curves), {"attach", "pull"})
        self.assertEqual([t.get_text() for t in fig.axes[0].get_legend().get_texts()],
                         ["forward (first f)", "reverse (last f)"])


if __name__ == "__main__":
    unittest.main()
