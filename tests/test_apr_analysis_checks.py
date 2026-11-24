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
    _halves,
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


if __name__ == "__main__":
    unittest.main()
