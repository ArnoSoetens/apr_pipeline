import sys
import tempfile
import unittest
from pathlib import Path

tests_dir = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(tests_dir / "apr_pipeline"))
from datapoints import (
    COMPLEX_PDB,
    DONE,
    GUEST_FRCMOD,
    GUEST_MOL2,
    GUEST_PDB,
    TOPOLOGY,
    WINDOWS,
    all_data_points,
    next_stage,
    placement_inputs,
)


class SetupStageTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.data_dir = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def make(self, name, *files):
        dp = self.data_dir / name
        dp.mkdir()
        for filename in files:
            (dp / filename).write_text("x\n")
        return dp

    def build_windows(self, dp):
        for window in WINDOWS:
            (dp / "windows" / window).mkdir(parents=True)
            (dp / "windows" / window / TOPOLOGY).write_text("x\n")

    def test_smiles_only_needs_structure(self):
        self.assertEqual(next_stage(self.make("iso_00001", "iso_00001.smi")), "structure")

    def test_any_raw_pdb_needs_placement(self):
        self.assertEqual(next_stage(self.make("PFBA_m_down", "pfba_m_down.pdb")), "placement")

    def test_raw_pdb_wins_over_smiles(self):
        dp = self.make("iso_00001", "iso_00001.smi", "iso_00001.pdb")
        self.assertEqual(next_stage(dp), "placement")

    def test_failed_placement_still_needs_placement(self):
        # a clash leaves the rejected guest pose behind, but no complex PDB
        dp = self.make("PFBA_m_up", "pfba_m_up.pdb", GUEST_PDB)
        self.assertEqual(next_stage(dp), "placement")
        self.assertEqual(placement_inputs(dp), [dp / "pfba_m_up.pdb"])

    def test_lone_guest_pdb_is_a_placement_input(self):
        dp = self.make("custom", GUEST_PDB)
        self.assertEqual(next_stage(dp), "placement")
        self.assertEqual(placement_inputs(dp), [dp / GUEST_PDB])

    def test_complex_and_sqm_pdbs_are_never_inputs(self):
        dp = self.make("odd", COMPLEX_PDB, "sqm.pdb")
        self.assertEqual(placement_inputs(dp), [])
        self.assertIsNone(next_stage(dp))

    def test_placed_guest_needs_parametrize(self):
        dp = self.make("PFBA_m_down", "pfba_m_down.pdb", GUEST_PDB, COMPLEX_PDB, "sqm.pdb")
        self.assertEqual(next_stage(dp), "parametrize")

    def test_half_parametrized_guest_needs_parametrize(self):
        dp = self.make("half", GUEST_PDB, COMPLEX_PDB, GUEST_MOL2)
        self.assertEqual(next_stage(dp), "parametrize")

    def test_parametrized_guest_needs_apr_setup(self):
        dp = self.make("ready", GUEST_PDB, COMPLEX_PDB, GUEST_MOL2, GUEST_FRCMOD)
        self.assertEqual(next_stage(dp), "apr_setup")
        (dp / "windows" / WINDOWS[0]).mkdir(parents=True)   # setup started, not finished
        self.assertEqual(next_stage(dp), "apr_setup")

    def test_built_windows_are_done(self):
        dp = self.make("built", GUEST_PDB, COMPLEX_PDB, GUEST_MOL2, GUEST_FRCMOD)
        self.build_windows(dp)
        self.assertEqual(next_stage(dp), DONE)

    def test_empty_folder_has_no_stage(self):
        self.assertIsNone(next_stage(self.make("empty")))

    def test_all_data_points_finds_every_recognized_folder(self):
        smiles = self.make("iso_00001", "iso_00001.smi")
        raw = self.make("PFBA_m_down", "pfba_m_down.pdb")
        built = self.make("built", GUEST_PDB, COMPLEX_PDB, GUEST_MOL2, GUEST_FRCMOD)
        self.build_windows(built)
        self.make("empty")
        (self.data_dir / "BCD.pdb").write_text("host\n")   # files in the data dir are not data points

        self.assertEqual(all_data_points(data_dir=self.data_dir), [raw, built, smiles])
        self.assertEqual(all_data_points(only=["built"], data_dir=self.data_dir), [built])

    def test_all_data_points_without_data_dir(self):
        self.assertEqual(all_data_points(data_dir=self.data_dir / "missing"), [])


if __name__ == "__main__":
    unittest.main()
