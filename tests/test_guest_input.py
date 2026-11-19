import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))
import guest_input
from guest_input import (
    GUEST_INPUT,
    anchors,
    net_charge,
    read_guest_input,
    save_anchors,
    update_guest_input,
)


class GuestInputTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.dp = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def write(self, values):
        (self.dp / GUEST_INPUT).write_text(json.dumps(values))

    def test_missing_file_gives_defaults(self):
        self.assertEqual(read_guest_input(self.dp),
                         {"smiles": None, "net_charge": None, "anchors": None})

    def test_update_keeps_other_values(self):
        update_guest_input(self.dp, smiles="CC", net_charge=0)
        update_guest_input(self.dp, anchors=[1, 2])
        self.assertEqual(json.loads((self.dp / GUEST_INPUT).read_text()),
                         {"smiles": "CC", "net_charge": 0, "anchors": [1, 2]})

    def test_update_rejects_unknown_keys(self):
        with self.assertRaises(KeyError):
            update_guest_input(self.dp, charge=-1)

    def test_invalid_json_names_the_file(self):
        (self.dp / GUEST_INPUT).write_text("{net_charge: -1")
        with self.assertRaisesRegex(ValueError, GUEST_INPUT):
            read_guest_input(self.dp)

    def test_auto_charge_comes_from_the_file(self):
        self.write({"net_charge": -1})
        with patch.dict(guest_input.CFG["parametrize"], {"net_charge": "auto"}):
            self.assertEqual(net_charge(self.dp), -1)

    def test_auto_charge_missing_asks_to_fill_in_the_file(self):
        update_guest_input(self.dp, anchors=[1, 5])     # e.g. created by placement of a PDB guest
        with patch.dict(guest_input.CFG["parametrize"], {"net_charge": "auto"}):
            with self.assertRaisesRegex(ValueError, "set net_charge"):
                net_charge(self.dp)
        for bad in ("-1", 1.5, True):
            self.write({"net_charge": bad})
            with patch.dict(guest_input.CFG["parametrize"], {"net_charge": "auto"}):
                with self.assertRaises(ValueError):
                    net_charge(self.dp)

    def test_config_integer_wins_over_the_file(self):
        self.write({"net_charge": -1})
        with patch.dict(guest_input.CFG["parametrize"], {"net_charge": 0}):
            self.assertEqual(net_charge(self.dp), 0)

    def test_anchors_are_stored_as_pdb_serials(self):
        save_anchors(self.dp, 0, 4)
        self.assertEqual(read_guest_input(self.dp)["anchors"], [1, 5])
        self.assertEqual(anchors(self.dp), (0, 4))

    def test_no_anchors_yet(self):
        self.assertIsNone(anchors(self.dp))

    def test_bad_anchors_are_rejected(self):
        for bad in ([1], [0, 5], [3, 3], ["C1", "C4"], 5):
            self.write({"anchors": bad})
            with self.assertRaisesRegex(ValueError, "column 2"):
                anchors(self.dp)


class OutputPdbConfigTests(unittest.TestCase):
    def test_output_pdb_must_be_a_plain_pdb_name(self):
        import copy
        import apr_config
        for bad in ("UNL", "sub/UNL.pdb", "BCD-UNL.pdb", "sqm.pdb"):
            cfg = copy.deepcopy(apr_config.CFG)
            cfg["structure"]["output_pdb"] = bad
            with self.subTest(bad=bad), self.assertRaisesRegex(ValueError, "structure.output_pdb"):
                apr_config._validate(cfg)
        cfg = copy.deepcopy(apr_config.CFG)
        cfg["structure"]["output_pdb"] = "raw_guest.pdb"
        apr_config._validate(cfg)


if __name__ == "__main__":
    unittest.main()
