import os
import unittest
import warnings

try:
    from importlib.resources import files as ilr_files
except ImportError:  # Python < 3.9
    from importlib_resources import files as ilr_files

from pymatgen.core.structure import Structure

from sumo.symmetry.brad_crack_kpath import BradCrackKpath
from sumo.symmetry.seekpath_kpath import SeekpathKpath


class BradCrackKpathTestCase(unittest.TestCase):
    def setUp(self):
        poscar = os.path.join(
            ilr_files("tests"), "data", "Cs2SnI6", "dos", "vasprun.xml.gz"
        )
        with warnings.catch_warnings():  # Not interested in Pymatgen warnings
            warnings.simplefilter("ignore")
            self.cs_sn_i_structure = Structure.from_file(poscar)

    def test_bravais_assignment(self):
        """Check BradCrack space group logic"""
        self.assertEqual(
            BradCrackKpath._get_bravais_lattice(
                None, "triclinic", None, None, None, None
            ),
            "triclinic",
        )
        self.assertEqual(
            BradCrackKpath._get_bravais_lattice("Fd-3m", "cubic", 5.66, 5.66, 5.66, 0),
            "cubic_f",
        )

    def test_bradcrack_data(self):
        """Check BradCrack data loading"""
        trig_p_c = BradCrackKpath._get_bradcrack_data("trig_p_c")
        self.assertEqual(trig_p_c["path"][0][4], r"\Gamma")
        self.assertEqual(trig_p_c["kpoints"]["K"][1], 0.667)

        mon_p_c = BradCrackKpath._get_bradcrack_data("mon_p_c")
        self.assertEqual(mon_p_c["path"][0][7], "E0")
        self.assertEqual(mon_p_c["kpoints"]["A0"][2], 0.0)

    def test_bradcrack_seekpath_consistent(self):
        """Check BradCrack k-point locations are consistent with Seekpath"""
        kpath_bradcrack = BradCrackKpath(self.cs_sn_i_structure)
        kpath_seekpath = SeekpathKpath(self.cs_sn_i_structure)

        for label, position in kpath_bradcrack.kpoints.items():
            self.assertIn(label, kpath_seekpath.kpoints)
            self.assertEqual(position, kpath_seekpath.kpoints[label])


if __name__ == "__main__":
    unittest.main()
