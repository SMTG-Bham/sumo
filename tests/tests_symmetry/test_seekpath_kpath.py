import os
import unittest
import warnings

try:
    from importlib.resources import files as ilr_files
except ImportError:  # Python < 3.9
    from importlib_resources import files as ilr_files

from pymatgen.core.structure import Structure

from sumo.symmetry.seekpath_kpath import SeekpathKpath


class SeekpathKpathTestCase(unittest.TestCase):
    def setUp(self):
        ge_poscar = os.path.join(ilr_files("tests"), "data", "Ge", "POSCAR")

        with warnings.catch_warnings():  # Not interested in Pymatgen warnings
            warnings.simplefilter("ignore")
            self.ge_structure = Structure.from_file(ge_poscar)

    def test_seekpath_conversion(self):
        """Check path format conversion from seekpath"""
        in_path = [("A", "B"), ("B", "C"), ("D", "GAMMA")]
        in_point_coords = {
            "A": [1.0, 0.0, 0.0],
            "B": [0.0, 1.0, 0.0],
            "C": [0.0, 0.0, 1.0],
            "D": [0.5, 0.5, 0.0],
            "E": [1.0, 0.0, 1.0],
            "GAMMA": [0.0, 0.0, 0.0],
        }

        kpath = SeekpathKpath.kpath_from_seekpath(in_path, in_point_coords)
        self.assertEqual(kpath["path"], [["A", "B", "C"], ["D", r"\Gamma"]])
        self.assertEqual(
            kpath["kpoints"],
            {
                "A": [1.0, 0.0, 0.0],
                "B": [0.0, 1.0, 0.0],
                "C": [0.0, 0.0, 1.0],
                "D": [0.5, 0.5, 0.0],
                r"\Gamma": [0.0, 0.0, 0.0],
            },
        )

    def test_seekpath_kpath_load(self):
        """Check path generation from structure with seekpath"""
        kpath = SeekpathKpath(self.ge_structure)
        self.assertEqual(kpath.spg_symbol, "Fd-3m")
        self.assertEqual(kpath.lattice_type, "cubic")
        self.assertEqual(kpath.spg_number, 227)
        self.assertEqual(kpath.kpoints["L"], [0.5, 0.5, 0.5])
        self.assertEqual(
            kpath.path, [[r"\Gamma", "X", "U"], ["K", r"\Gamma", "L", "W", "X"]]
        )
        self.assertEqual(
            kpath.path_string, r"\Gamma -> X -> U | K -> \Gamma -> L -> W -> X"
        )


if __name__ == "__main__":
    unittest.main()
