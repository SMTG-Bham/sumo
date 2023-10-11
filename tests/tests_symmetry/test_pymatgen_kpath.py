import os
import unittest
import warnings

try:
    from importlib.resources import files as ilr_files
except ImportError:  # Python < 3.9
    from importlib_resources import files as ilr_files

from pymatgen.core.structure import Structure

from sumo.symmetry.brad_crack_kpath import BradCrackKpath
from sumo.symmetry.pymatgen_kpath import PymatgenKpath


class SeekpathKpathTestCase(unittest.TestCase):
    def setUp(self):
        zno_poscar = os.path.join(ilr_files("tests"), "data", "ZnO", "POSCAR")
        hgs_poscar = os.path.join(ilr_files("tests"), "data", "Ge", "POSCAR")

        with warnings.catch_warnings():  # Not interested in Pymatgen warnings
            warnings.simplefilter("ignore")
            self.zno_structure = Structure.from_file(zno_poscar)
            self.hgs_structure = Structure.from_file(hgs_poscar)

    def test_pymatgen_path(self):
        """Check kpoint sequence from pymatgen: Ge example"""
        kpath = PymatgenKpath(self.zno_structure)
        self.assertEqual(
            kpath.path,
            [
                [r"\Gamma", "M", "K", r"\Gamma", "A", "L", "H", "A"],
                ["L", "M"],
                ["K", "H"],
            ],
        )

    def test_pymatgen_points(self):
        """Check special points agree between pymatgen and Bradley-Cracknell"""
        kpath_pymatgen = PymatgenKpath(self.hgs_structure)
        kpath_bradcrack = BradCrackKpath(self.hgs_structure)

        # Bradcrack kpoints should be a subset of pymatgen kpoints
        for label, position in kpath_bradcrack.kpoints.items():
            self.assertIn(label, kpath_pymatgen.kpoints)
            self.assertEqual(list(position), list(kpath_pymatgen.kpoints[label]))
