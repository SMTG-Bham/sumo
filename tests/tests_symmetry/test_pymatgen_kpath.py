import pkg_resources
import unittest
import warnings

from os.path import join as path_join

from pymatgen.core.structure import Structure

from sumo.symmetry.pymatgen_kpath import PymatgenKpath
from sumo.symmetry.brad_crack_kpath import BradCrackKpath


class SeekpathKpathTestCase(unittest.TestCase):

    def setUp(self):
        zno_poscar = pkg_resources.resource_filename(
            __name__, path_join('..', 'data', 'ZnO', 'POSCAR'))

        hgs_poscar = pkg_resources.resource_filename(
            __name__, path_join('..', 'data', 'Ge', 'POSCAR'))

        with warnings.catch_warnings():  # Not interested in Pymatgen warnings
            warnings.simplefilter("ignore")
            self.zno_structure = Structure.from_file(zno_poscar)
            self.hgs_structure = Structure.from_file(hgs_poscar)

    def test_pymatgen_path(self):
        """Check kpoint sequence from pymatgen: Ge example"""
        kpath = PymatgenKpath(self.zno_structure)
        self.assertEqual(kpath.path,
                         [[r'\Gamma', 'M', 'K', '\Gamma', 'A', 'L', 'H', 'A'],
                          ['L', 'M'], ['K', 'H']])

    def test_pymatgen_points(self):
        """Check special points agree between pymatgen and Bradley-Cracknell"""
        kpath_pymatgen = PymatgenKpath(self.hgs_structure)
        kpath_bradcrack = BradCrackKpath(self.hgs_structure)

        # Bradcrack kpoints should be a subset of pymatgen kpoints
        for label, position in kpath_bradcrack.kpoints.items():
            self.assertIn(label, kpath_pymatgen.kpoints)
            self.assertEqual(list(position),
                             list(kpath_pymatgen.kpoints[label]))
