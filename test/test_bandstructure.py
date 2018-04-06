from os.path import join as path_join
import pkg_resources
import unittest

from pymatgen.core.structure import Structure

from vaspy.electronic_structure.bandstructure import SeekpathKpath

class SeekpathKpathTestCase(unittest.TestCase):

    def setUp(self):
        ge_poscar = pkg_resources.resource_filename(
            __name__, path_join('Ge', 'POSCAR')
        )
        self.ge_structure = Structure.from_file(ge_poscar)

    def test_seekpath_conversion(self):
        """Check path format conversion from seekpath"""
        in_path = [('A', 'B'), ('B', 'C'), ('D', 'GAMMA')]
        in_point_coords = {'A': [1., 0., 0.], 'B': [0., 1., 0.],
                           'C': [0., 0., 1.], 'D': [0.5, 0.5, 0.],
                           'E': [1., 0., 1.], 'GAMMA': [0., 0., 0.]}

        kpath = SeekpathKpath.kpath_from_seekpath(in_path, in_point_coords)
        self.assertEqual(kpath['path'], [['A', 'B', 'C'], ['D', '\Gamma']])
        self.assertEqual(kpath['kpoints'],
                         {'A': [1., 0., 0.], 'B': [0., 1., 0.],
                          'C': [0., 0., 1.], 'D': [0.5, 0.5, 0.],
                          '\Gamma': [0., 0., 0.]})

    def test_seekpath_kpath_load(self):
        """Check path generation from structure with seekpath"""
        kpath = SeekpathKpath(self.ge_structure)
        self.assertEqual(kpath.spg_symbol, 'Fd-3m')
        self.assertEqual(kpath.lattice_type, 'cubic')
        self.assertEqual(kpath.spg_number, 227)
        self.assertEqual(kpath.kpoints['L'], [0.5, 0.5, 0.5])
        self.assertEqual(kpath.path,
            [['\Gamma', 'X', 'U'], ['K', '\Gamma', 'L', 'W', 'X']])
        self.assertEqual(kpath.path_string,
                         '\Gamma -> X -> U | K -> \Gamma -> L -> W -> X')

if __name__ == '__main__':
    unittest.main()
