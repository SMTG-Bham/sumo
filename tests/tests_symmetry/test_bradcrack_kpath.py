import pkg_resources
import unittest
import warnings

from os.path import join as path_join

from pymatgen.core.structure import Structure

from vaspy.symmetry.brad_crack_kpath import BradCrackKpath
from vaspy.symmetry.seekpath_kpath import SeekpathKpath


class BradCrackKpathTestCase(unittest.TestCase):

    def setUp(self):
        poscar = pkg_resources.resource_filename(
            __name__,
            path_join('..', 'data', 'Cs2SnI6', 'dos', 'vasprun.xml.gz'))
        with warnings.catch_warnings():  # Not interested in Pymatgen warnings
            warnings.simplefilter("ignore")
            self.cs_sn_i_structure = Structure.from_file(poscar)

    def test_bravais_assignment(self):
        """Check BradCrack space group logic"""
        self.assertEqual(BradCrackKpath._get_bravais_lattice(
            None, 'triclinic', None, None, None, None),
            'triclinic')
        self.assertEqual(BradCrackKpath._get_bravais_lattice(
            'Fd-3m', 'cubic', 5.66, 5.66, 5.66, 0),
            'cubic_f')

    def test_bradcrack_data(self):
        """Check BradCrack data loading"""
        trig_p_c = BradCrackKpath._get_bradcrack_data('trig_p_c')
        self.assertEqual(trig_p_c['path'][0][4], r'\Gamma')
        self.assertEqual(trig_p_c['kpoints']['K'][1], 0.667)

        mon_p_c = BradCrackKpath._get_bradcrack_data('mon_p_c')
        self.assertEqual(mon_p_c['path'][0][7], 'E0')
        self.assertEqual(mon_p_c['kpoints']['A0'][2], 0.)

    def test_bradcrack_seekpath_consistent(self):
        """Check BradCrack k-point locations are consistent with Seekpath"""
        bckp = BradCrackKpath(self.cs_sn_i_structure)
        spkp = SeekpathKpath(self.cs_sn_i_structure)

        for label, position in bckp.kpoints.items():
            self.assertIn(label, spkp.kpoints)
            self.assertEqual(position, spkp.kpoints[label])

if __name__ == '__main__':
    unittest.main()
