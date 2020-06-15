import unittest
from os.path import join as path_join
from pkg_resources import resource_filename
import json
from numpy.testing import assert_array_almost_equal

from pymatgen import Structure
from pymatgen.electronic_structure.core import Spin
from sumo.io.castep import (read_bands_header,
                            read_bands_eigenvalues,
                            labels_from_cell,
                            CastepCell,
                            CastepPhonon)

_ry_to_ev = 13.605693009


class CastepCellTestCase(unittest.TestCase):
    def setUp(self):
        self.si_cell = resource_filename(
            __name__,
            path_join('..', 'data', 'Si', 'Si2.cell'))
        self.si_cell_alt = resource_filename(
            __name__,
            path_join('..', 'data', 'Si', 'Si2-alt.cell'))
        self.zns_band_cell = resource_filename(
            __name__,
            path_join('..', 'data', 'ZnS', 'zns.cell'))
        self.zns_singlepoint_cell = resource_filename(
            __name__,
            path_join('..', 'data', 'ZnS', 'zns-sp.cell'))
        si_structure_file = resource_filename(
            __name__,
            path_join('..', 'data', 'Si', 'Si8.json'))
        self.si_structure = Structure.from_file(si_structure_file)

    def test_castep_cell_null_init(self):
        null_cell = CastepCell()
        self.assertEqual(null_cell.blocks, {})
        self.assertEqual(null_cell.tags, {})

        with self.assertRaises(ValueError):
            structure = null_cell.structure

    def test_castep_cell_from_singlepoint_file(self):
        cc = CastepCell.from_file(self.zns_singlepoint_cell)
        self.assertEqual(set(cc.blocks.keys()),
                         set(('lattice_cart', 'positions_abs', 'species_pot')))
        self.assertEqual({k: v.value for k, v in cc.tags.items()},
                         {'fix_all_cell': ['true'],
                          'fix_all_ions': ['true'],
                          'symmetry_generate': ['true'],
                          'kpoint_mp_grid': ['4', '4', '4'],
                          'snap_to_symmetry': ['true']})
        self.assertEqual(cc.blocks['species_pot'].values[1],
                         ['S', 'NCP'])
        self.assertEqual(cc.blocks['species_pot'].comments,
                         ['', ''])

        structure = cc.structure
        self.assertIsInstance(structure, Structure)
        assert_array_almost_equal(structure.lattice.matrix,
                                  [[0., 2.71, 2.71],
                                   [2.71, 0., 2.71],
                                   [2.71, 2.71, 0.]])

    def test_castep_cell_from_structure(self):
        cell = CastepCell.from_structure(self.si_structure)
        self.assertEqual(cell.blocks['lattice_cart'].values,
                         [['ang'],
                          ['5.43', '0.0', '0.0'],
                          ['0.0', '5.43', '0.0'],
                          ['0.0', '0.0', '5.43']])
        self.assertEqual(cell.blocks['positions_frac'].values[0],
                         ['Si', '0.0', '0.0', '0.0'])
        self.assertEqual(cell.blocks['positions_frac'].values[7],
                         ['Si', '0.75', '0.75', '0.25'])


class CastepBandStructureTestCaseNoSpin(unittest.TestCase):
    def setUp(self):
        self.si_bands = resource_filename(
            __name__,
            path_join('..', 'data', 'Si', 'Si2.bands'))
        self.si_cell = resource_filename(
            __name__,
            path_join('..', 'data', 'Si', 'Si2.cell'))
        self.si_cell_alt = resource_filename(
            __name__,
            path_join('..', 'data', 'Si', 'Si2-alt.cell'))
        self.si_header_ref = resource_filename(
            __name__,
            path_join('..', 'data', 'Si', 'Si2.bands_header.json'))

        self.ref_labels = {r'\Gamma': (0.0, 0.0, 0.0),
                           'W': (0.5, 0.25, 0.75),
                           'L': (0.5, 0.5, 0.5),
                           'X': (0.5, 0.0, 0.5),
                           'K': (0.375, 0.375, 0.75)}

    def test_castep_bands_read_header(self):
        header = read_bands_header(self.si_bands)
        with open(self.si_header_ref, 'r') as f:
            ref_header = json.load(f)
        self.assertEqual(header, ref_header)

    def test_castep_bands_read_eigenvalues(self):
        with open(self.si_header_ref, 'r') as f:
            ref_header = json.load(f)
        kpoints, weights, eigenvals = read_bands_eigenvalues(
            self.si_bands, ref_header)

        for i, k in enumerate([0.5,
                               0.36111111,
                               0.63888889]):
            self.assertAlmostEqual(kpoints[4][i], k)

        self.assertAlmostEqual(eigenvals[Spin.up][2, 4],
                               0.09500443 * _ry_to_ev * 2)

        for weight in weights:
            self.assertAlmostEqual(weight, 0.02272727)

    def test_castep_cell_read_labels(self):
        labels = labels_from_cell(self.si_cell)
        self.assertEqual(labels,
                         self.ref_labels)

    def test_castep_cell_read_labels_alt_spelling(self):
        # CASTEP input is case-insensitive anf allows BS_KPOINTS_LIST
        # as well as BS_KPOINT_LIST                            ^
        labels = labels_from_cell(self.si_cell_alt)
        self.assertEqual(labels,
                         self.ref_labels)


class CastepBandStructureTestCaseNickel(unittest.TestCase):
    def setUp(self):
        self.ni_cell = resource_filename(
            __name__,
            path_join('..', 'data', 'Ni', 'ni-band.cell'))

        self.ref_labels = {r'\Gamma': (0.0, 0.0, 0.0),
                           'L': (0.5, 0.5, 0.5),
                           'W': (0.5, 0.25, 0.75),
                           'X': (0.5, 0.0, 0.5)}

    def test_castep_cell_read_labels_from_list(self):
        # The Si example uses handwritten .cell files in line-mode.
        # This example using a kgen-written .cell file with all k-points listed
        labels = labels_from_cell(self.ni_cell)
        self.assertEqual(labels, self.ref_labels)


class CastepBandStructureTestCaseWithSpin(unittest.TestCase):
    def setUp(self):
        self.fe_bands = resource_filename(
            __name__,
            path_join('..', 'data', 'Fe', 'Fe.bands'))
        self.fe_cell = resource_filename(
            __name__,
            path_join('..', 'data', 'Fe', 'Fe.cell'))
        self.fe_header_ref = resource_filename(
            __name__,
            path_join('..', 'data', 'Fe', 'Fe.bands_header.json'))

    def test_castep_bands_read_header(self):
        header = read_bands_header(self.fe_bands)
        with open(self.fe_header_ref, 'r') as f:
            ref_header = json.load(f)
        self.assertEqual(header, ref_header)


# To generate reference file::
#
#   sumo-phonon-bandplot -f zns.phonon --units cm-1 --to-json zns_phonon.json
class CastepPhononTestCaseZincblende(unittest.TestCase):
    def setUp(self):
        self.zns_phonon = resource_filename(
            __name__,
            path_join('..', 'data', 'ZnS', 'zns.phonon'))

        self.zns_cell = resource_filename(
            __name__,
            path_join('..', 'data', 'ZnS', 'zns.cell'))
        self.zns_phonon_ref = resource_filename(
            __name__,
            path_join('..', 'data', 'ZnS', 'zns_phonon.json'))

    def test_castep_phonon_read_bands(self):
        castep_phonon = CastepPhonon.from_file(self.zns_phonon)
        castep_phonon.set_labels_from_file(self.zns_cell)

        bs = castep_phonon.get_band_structure()
        bs_dict = bs.as_dict()

        with open(self.zns_phonon_ref, 'r') as f:
            ref_dict = json.load(f)

        for key in bs_dict.keys():
            self.assertIn(key, ref_dict.keys())

        # Comparing the band frequencies and displacements is painfully slow;
        # compare other stuff then use numpy for band data
        for key in ('lattice_rec', 'qpoints', 'labels_dict', 'structure'):
            self.assertEqual(bs_dict[key], ref_dict[key])

        assert_array_almost_equal(bs_dict['bands'], ref_dict['bands'])
        assert_array_almost_equal(bs_dict['eigendisplacements']['real'],
                                  ref_dict['eigendisplacements']['real'])
        assert_array_almost_equal(bs_dict['eigendisplacements']['imag'],
                                  ref_dict['eigendisplacements']['imag'])

        # May be an idea to check has_nac if that logic ever gets sorted out
