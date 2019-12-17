import unittest
from os.path import join as path_join
from pkg_resources import resource_filename
import json
import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.electronic_structure.core import Spin
from sumo.io.castep import (read_bands_header,
                            read_bands_eigenvalues,
                            labels_from_cell,
                            band_structure)

_ry_to_ev = 13.605693009

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

        self.ref_labels = {'\Gamma': (0.0, 0.0, 0.0),
                           'W': (0.5, 0.25, 0.75),
                           'L': (0.5, 0.5, 0.5),
                           'X': (0.5, 0.0, 0.5),
                           'K': (0.375, 0.375, 0.75)}

    def test_castep_bands_read_header(self):
        header = read_bands_header(self.si_bands)
        with open(self.si_header_ref, 'r') as f:
            ref_header = json.load(f)
        self.assertEquals(header, ref_header)

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
        self.assertEquals(labels,
                          self.ref_labels)

    def test_castep_cell_read_labels_alt_spelling(self):
        # CASTEP input is case-insensitive anf allows BS_KPOINTS_LIST
        # as well as BS_KPOINT_LIST                            ^
        labels = labels_from_cell(self.si_cell_alt)
        self.assertEquals(labels,
                          self.ref_labels)

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
        self.assertEquals(header, ref_header)

    
