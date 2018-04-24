# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Module containing class for generating Bradley and Cracknell k-point paths.
"""

import pkg_resources
from json import load as load_json
import numpy as np

from sumo.symmetry import Kpath


class BradCrackKpath(Kpath):
    r"""Class to generate k-points along paths from Bradley and Cracknell.

    The paths used are based on Brillouin zones depicted in reference [brad]_.

    .. [brad] "The Mathematical Theory of Symmetry in Solids", C. J. Bradley
              and A. P. Cracknell, Clarendon Press, 1972.

    These paths represent only one particular route through the Brillouin
    zone and do not cover every possible path (though they do visit every
    high-symmetry k-point at least once).

    These paths should be used with the correct primitive structure. This
    structure can be accessed using the ``prim`` attribute. Compliance
    between the provided structure and standardised structure checked using the
    ``correct_structure()`` method.

    Args:
        structure (:obj:`~pymatgen.core.structure.Structure`): The structure.
        symprec (:obj:`float`, optional): The tolerance for determining the
            crystal symmetry.
        spg (:obj:`~pymatgen.symmetry.groups.SpaceGroup`, optional): Space
            group used to override the symmetry determined by spglib. This is
            not recommended and only provided for testing purposes.

    Attributes:
        prim (:obj:`~pymatgen.core.structure.Structure`): The standardised
            primitive cell structure for the generated k-point path.
        conv (:obj:`~pymatgen.core.structure.Structure`): The standardised
            conventional cell structure.
    """

    def __init__(self, structure, symprec=1e-3, spg=None):
        Kpath.__init__(self, structure, symprec=symprec)

        angles = self.conv.lattice.angles
        unique = angles.index(min(angles, key=angles.count))
        a = self.conv.lattice.abc[0]
        b = self.conv.lattice.abc[1]
        c = self.conv.lattice.abc[2]

        if spg:
            spg_symbol = spg.symbol
            lattice_type = self.get_lattice_type(spg.int_number)
        else:
            spg_symbol = self.spg_symbol
            lattice_type = self.lattice_type

        bravais = self._get_bravais_lattice(spg_symbol, lattice_type,
                                            a, b, c, unique)
        self._kpath = self._get_bradcrack_data(bravais)

    @staticmethod
    def _get_bradcrack_data(bravais):
        r"""Read Bradley--Cracknell k-points path from data file

        Args:
            bravais (str): Lattice code including orientation e.g. 'trig_p_c'

        Returns:
            dict: kpoint path and special point locations, formatted as e.g.::

              {'kpoints': {'\Gamma': [0., 0., 0.], 'X': [0., 0.5, 0.], ...},
               'path': [['\Gamma', 'X', ..., 'P'], ['H', 'N', ...]]}

        """
        json_file = pkg_resources.resource_filename(__name__, 'bradcrack.json')
        with open(json_file, 'r') as f:
            bradcrack_data = load_json(f)
            return bradcrack_data[bravais]

    @staticmethod
    def _get_bravais_lattice(spg_symbol, lattice_type, a, b, c, unique):
        """Get Bravais lattice symbol from symmetry data"""

        if lattice_type == 'triclinic':
            return('triclinic')

        elif lattice_type == 'monoclinic':
            if 'P' in spg_symbol:
                if unique == 0:
                    return('mon_p_a')
                elif unique == 1:
                    return('mon_p_b')
                elif unique == 2:
                    return('mon_p_c')

            elif 'C' in spg_symbol:
                if unique == 0:
                    return('mon_c_a')
                elif unique == 1:
                    return('mon_c_b')
                elif unique == 2:
                    return('mon_c_c')

        elif lattice_type == 'orthorhombic':
            if 'P' in spg_symbol:
                return('orth_p')

            elif 'A' in spg_symbol or 'C' in spg_symbol:
                if a > b:
                    return('orth_c_a')
                elif b > a:
                    return('orth_c_b')

            elif 'F' in spg_symbol:
                if (1/a**2 < 1/b**2 + 1/c**2 and 1/b**2 < 1/c**2 + 1/a**2 and
                        1/c**2 < 1/a**2 + 1/b**2):
                    return('orth_f_1')
                elif 1/c**2 > 1/a**2 + 1/b**2:
                    return('orth_f_2')
                elif 1/b**2 > 1/a**2 + 1/c**2:
                    return('orth_f_3')
                elif 1/a**2 > 1/c**2 + 1/b**2:
                    return('orth_f_4')

            elif 'I' in spg_symbol:
                if a > b and a > c:
                    return('orth_i_a')
                elif b > a and b > c:
                    return('orth_i_b')
                elif c > a and c > b:
                    return('orth_i_c')

        elif lattice_type == 'tetragonal':
            if 'P' in spg_symbol:
                return('tet_p')

            elif 'I' in spg_symbol:
                if a > c:
                    return('tet_i_a')
                else:
                    return('tet_i_c')

        elif (lattice_type == 'trigonal' or lattice_type == 'hexagonal'
                or lattice_type == 'rhombohedral'):
            if 'R' in spg_symbol:
                if a > np.sqrt(2) * c:
                    return('trig_r_a')
                else:
                    return('trig_r_c')

            elif 'P' in spg_symbol:
                if unique == 0:
                    return('trig_p_a')
                elif unique == 2:
                    return('trig_p_c')

        elif lattice_type == "cubic":
            if 'P' in spg_symbol:
                return('cubic_p')
            elif 'I' in spg_symbol:
                return('cubic_i')
            elif 'F' in spg_symbol:
                return('cubic_f')
