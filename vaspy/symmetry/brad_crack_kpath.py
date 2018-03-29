# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Module containing class for generating Bradley and Cracknell k-point paths.
"""

import numpy as np

from vaspy.symmetry import Kpath


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

        if lattice_type == 'triclinic':
            self._kpath = self._triclinic()

        elif lattice_type == 'monoclinic':
            if 'P' in spg_symbol:
                if unique == 0:
                    self._kpath = self._mon_p_a()
                elif unique == 1:
                    self._kpath = self._mon_p_b()
                elif unique == 2:
                    self._kpath = self._mon_p_c()

            elif 'C' in spg_symbol:
                if unique == 0:
                    self._kpath = self._mon_c_a()
                elif unique == 1:
                    self._kpath = self._mon_c_b()
                elif unique == 2:
                    self._kpath = self._mon_c_c()

        elif lattice_type == 'orthorhombic':
            if 'P' in spg_symbol:
                self._kpath = self._orth_p()

            elif 'C' in spg_symbol:
                if a > b:
                    self._kpath = self._orth_c_a()
                elif b > a:
                    self._kpath = self._orth_c_b()

            elif 'F' in spg_symbol:
                if (1/a**2 < 1/b**2 + 1/c**2 and 1/b**2 < 1/c**2 + 1/a**2 and
                        1/c**2 < 1/a**2 + 1/b**2):
                    self._kpath = self._orth_f_1()
                elif 1/c**2 > 1/a**2 + 1/b**2:
                    self._kpath = self._orth_f_2()
                elif 1/b**2 > 1/a**2 + 1/c**2:
                    self._kpath = self._orth_f_3()
                elif 1/a**2 > 1/c**2 + 1/b**2:
                    self._kpath = self._orth_f_4()

            elif 'I' in spg_symbol:
                if a > b and a > c:
                    self._kpath = self._orth_i_a()
                elif b > a and b > c:
                    self._kpath = self._orth_i_b()
                elif c > a and c > b:
                    self._kpath = self._orth_i_c()

        elif lattice_type == 'tetragonal':
            if 'P' in spg_symbol:
                self._kpath = self._tet_p()

            elif 'I' in spg_symbol:
                if a > c:
                    self._kpath = self._tet_i_a()
                else:
                    self._kpath = self._tet_i_c()

        elif (lattice_type == 'trigonal' or lattice_type == 'hexagonal'
                or lattice_type == 'rhombohedral'):
            if 'R' in spg_symbol:
                if a > np.sqrt(2) * c:
                    self._kpath = self._trig_r_a()
                else:
                    self._kpath = self._trig_r_c()

            elif 'P' in spg_symbol:
                if unique == 0:
                    self._kpath = self._trig_p_a()
                elif unique == 2:
                    self._kpath = self._trig_p_c()

        elif lattice_type == "cubic":
            if 'P' in spg_symbol:
                self._kpath = self._cubic_p()
            elif 'I' in spg_symbol:
                self._kpath = self._cubic_i()
            elif 'F' in spg_symbol:
                self._kpath = self._cubic_f()

    def _triclinic(self):
        path = [[r"\Gamma", "Z", "T", "Y", r"\Gamma", "X", "V", "R", "U"]]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5]),
                   'T': np.array([0.0, 0.5, 0.5]),
                   'Y': np.array([0.0, 0.5, 0.0]),
                   'X': np.array([0.5, 0.0, 0.0]),
                   'V': np.array([0.5, 0.5, 0.0]),
                   'R': np.array([0.5, 0.5, 0.5]),
                   'U': np.array([0.5, 0.0, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def _tet_p(self):
        path = [[r"\Gamma", "X", "M", r"\Gamma", "Z", "R", "A", "Z"]]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'X': np.array([0.0, 0.5, 0.0]),
                   'M': np.array([0.5, 0.5, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5]),
                   'R': np.array([0.0, 0.5, 0.5]),
                   'A': np.array([0.5, 0.5, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def _tet_i_a(self):
        path = [[r"\Gamma", "X", "P", "N", r"\Gamma", "Z"]]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'X': np.array([0.0, 0.0, 0.5]),
                   'P': np.array([0.25, 0.25, 0.25]),
                   'N': np.array([0.0, 0.5, 0.0]),
                   'Z': np.array([-0.5, 0.5, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def _tet_i_c(self):
        path = [[r"\Gamma", "X", "P", "N", r"\Gamma", "Z"]]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'X': np.array([0.0, 0.0, 0.5]),
                   'P': np.array([0.25, 0.25, 0.25]),
                   'N': np.array([0.0, 0.5, 0.0]),
                   'Z': np.array([0.5, 0.5, -0.5])}
        return {'kpoints': kpoints, 'path': path}

    def _cubic_p(self):
        path = [[r"\Gamma", "M", "R", "X", r"\Gamma", "R"]]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'M': np.array([0.5, 0.5, 0.0]),
                   'R': np.array([0.5, 0.5, 0.5]),
                   'X': np.array([0.0, 0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def _cubic_f(self):
        path = [[r"\Gamma", "L", "W", "X", r"\Gamma"]]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'W': np.array([0.5, 0.25, 0.75]),
                   'X': np.array([0.5, 0.0, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def _cubic_i(self):
        path = [[r"\Gamma", "P", "N", r"\Gamma", "H", "P"], ["H", "N"]]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'P': np.array([0.25, 0.25, 0.25]),
                   'N': np.array([0.0, 0.0, 0.5]),
                   'H': np.array([0.5, -0.5, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def _trig_r_a(self):
        path = [[r'\Gamma', 'L', 'F', r'\Gamma', 'Z']]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'L': np.array([0.0, 0.5, 0.0]),
                   'F': np.array([0.5, 0.5, 0.0]),
                   'Z': np.array([0.5, 0.5, -0.5])}
        return {'kpoints': kpoints, 'path': path}

    def _trig_r_c(self):
        path = [[r'\Gamma', 'L', 'F', r'\Gamma', 'Z']]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'L': np.array([0.0, 0.5, 0.0]),
                   'F': np.array([0.5, 0.5, 0.0]),
                   'Z': np.array([0.5, 0.5, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def _trig_p_a(self):
        path = [[r'\Gamma', 'A', 'L', 'M', r'\Gamma', 'K', 'H', 'A']]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([0.5, 0.0, 0.0]),
                   'L': np.array([0.5, 0.5, 0.0]),
                   'M': np.array([0.0, 0.5, 0.0]),
                   'K': np.array([0.0, 0.333, 0.333]),
                   'H': np.array([0.5, 0.333, 0.333])}
        return {'kpoints': kpoints, 'path': path}

    def _trig_p_c(self):
        path = [[r'\Gamma', 'A', 'L', 'M', r'\Gamma', 'K', 'H', 'A']]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([0.0, 0.0, 0.5]),
                   'L': np.array([0.0, 0.5, 0.5]),
                   'M': np.array([0.0, 0.5, 0.0]),
                   'K': np.array([-0.333, 0.667, 0.0]),
                   'H': np.array([-0.333, 0.667, 0.6])}
        return {'kpoints': kpoints, 'path': path}

    def _orth_p(self):
        path = [[r'\Gamma', 'Z', 'T', 'Y', 'S', 'R', 'U', 'X', r'\Gamma', 'Y']]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5]),
                   'T': np.array([-0.5, 0.0, 0.5]),
                   'Y': np.array([-0.5, 0.0, 0.0]),
                   'S': np.array([-0.5, 0.5, 0.0]),
                   'R': np.array([-0.5, 0.5, 0.5]),
                   'U': np.array([0.0, 0.5, 0.5]),
                   'X': np.array([0.0, 0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def _orth_c_a(self):
        path = [['R', 'S', r'\Gamma', 'Z', 'T', 'Y', r'\Gamma']]
        kpoints = {'R': np.array([0.0, 0.5, 0.5]),
                   'S': np.array([0.0, 0.5, 0.0]),
                   r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5]),
                   'T': np.array([0.5, 0.5, 0.5]),
                   'Y': np.array([0.5, 0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def _orth_c_b(self):
        path = [['R', 'S', r'\Gamma', 'Z', 'T', 'Y', r'\Gamma']]
        kpoints = {'R': np.array([0.0, 0.5, 0.5]),
                   'S': np.array([0.0, 0.5, 0.0]),
                   r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5]),
                   'T': np.array([-0.5, 0.5, 0.5]),
                   'Y': np.array([-0.5, 0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def _orth_f_1(self):
        path = [[r'\Gamma', 'Y', 'X', 'Z', r'\Gamma', 'L']]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Y': np.array([0.0, -0.5, -0.5]),
                   'X': np.array([0.5, 0.0, 0.5]),
                   'Z': np.array([0.5, 0.5, 0.0]),
                   'L': np.array([0.5, 0.0, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def _orth_f_2(self):
        path = [[r'\Gamma', 'Y', 'X', 'Z', r'\Gamma', 'L']]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Y': np.array([0.0, -0.5, -0.5]),
                   'X': np.array([0.5, 0.0, 0.5]),
                   'Z': np.array([0.5, -0.5, 0.0]),
                   'L': np.array([0.5, 0.0, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def _orth_f_3(self):
        path = [[r'\Gamma', 'Y', 'X', 'Z', r'\Gamma', 'L']]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Y': np.array([1.0, 0.5, 0.5]),
                   'X': np.array([0.5, 0.0, 0.5]),
                   'Z': np.array([0.5, 0.5, 0.0]),
                   'L': np.array([0.5, 0.0, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def _orth_f_4(self):
        path = [[r'\Gamma', 'Y', 'X', 'Z', r'\Gamma', 'L']]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Y': np.array([0.0, -0.5, -0.5]),
                   'X': np.array([0.5, 0.0, -0.5]),
                   'Z': np.array([0.5, 0.5, 0.0]),
                   'L': np.array([0.5, 0.0, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def _orth_i_a(self):
        path = [['R', r'\Gamma', 'X', 'S', 'W', 'T']]
        kpoints = {'R': np.array([0.5, 0.0, 0.0]),
                   r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'X': np.array([0.5, -0.5, 0.5]),
                   'S': np.array([0.5, 0.0, -0.5]),
                   'W': np.array([0.75, -0.25, -0.25]),
                   'T': np.array([0.5, -0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def _orth_i_b(self):
        path = [['R', r'\Gamma', 'X', 'S', 'W', 'T']]
        kpoints = {'R': np.array([0.5, 0.0, 0.0]),
                   r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'X': np.array([0.5, -0.5, -0.5]),
                   'S': np.array([0.5, 0.0, -0.5]),
                   'W': np.array([0.75, -0.25, -0.25]),
                   'T': np.array([0.5, -0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def _orth_i_c(self):
        path = [['R', r'\Gamma', 'X', 'S', 'W', 'T']]
        kpoints = {'R': np.array([0.5, 0.0, 0.0]),
                   r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'X': np.array([0.5, 0.5, -0.5]),
                   'S': np.array([0.5, 0.0, -0.5]),
                   'W': np.array([0.75, -0.25, -0.25]),
                   'T': np.array([0.5, -0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def _mon_p_a(self):
        path = [[r'\Gamma', 'Z', 'C', 'Y', r'\Gamma', 'B', 'D', 'E0', 'A0',
                 'Y']]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Z': np.array([0.5, 0.0, 0.0]),
                   'C': np.array([0.5, 0.0, 0.5]),
                   'Y': np.array([0.0, 0.0, 0.5]),
                   'B': np.array([0.0, 0.5, 0.0]),
                   'D': np.array([0.5, 0.5, 0.0]),
                   'E0': np.array([0.5, 0.5, 0.5]),
                   'A0': np.array([0.0, 0.5, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def _mon_p_b(self):
        path = [[r'\Gamma', 'Z', 'C', 'Y', r'\Gamma', 'B', 'D', 'E0', 'A0',
                 'Y']]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Z': np.array([0.5, 0.0, 0.0]),
                   'C': np.array([0.5, 0.5, 0.0]),
                   'Y': np.array([0.5, 0.0, 0.0]),
                   'B': np.array([0.0, 0.0, 0.5]),
                   'D': np.array([0.0, 0.5, 0.5]),
                   'E0': np.array([0.5, 0.5, 0.5]),
                   'A0': np.array([0.5, 0.0, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def _mon_p_c(self):
        path = [[r'\Gamma', 'Z', 'C', 'Y', r'\Gamma', 'B', 'D', 'E0', 'A0',
                 'Y']]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5]),
                   'C': np.array([0.0, 0.5, 0.5]),
                   'Y': np.array([0.0, 0.5, 0.0]),
                   'B': np.array([0.5, 0.0, 0.0]),
                   'D': np.array([0.5, 0.0, 0.5]),
                   'E0': np.array([0.5, 0.5, 0.5]),
                   'A0': np.array([0.5, 0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def _mon_c_a(self):
        path = [[r'\Gamma', 'Y', 'V', r'\Gamma', 'A', 'M', 'L', 'V']]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Y': np.array([0.5, 0.0, 0.5]),
                   'V': np.array([0.0, 0.0, 0.5]),
                   'A': np.array([0.0, 0.5, 0.0]),
                   'M': np.array([0.5, 0.5, 0.5]),
                   'L': np.array([0.0, 0.5, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def _mon_c_b(self):
        path = [[r'\Gamma', 'Y', 'V', r'\Gamma', 'A', 'M', 'L', 'V']]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Y': np.array([0.5, 0.5, 0.0]),
                   'V': np.array([0.5, 0.0, 0.0]),
                   'A': np.array([0.0, 0.0, 0.5]),
                   'M': np.array([0.5, 0.5, 0.5]),
                   'L': np.array([0.5, 0.0, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def _mon_c_c(self):
        path = [[r'\Gamma', 'Y', 'V', r'\Gamma', 'A', 'M', 'L', 'V']]
        kpoints = {r'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Y': np.array([0.0, 0.5, 0.5]),
                   'V': np.array([0.0, 0.5, 0.0]),
                   'A': np.array([0.5, 0.0, 0.0]),
                   'M': np.array([0.5, 0.5, 0.5]),
                   'L': np.array([0.5, 0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}
