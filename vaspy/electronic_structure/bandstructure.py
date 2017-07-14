import math
import spglib
import seekpath

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.groups import SpaceGroup
from pymatgen.symmetry.bandstructure import HighSymmKpath


class Kpath(object):
    """
    Dummy class providing helper functions
    """

    def __init__(self, structure, symprec=1e-3):
        self.structure = structure

        # use sym as a quick way to access the cell data
        sym = SpacegroupAnalyzer(structure, symprec=symprec)
        self._spg_data = spglib.get_symmetry_dataset(sym._cell)

        # make primitive and conventional cell from seekpath output
        std = spglib.standardize_cell(sym._cell, symprec=symprec,
                                      to_primitive=True)
        self._seek_data = seekpath.get_path(std)

        prim_lattice = self._seek_data['primitive_lattice']
        prim_scaled_positions = self._seek_data['primitive_positions']
        prim_numbers = self._seek_data['primitive_types']
        prim_atoms = [sym._unique_species[i - 1] for i in prim_numbers]
        self.prim = Structure(prim_lattice, prim_atoms, prim_scaled_positions)

        conv_lattice = self._seek_data['conv_lattice']
        conv_scaled_positions = self._seek_data['conv_positions']
        conv_numbers = self._seek_data['conv_types']
        conv_atoms = [sym._unique_species[i - 1] for i in conv_numbers]
        self.conv = Structure(conv_lattice, conv_atoms, conv_scaled_positions)

    @property
    def kpoints(self):
        return self._kpath['kpoints']

    @property
    def path(self):
        return self._kpath['path']

    @property
    def lattice_type(self):
        return get_lattice_type(self._spg_data['number'])

    @property
    def spg_symbol(self):
        return self._spg_data['international']

    @property
    def path_string(self):
        path_str = ' | '.join([' -> '.join(subpath) for subpath in self.path])
        return path_str

    def correct_structure(self, atol=1e-8):
        return np.allclose(self.structure.lattice.matrix,
                           self.prim.lattice.matrix, atol=atol)

    def get_kpoints(self, line_density=20, cart_coords=False):
        return get_kpoints(self.structure, self.kpoints, self.path,
                           line_density=line_density, cart_coords=cart_coords)

class PymatgenKpath(Kpath):

    def __init__(self, structure, symprec=1e-3):
        Kpath.__init__(self, structure, symprec=symprec)
        pmg_path = HighSymmKpath(structure, symprec=symprec)
        self._kpath = pmg_path._kpath
        self.prim = pmg_path.prim
        self.conv = pmg_path.conventional


class SeekpathKpath(Kpath):

    def __init__(self, structure, symprec=1e-3):
        Kpath.__init__(self, structure, symprec=symprec)

        # need to convert from seekpath format to something useable
        path = [[self._seek_data['path'][0][0]]]
        for (k1, k2) in self._seek_data['path']:
            if path[-1] and path[-1][-1] == k1:
                path[-1].append(k2)
            else:
                path.append([k1, k2])

        # change gamma label to \Gamma
        kpoints = self._seek_data['point_coords']
        kpoints['\Gamma'] = kpoints.pop('GAMMA')
        path = [[label.replace('GAMMA', '\Gamma') for label in subpath]
                for subpath in path]

        # remove unused k-points
        # TODO: this but better
        pts = []
        for subpath in path:
            pts += subpath
        pts = list(set(pts))
        pts_coords = [kpoints[p] for p in pts]
        kpoints = dict(zip(pts, pts_coords))
        self._kpath = {'kpoints': kpoints, 'path': path}


class BradCrackKpath(Kpath):

    def __init__(self, structure, symprec=1e-3, spg=None):
        Kpath.__init__(self, structure, symprec=symprec)

        angles = self.conv.lattice.angles
        unique = angles.index(min(angles, key=angles.count))
        a = self.conv.lattice.abc[0]
        b = self.conv.lattice.abc[1]
        c = self.conv.lattice.abc[2]

        if spg:
            spg_symbol = spg.symbol
            lattice_type = get_lattice_type(spg.int_number)
        else:
            spg_symbol = self.spg_symbol
            lattice_type = self.lattice_type

        if lattice_type == 'triclinic':
            self._kpath == self.triclinic()

        elif lattice_type == 'monoclinic':
            if 'P' in spg_symbol:
                if unique == 0:
                    self._kpath = self.mon_p_a()
                elif unique == 1:
                    self._kpath = self.mon_p_b()
                elif unique == 2:
                    self._kpath = self.mon_p_c()

            elif 'C' in spg_symbol:
                if unique == 0:
                    return self.mon_c_a()
                elif unique == 1:
                    return self.mon_c_b()
                elif unique == 2:
                    return self.mon_c_c()

        elif lattice_type == 'orthorhombic':
            if 'P' in spg_symbol:
                self._kpath = self.orth_p()

            elif 'C' in spg_symbol:
                if a > b:
                    self._kpath = self.orth_c_a()
                elif b > a:
                    self._kpath = self.orth_c_b()

            elif 'F' in spg_symbol:
                if (1/a**2 < 1/b**2 + 1/c**2 and 1/b**2 < 1/c**2 + 1/a**2 and
                        1/c**2 < 1/a**2 + 1/b**2):
                    self._kpath = self.orth_f_1()
                elif 1/c**2 > 1/a**2 + 1/b**2:
                    self._kpath = self.orth_f_2()
                elif 1/b**2 > 1/a**2 + 1/c**2:
                    self._kpath = self.orth_f_3()
                elif 1/a**2 > 1/c**2 + 1/b**2:
                    self._kpath = self.orth_f_4()

            elif 'I' in spg_symbol:
                if a > b and a > c:
                    self._kpath = self.orth_i_a()
                elif b > a and b > c:
                    self._kpath = self.orth_i_b()
                elif c > a and c > b:
                    self._kpath = self.orth_i_c()

        elif lattice_type == 'tetragonal':
            if 'P' in spg_symbol:
                self._kpath = self.tet_p()

            elif 'I' in spg_symbol:
                if a > c:
                    self._kpath = self.tet_i_a()
                else:
                    self._kpath = self.tet_i_c()

        elif lattice_type == 'trigonal' or lattice_type == 'hexagonal':
            if 'R' in spg_symbol:
                if a > math.sqrt(2) * c:
                    self._kpath = self.trig_r_a()
                else:
                    self._kpath = self.trig_r_c()

            elif 'P' in spg_symbol:
                if unique == 0:
                    self._kpath = self.trig_p_a()
                elif unique == 2:
                    self._kpath = self.trig_p_c()

        elif lattice_type == "cubic":
            if 'P' in spg_symbol:
                self._kpath = self.cubic_p()
            elif 'I' in spg_symbol:
                self._kpath = self.cubic_i()
            elif 'F' in spg_symbol:
                self._kpath = self.cubic_f()

    def triclinic(self):
        path = [["\Gamma", "Z", "T", "Y", "\Gamma", "X", "V", "R", "U"]]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5]),
                   'T': np.array([0.0, 0.5, 0.5]),
                   'Y': np.array([0.0, 0.5, 0.0]),
                   'X': np.array([0.5, 0.0, 0.0]),
                   'V': np.array([0.5, 0.5, 0.0]),
                   'R': np.array([0.5, 0.5, 0.5]),
                   'U': np.array([0.5, 0.0, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def tet_p(self):
        path = [["\Gamma", "X", "M", "\Gamma", "Z", "R", "A", "Z"]]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'X': np.array([0.0, 0.5, 0.0]),
                   'M': np.array([0.5, 0.5, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5]),
                   'R': np.array([0.0, 0.5, 0.5]),
                   'A': np.array([0.5, 0.5, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def tet_i_a(self):
        path = [["\Gamma", "X", "P", "N", "\Gamma", "Z"]]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'X': np.array([0.0, 0.0, 0.5]),
                   'P': np.array([0.25, 0.25, 0.25]),
                   'N': np.array([0.0, 0.5, 0.0]),
                   'Z': np.array([-0.5, 0.5, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def tet_i_c(self):
        path = [["\Gamma", "X", "P", "N", "\Gamma", "Z"]]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'X': np.array([0.0, 0.0, 0.5]),
                   'P': np.array([0.25, 0.25, 0.25]),
                   'N': np.array([0.0, 0.5, 0.0]),
                   'Z': np.array([0.5, 0.5, -0.5])}
        return {'kpoints': kpoints, 'path': path}

    def cubic_p(self):
        path = [["\Gamma", "M", "R", "X", "\Gamma", "R"]]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'M': np.array([0.5, 0.5, 0.0]),
                   'R': np.array([0.5, 0.5, 0.5]),
                   'X': np.array([0.0, 0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def cubic_f(self):
        path = [["\Gamma", "L", "W", "X", "\Gamma"]]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'W': np.array([0.5, 0.25, 0.75]),
                   'X': np.array([0.5, 0.0, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def cubic_i(self):
        path = [["\Gamma", "P", "N", "\Gamma", "H", "P"], ["H", "N"]]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'P': np.array([0.25, 0.25, 0.25]),
                   'N': np.array([0.0, 0.0, 0.5]),
                   'H': np.array([0.5, -0.5, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def trig_r_a(self):
        path = [['\Gamma', 'L', 'F', '\Gamma', 'Z']]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'L': np.array([0.0, 0.5, 0.0]),
                   'F': np.array([0.5, 0.5, 0.0]),
                   'Z': np.array([0.5, 0.5, -0.5])}
        return {'kpoints': kpoints, 'path': path}

    def trig_r_c(self):
        path = [['\Gamma', 'L', 'F', '\Gamma', 'Z']]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'L': np.array([0.0, 0.5, 0.0]),
                   'F': np.array([0.5, 0.5, 0.0]),
                   'Z': np.array([0.5, 0.5, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def trig_p_a(self):
        path = [['\Gamma', 'A', 'L', 'M', '\Gamma', 'K', 'H', 'A']]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([0.5, 0.0, 0.0]),
                   'L': np.array([0.5, 0.5, 0.0]),
                   'M': np.array([0.0, 0.5, 0.0]),
                   'K': np.array([0.0, 0.333, 0.333]),
                   'H': np.array([0.5, 0.333, 0.333])}
        return {'kpoints': kpoints, 'path': path}

    def trig_p_c(self):
        path = [['\Gamma', 'A', 'L', 'M', '\Gamma', 'K', 'H', 'A']]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([0.0, 0.0, 0.5]),
                   'L': np.array([0.0, 0.5, 0.5]),
                   'M': np.array([0.0, 0.5, 0.0]),
                   'K': np.array([-0.333, 0.667, 0.0]),
                   'H': np.array([-0.333, 0.667, 0.6])}
        return {'kpoints': kpoints, 'path': path}

    def orth_p(self):
        path = [['\Gamma', 'Z', 'T', 'Y', 'S', 'R', 'U', 'X', '\Gamma', 'Y']]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5]),
                   'T': np.array([-0.5, 0.0, 0.5]),
                   'Y': np.array([-0.5, 0.0, 0.0]),
                   'S': np.array([-0.5, 0.5, 0.0]),
                   'R': np.array([-0.5, 0.5, 0.5]),
                   'U': np.array([0.0, 0.5, 0.5]),
                   'X': np.array([0.0, 0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def orth_c_a(self):
        path = [['R', 'S', '\Gamma', 'Z', 'T', 'Y', '\Gamma']]
        kpoints = {'R': np.array([0.0, 0.5, 0.5]),
                   'S': np.array([0.0, 0.5, 0.0]),
                   '\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5]),
                   'T': np.array([0.5, 0.5, 0.5]),
                   'Y': np.array([0.5, 0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def orth_c_b(self):
        path = [['R', 'S', '\Gamma', 'Z', 'T', 'Y', '\Gamma']]
        kpoints = {'R': np.array([0.0, 0.5, 0.5]),
                   'S': np.array([0.0, 0.5, 0.0]),
                   '\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5]),
                   'T': np.array([-0.5, 0.5, 0.5]),
                   'Y': np.array([-0.5, 0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def orth_f_1(self):
        path = [['\Gamma', 'Y', 'X', 'Z', '\Gamma', 'L']]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Y': np.array([0.0, -0.5, -0.5]),
                   'X': np.array([0.5, 0.0, 0.5]),
                   'Z': np.array([0.5, 0.5, 0.0]),
                   'L': np.array([0.5, 0.0, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def orth_f_2(self):
        path = [['\Gamma', 'Y', 'X', 'Z', '\Gamma', 'L']]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Y': np.array([0.0, -0.5, -0.5]),
                   'X': np.array([0.5, 0.0, 0.5]),
                   'Z': np.array([0.5, -0.5, 0.0]),
                   'L': np.array([0.5, 0.0, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def orth_f_3(self):
        path = [['\Gamma', 'Y', 'X', 'Z', '\Gamma', 'L']]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Y': np.array([1.0, 0.5, 0.5]),
                   'X': np.array([0.5, 0.0, 0.5]),
                   'Z': np.array([0.5, 0.5, 0.0]),
                   'L': np.array([0.5, 0.0, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def orth_f_4(self):
        path = [['\Gamma', 'Y', 'X', 'Z', '\Gamma', 'L']]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Y': np.array([0.0, -0.5, -0.5]),
                   'X': np.array([0.5, 0.0, -0.5]),
                   'Z': np.array([0.5, 0.5, 0.0]),
                   'L': np.array([0.5, 0.0, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def orth_i_a(self):
        path = [['R', '\Gamma', 'X', 'S', 'W', 'T']]
        kpoints = {'R': np.array([0.5, 0.0, 0.0]),
                   '\Gamma': np.array([0.0, 0.0, 0.0]),
                   'X': np.array([0.5, -0.5, 0.5]),
                   'S': np.array([0.5, 0.0, -0.5]),
                   'W': np.array([0.75, -0.25, -0.25]),
                   'T': np.array([0.5, -0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def orth_i_b(self):
        path = [['R', '\Gamma', 'X', 'S', 'W', 'T']]
        kpoints = {'R': np.array([0.5, 0.0, 0.0]),
                   '\Gamma': np.array([0.0, 0.0, 0.0]),
                   'X': np.array([0.5, -0.5, -0.5]),
                   'S': np.array([0.5, 0.0, -0.5]),
                   'W': np.array([0.75, -0.25, -0.25]),
                   'T': np.array([0.5, -0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def orth_i_c(self):
        path = [['R', '\Gamma', 'X', 'S', 'W', 'T']]
        kpoints = {'R': np.array([0.5, 0.0, 0.0]),
                   '\Gamma': np.array([0.0, 0.0, 0.0]),
                   'X': np.array([0.5, 0.5, -0.5]),
                   'S': np.array([0.5, 0.0, -0.5]),
                   'W': np.array([0.75, -0.25, -0.25]),
                   'T': np.array([0.5, -0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def mon_p_a(self):
        path = [['\Gamma', 'Z', 'C', 'Y', '\Gamma', 'B', 'D', 'E0', 'A0', 'Y']]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Z': np.array([0.5, 0.0, 0.0]),
                   'C': np.array([0.5, 0.0, 0.5]),
                   'Y': np.array([0.0, 0.0, 0.5]),
                   'B': np.array([0.0, 0.5, 0.0]),
                   'D': np.array([0.5, 0.5, 0.0]),
                   'E0': np.array([0.5, 0.5, 0.5]),
                   'A0': np.array([0.0, 0.5, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def mon_p_b(self):
        path = [['\Gamma', 'Z', 'C', 'Y', '\Gamma', 'B', 'D', 'E0', 'A0', 'Y']]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Z': np.array([0.5, 0.0, 0.0]),
                   'C': np.array([0.5, 0.5, 0.0]),
                   'Y': np.array([0.5, 0.0, 0.0]),
                   'B': np.array([0.0, 0.0, 0.5]),
                   'D': np.array([0.0, 0.5, 0.5]),
                   'E0': np.array([0.5, 0.5, 0.5]),
                   'A0': np.array([0.5, 0.0, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def mon_p_c(self):
        path = [['\Gamma', 'Z', 'C', 'Y', '\Gamma', 'B', 'D', 'E0', 'A0', 'Y']]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5]),
                   'C': np.array([0.0, 0.5, 0.5]),
                   'Y': np.array([0.0, 0.5, 0.0]),
                   'B': np.array([0.5, 0.0, 0.0]),
                   'D': np.array([0.5, 0.0, 0.5]),
                   'E0': np.array([0.5, 0.5, 0.5]),
                   'A0': np.array([0.5, 0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}

    def mon_c_a(self):
        path = [['\Gamma', 'Y', 'V', '\Gamma', 'A', 'M', 'L', 'V']]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Y': np.array([0.5, 0.0, 0.5]),
                   'V': np.array([0.0, 0.0, 0.5]),
                   'A': np.array([0.0, 0.5, 0.0]),
                   'M': np.array([0.5, 0.5, 0.5]),
                   'L': np.array([0.0, 0.5, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def mon_c_b(self):
        path = [['\Gamma', 'Y', 'V', '\Gamma', 'A', 'M', 'L', 'V']]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Y': np.array([0.5, 0.5, 0.0]),
                   'V': np.array([0.5, 0.0, 0.0]),
                   'A': np.array([0.0, 0.0, 0.5]),
                   'M': np.array([0.5, 0.5, 0.5]),
                   'L': np.array([0.5, 0.0, 0.5])}
        return {'kpoints': kpoints, 'path': path}

    def mon_c_c(self):
        path = [['\Gamma', 'Y', 'V', '\Gamma', 'A', 'M', 'L', 'V']]
        kpoints = {'\Gamma': np.array([0.0, 0.0, 0.0]),
                   'Y': np.array([0.0, 0.5, 0.5]),
                   'V': np.array([0.0, 0.5, 0.0]),
                   'A': np.array([0.5, 0.0, 0.0]),
                   'M': np.array([0.5, 0.5, 0.5]),
                   'L': np.array([0.5, 0.5, 0.0])}
        return {'kpoints': kpoints, 'path': path}


def get_kpoints(structure, kpoints, path, line_density=20, cart_coords=False):
    list_k_points = []
    sym_point_labels = []
    recip_lattice = structure.lattice.reciprocal_lattice
    for b in path:
        for i in range(1, len(b)):
            start = np.array(kpoints[b[i - 1]])
            end = np.array(kpoints[b[i]])
            distance = np.linalg.norm(
                recip_lattice.get_cartesian_coords(start) -
                recip_lattice.get_cartesian_coords(end))
            nb = int(math.ceil(distance * line_density))
            sym_point_labels.extend([b[i - 1]] + [''] * (nb - 1))
            list_k_points.extend(
                [recip_lattice.get_cartesian_coords(start)
                    + float(i) / float(nb) *
                    (recip_lattice.get_cartesian_coords(end)
                     - recip_lattice.get_cartesian_coords(start))
                    for i in range(0, nb)])
        # append last k-point to avoid repition as in pymatgen
        sym_point_labels.append(b[-1])
        list_k_points.append(recip_lattice.get_cartesian_coords(end))
    if cart_coords:
        return list_k_points, sym_point_labels
    else:
        frac_k_points = [recip_lattice.get_fractional_coords(k)
                         for k in list_k_points]
        return frac_k_points, sym_point_labels


def get_lattice_type(number):
    """
    Adapted from pymatgen.symmetry.analyzer.SpacegroupAnalyzer
    """
    f = lambda i, j: i <= number <= j
    cs = {'triclinic': (1, 2), 'monoclinic': (3, 15),
          'orthorhombic': (16, 74), 'tetragonal': (75, 142),
          'trigonal': (143, 167), 'hexagonal': (168, 194),
          'cubic': (195, 230)}

    crystal_system = None
    for k, v in cs.items():
        if f(*v):
            crystal_system = k
            break

    if number in [146, 148, 155, 160, 161, 166, 167]:
        return "rhombohedral"
    elif crystal_system == "trigonal":
        return "hexagonal"
    else:
        return crystal_system
