import os
import sys
import subprocess
import re
import math

from subprocess import PIPE

import vaspy

import spglib as spg
import numpy as np

from vaspy.io.vasp_input import Kpoints
from vaspy.misc.math_func import addl, minl, mag, floor_int
from vaspy.misc.math_func import divl, divl_const, multl_const

class SymException(Exception):
    """Symmetry error"""


class Symmetry(object):

    def __init__(self, symmetry, conv_cell):
        self._symmetry = symmetry
        self._conv_cell = conv_cell

    @staticmethod
    def from_poscar(poscar, tolerance=3):
        prec = 1.0/(10**tolerance)

        numbers = []
        for i, n in  enumerate(poscar.natoms):
            numbers.extend([i] * n)
        sym_data = spg.get_symmetry_dataset((poscar.lattice, poscar.coords,
                                             numbers), symprec=prec)
        num = sym_data['number']
        symbol = sym_data['international']

        if 1 <= num <= 2:
            crystal_system = 'triclinic'
        elif 3 <= num <= 15:
            crystal_system = 'monoclinic'
        elif 16 <= num <= 74:
            crystal_system = 'orthorhombic'
        elif 75 <= num <= 142:
            crystal_system = 'tetragonal'
        elif 143 <= num <= 167:
            crystal_system = 'trigonal'
        elif 168 <= num <= 194:
            crystal_system = 'hexagonal'
        elif 195 <= num <= 230:
            crystal_system = 'cubic'

        sym = {'spacegroup_number': num, 'spacegroup_symbol': symbol,
               'crystal_system': crystal_system}

        def abs_cap(val, max_abs_val=1):
            return max(min(val, max_abs_val), -max_abs_val)

        m = sym_data['std_lattice']
        lengths = np.sqrt(np.sum(m ** 2, axis=1))
        angles = np.zeros(3)
        for i in range(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            angles[i] = abs_cap(np.dot(m[j], m[k]) / (lengths[j] * lengths[k]))
        angles = np.arccos(angles) * 180. / np.pi

        cell = np.append(lengths, angles).tolist()
        return Symmetry(sym, cell)

    @property
    def symbol(self):
        return self._symmetry['spacegroup_symbol']

    @property
    def centre(self):
        return self._symmetry['spacegroup_symbol'][0]

    @property
    def number(self):
        return self._symmetry['spacegroup_number']

    @property
    def crystal_system(self):
        return self._symmetry['crystal_system']

    @property
    def conv_cell(self):
        return self._conv_cell  # in the form a,b,c,alpha,beta,gamma

    @property
    def unique_axis(self):
        angles = self.conv_cell[3:]
        return angles.index(min(angles, key=angles.count))

    def str_conv_cell(self):
        s = ", ".join(map(str, self.conv_cell))
        return "a, b, c, alpha, beta, gamma:\n%s" % s

    def __str__(self):
        return 'space group %d: %s centred %s lattice' %\
            (self.number, self.centre, self.crystal_system)

    def __repr__(self):
        return 'space group %d: %s centred %s lattice' %\
            (self.number, self.centre, self.crystal_system)


class PathGen:

    def __init__(self, poscar, symmetry, density=400, old=False):
        self.old_orth_p = old
        self._labels_for_paths, self._special_kpoints_for_paths =\
            self._gen_band_path(symmetry)
        self._kpoints_for_paths = []
        self._nkpoints_per_step_for_paths = []
        for labels, special_kpoints in zip(self._labels_for_paths,
                                           self._special_kpoints_for_paths):
            kpoints, nkpoints_per_step = self._gen_kpoints(poscar.lattice_abc,
                                                           labels,
                                                           special_kpoints,
                                                           density)
            self._kpoints_for_paths.append(kpoints)
            self._nkpoints_per_step_for_paths.append(nkpoints_per_step)

        self._nspecial_kpoints = map(sum, [map(len, self._special_kpoints_for_paths)])
        self._nkpoints = map(sum, [map(len, self._kpoints_for_paths)])

    def _gen_band_path(self, symmetry):
        a, b, c, alpha, beta, gamma = symmetry.conv_cell
        unique = symmetry.unique_axis
        system = symmetry.crystal_system
        centre = symmetry.centre
        self._message = ""

        if system == "triclinic":
            return self.triclinic()

        elif system == 'monoclinic':
            if centre in 'pP':
                if unique == 0:
                    self._message = "a is unique axis"
                    return self.mon_p_a()
                elif unique == 1:
                    self._message = "b is unique axis"
                    return self.mon_p_b()
                elif unique == 2:
                    self._message = "c is unique axis"
                    return self.mon_p_c()
                else:
                    raise SymException("Could not find unique axis for %s" %
                                       symmetry)
            elif centre in 'abcABC':
                if unique == 0:
                    self._message = "a is unique axis"
                    return self.mon_c_a()
                elif unique == 1:
                    self._message = "b is unique axis"
                    return self.mon_c_b()
                elif unique == 2:
                    self._message = "c is unique axis"
                    return self.mon_c_c()
                else:
                    raise SymException("Could not find unique axis for %s" %
                                       symmetry)

        elif system == 'orthorhombic':
            if centre in "pP":
                return self.orth_p()
            elif centre in "abcABC":
                if a > b:
                    self._message = "a > b"
                    return self.orth_c_a()
                elif b > a:
                    self._message = "a > b"
                    return self.orth_c_b()
            elif centre in "fF":
                if (1/a**2 < 1/b**2 + 1/c**2 and 1/b**2 < 1/c**2 + 1/a**2 and
                        1/c**2 < 1/a**2 + 1/b**2):
                    self._message = "1/a^2 < 1/b^2 + 1/c^2 etc..."
                    return self.orth_f_1()
                elif 1/c**2 > 1/a**2 + 1/b**2:
                    self.message = "1/c^2 > 1/a^2 + 1/b^2"
                    return self.orth_f_2()
                elif 1/b**2 > 1/a**2 + 1/c**2:
                    self.message = "1/b^2 > 1/a^2 + 1/c^2"
                    return self.orth_f_3()
                elif 1/a**2 > 1/c**2 + 1/b**2:
                    self.message = "1/a^2 > 1/c^2 + 1/b^2"
                    return self.orth_f_4()
            elif centre in "iI":
                if a > b and a > c:
                    self._message = "a is longest axis"
                    return self.orth_i_a()
                elif b > a and b > c:
                    self._message = "b is longest axis"
                    return self.orth_i_b()
                elif c > a and c > b:
                    self._message = "c is longest axis"
                    return self.orth_i_c()

        elif system == "tetragonal":
            if centre in 'pP':
                return self.tet_p()
            elif centre in 'iI':
                if a > c:
                    self._message = "a > c"
                    return self.tet_i_a()
                else:
                    self._message = "a < c"
                    return self.tet_i_c()

        elif system == "trigonal" or system == "hexagonal":
            if centre in 'rR':
                if a > math.sqrt(2)*c:
                    self._message = "a > sqrt(2)*c"
                    return self.trig_r_a()
                else:
                    self._message = "a < sqrt(2)*c"
                    return self.trig_r_c()
            elif centre in 'pP':
                if unique == 0:
                    self._message = "a is unique axis"
                    return self.trig_p_a()
                elif unique == 2:
                    self._message = "c is unique axis"
                    return self.trig_p_c()
                else:
                    raise SymException("Could not find unique axis for %s" %
                                       symmetry)

        elif system == "cubic":
            if centre in 'pP':
                return self.cubic_p()
            elif centre in 'iI':
                return self.cubic_i()
            elif centre in 'fF':
                return self.cubic_f()
        raise SymException('could not find path for %s' % symmetry)

    def _gen_kpoints(self, lattice_abc, labels, special_kpoints, density):
        # Uses the fact zip generates a list of the smaller of the two input
        # lists to make a new list of lists [(k1, k2), (k2, k3)...]
        # which we can then minus to get the path through kspace
        offset_list = list(special_kpoints)
        offset_list.pop(0)
        steps = [minl(y, x) for x, y in zip(special_kpoints, offset_list)]

        step_lengths = map(mag, [divl(step, lattice_abc) for step in steps])
        nkpoints_per_step = map(floor_int, multl_const(density, step_lengths))
        incs = [divl_const(nkpoints, step)
                for nkpoints, step in zip(nkpoints_per_step, steps)]

        kpoints = []
        for special_kpoint, nkpoints, inc in zip(special_kpoints,
                                                 nkpoints_per_step, incs):
            for i in range(nkpoints):
                if i == (0):
                    # remove any rounding errors that might have cropped up
                    kpoints.append(special_kpoint)
                else:
                    kpoints.append(addl(kpoints[-1], inc))
        # Need this due to the above rounding error avoidance
        kpoints.append(special_kpoints[-1])

        return kpoints, nkpoints_per_step

    @property
    def labels(self):
        return self._labels_for_paths

    @property
    def special_kpoints(self):
        return self._special_kpoints_for_paths

    @property
    def nkpoints_per_step(self):
        return self._nkpoints_per_step_for_paths

    @property
    def nspecial_kpoints(self):
        return self._nspecial_kpoints

    @property
    def nkpoints(self):
        return self._nkpoints

    @property
    def message(self):
        return self._message

    def str_labels(self):
        labels_for_paths = []
        for labels in self._labels_for_paths:
            labels_for_paths.append(" -> ".join(labels))
        return " // ".join(labels_for_paths).rstrip()  # remove trailing newline

    def str_special_kpoints(self):
        special_kpoints_for_paths = []
        for special_kpoints in self._special_kpoints_for_paths:
            s = ""
            for special_kpoint in special_kpoints:
                s += "\t"
                s += "  ".join(str(float(x)) for x in special_kpoint) + "\n"
            special_kpoints_for_paths.append(s)
        return "//\n".join(special_kpoints_for_paths).rstrip()

    def str_nkpoints_per_step(self):
        nkpoints_per_step_for_paths = []
        for nkpoints_per_step in self._nkpoints_per_step_for_paths:
            s = ""
            for step in range(len(nkpoints_per_step)):
                s += "point  %d  to  %d:   number of steps\t%d\n" %\
                    (step+1, step+2, nkpoints_per_step[step])
            nkpoints_per_step_for_paths.append(s)
        return "//\n".join(nkpoints_per_step_for_paths).rstrip()

    def str_nspecial_kpoints(self):
        return " and ".join(str(n) for n in self.nspecial_kpoints)

    def str_nkpoints(self):
        return " and ".join(str(n) for n in self.nkpoints)

    # Should do this like they do in pymatgen, namely:
    # Will do this later
    # elif lattice_type == "triclinic":
    #     kalpha = self._prim_rec.lengths_and_angles[1][0]
    #     kbeta = self._prim_rec.lengths_and_angles[1][1]
    #     kgamma = self._prim_rec.lengths_and_angles[1][2]
    #     if kalpha > 90 and kbeta > 90 and kgamma > 90:
    #         self._kpath = self.tria()
    #     if kalpha < 90 and kbeta < 90 and kgamma < 90:
    #         self._kpath = self.trib()
    #     if kalpha > 90 and kbeta > 90 and kgamma == 90:
    #         self._kpath = self.tria()
    #     if kalpha < 90 and kbeta < 90 and kgamma == 90:
    #         self._kpath = self.trib()

    def triclinic(self):
        labels = [["GP", "Z", "T", "Y", "GP", "X", "V", "R", "U"]]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.0, 0.0, 0.5],
                            [0.0, 0.5, 0.5], [0.0, 0.5, 0.0],
                            [0.0, 0.0, 0.0], [0.5, 0.0, 0.0],
                            [0.5, 0.5, 0.0], [0.5, 0.5, 0.5],
                            [0.5, 0.0, 0.5]]]
        return labels, special_kpoints

    def tet_p(self):
        labels = [["GP", "X", "M", "GP", "Z", "R", "A", "Z"]]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.0, 0.5, 0.0],
                            [0.5, 0.5, 0.0], [0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.5], [0.0, 0.5, 0.5],
                            [0.5, 0.5, 0.5], [0.0, 0.0, 0.5]]]
        return labels, special_kpoints

    def tet_i_a(self):
        labels = [["GP", "X", "P", "N", "GP", "Z"]]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.0, 0.0, 0.5],
                            [0.25, 0.25, 0.25], [0.0, 0.5, 0.0],
                            [0.0, 0.0, 0.0], [-0.5, 0.5, 0.5]]]
        return labels, special_kpoints

    def tet_i_c(self):
        labels = [["GP", "X", "P", "N", "GP", "Z"]]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.0, 0.0, 0.5],
                            [0.25, 0.25, 0.25], [0.0, 0.5, 0.0],
                            [0.0, 0.0, 0.0], [0.5, 0.5, -0.5]]]
        return labels, special_kpoints

    def cubic_p(self):
        # labels = [["GP", "M", "X", "GP", "R", "X"], ["M", "R"]]
        # special_kpoints = [[[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0],
        #                    [0, 0, 0], [0.5, 0.5, 0.5], [0, 0.5, 0]],
        #                    [[0.5, 0.5, 0],[0.5, 0.5, 0.5]]]
        labels = [["GP", "M", "R", "X", "GP", "R"]]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.5, 0.5, 0.0],
                            [0.5, 0.5, 0.5], [0.0, 0.5, 0.0],
                            [0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]]
        return labels, special_kpoints

    def cubic_f(self):
        labels = [["GP", "L", "W", "X", "GP"]]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5],
                            [0.5, 0.25, 0.75], [0.5, 0.0, 0.5],
                            [0.0, 0.0, 0.0]]]
        return labels, special_kpoints

    def cubic_i(self):
        labels = [["GP", "P", "N", "GP", "H", "P"],
                  ["H", "N"]]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.25, 0.25, 0.25],
                            [0.0, 0.0, 0.5], [0.0, 0.0, 0.0],
                            [0.5, -0.5, 0.5], [0.25, 0.25, 0.25]],
                           [[0.5, -0.5, 0.5], [0.0, 0.0, 0.5]]]
        return labels, special_kpoints

    def trig_r_a(self):
        labels = [['GP', 'L', 'F', 'GP', 'Z']]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.0, 0.5, 0.0],
                            [0.5, 0.5, 0.0], [0.0, 0.0, 0.0],
                            [0.5, 0.5, -0.5]]]
        return labels, special_kpoints

    def trig_r_c(self):
        labels = [['GP', 'L', 'F', 'GP', 'Z']]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.0, 0.5, 0.0],
                            [0.5, 0.5, 0.0], [0.0, 0.0, 0.0],
                            [0.5, 0.5, 0.5]]]
        return labels, special_kpoints

    def trig_p_a(self):
        labels = [['GP', 'A', 'L', 'M', 'GP', 'K', 'H', 'A']]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0],
                            [0.5, 0.5, 0.0], [0.0, 0.5, 0.0],
                            [0.0, 0.0, 0.0], [0.0, 0.333, 0.333],
                            [0.5, 0.333, 0.333], [0.5, 0.0, 0.0]]]
        return labels, special_kpoints

#     def trig_p_b(self):
#         labels = [['GP', 'A', 'L', 'M', 'GP', 'K', 'H', 'A']]
#         special_kpoints =
#         return labels, special_kpoints

    def trig_p_c(self):
        labels = [['GP', 'A', 'L', 'M', 'GP', 'K', 'H', 'A']]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.0, 0.0, 0.5],
                            [0.0, 0.5, 0.5], [0.0, 0.5, 0.0],
                            [0.0, 0.0, 0.0], [-0.333, 0.667, 0.0],
                            [-0.333, 0.667, 0.5], [0.0, 0.0, 0.5]]]
        return labels, special_kpoints

    def orth_p(self):
        if self.old_orth_p:
            labels = [['GP', 'Z', 'U', 'X', 'S', 'R', 'T', 'Y', 'GP']]
            special_kpoints = [[[0.0, 0.0, 0.0], [0.0, 0.0, 0.5],
                                [0.0, 0.5, 0.5], [0.0, 0.5, 0.0],
                                [-0.5, 0.5, 0.0], [-0.5, 0.5, 0.5],
                                [-0.5, 0.0, 0.5], [-0.5, 0.0, 0.0],
                                [0.0, 0.0, 0.0]]]
        else:
            labels = [['GP', 'Z', 'T', 'Y', 'S', 'R', 'U', 'X', 'GP', 'Y']]
            special_kpoints = [[[0.0, 0.0, 0.0], [0.0, 0.0, 0.5],
                                [-0.5, 0.0, 0.5], [-0.5, 0.0, 0.0],
                                [-0.5, 0.5, 0.0], [-0.5, 0.5, 0.5],
                                [0.0, 0.5, 0.5], [0.0, 0.5, 0.0],
                                [0.0, 0.0, 0.0], [-0.5, 0.0, 0.0]]]
        return labels, special_kpoints

    def orth_c_a(self):
        labels = [['R', 'S', 'GP', 'Z', 'T', 'Y', 'GP']]
        special_kpoints = [[[0.0, 0.5, 0.5], [0.0, 0.5, 0.0],
                            [0.0, 0.0, 0.0], [0.0, 0.0, 0.5],
                            [0.5, 0.5, 0.5], [0.5, 0.5, 0.0],
                            [0.0, 0.0, 0.0]]]
        return labels, special_kpoints

    def orth_c_b(self):
        labels = [['R', 'S', 'GP', 'Z', 'T', 'Y', 'GP']]
        special_kpoints = [[[0.0, 0.5, 0.5], [0.0, 0.5, 0.0],
                            [0.0, 0.0, 0.0], [0.0, 0.0, 0.5],
                            [-0.5, 0.5, 0.5], [-0.5, 0.5, 0.0],
                            [0.0, 0.0, 0.0]]]
        return labels, special_kpoints

    def orth_f_1(self):
        labels = [['GP', 'Y', 'X', 'Z', 'GP', 'L']]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.0, -0.5, -0.5],
                            [0.5, 0.0, 0.5], [0.5, 0.5, 0.0],
                            [0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]]
        return labels, special_kpoints

    def orth_f_2(self):
        labels = [['GP', 'Y', 'X', 'Z', 'GP', 'L']]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.0, -0.5, -0.5],
                            [0.5, 0.0, 0.5], [0.5, -0.5, 0.0],
                            [0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]]
        return labels, special_kpoints

    def orth_f_3(self):
        labels = [['GP', 'Y', 'X', 'Z', 'GP', 'L']]
        special_kpoints = [[[0.0, 0.0, 0.0], [1.0, 0.5, 0.5],
                            [0.5, 0.0, 0.5], [0.5, 0.5, 0.0],
                            [0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]]
        return labels, special_kpoints

    def orth_f_4(self):
        labels = [['GP', 'Y', 'X', 'Z', 'GP', 'L']]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.0, -0.5, -0.5],
                            [0.5, 0.0, -0.5], [0.5, 0.5, 0.0],
                            [0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]]
        return labels, special_kpoints

    def orth_i_a(self):
        labels = [['R', 'GP', 'X', 'S', 'W', 'T']]
        special_kpoints = [[[0.5, 0.0, 0.0], [0.0, 0.0, 0.0],
                            [0.5, -0.5, 0.5], [0.5, 0.0, -0.5],
                            [0.75, -0.25, -0.25], [0.5, -0.5, 0.0]]]
        return labels, special_kpoints

    def orth_i_b(self):
        labels = [['R', 'GP', 'X', 'S', 'W', 'T']]
        special_kpoints = [[[0.5, 0.0, 0.0], [0.0, 0.0, 0.0],
                            [0.5, -0.5, -0.5], [0.5, 0.0, 0.5],
                            [0.75, -0.25, -0.25], [0.5, -0.5, 0.0]]]
        return labels, special_kpoints

    def orth_i_c(self):
        labels = [['R', 'GP', 'X', 'S', 'W', 'T']]
        special_kpoints = [[[0.5, 0.0, 0.0], [0.0, 0.0, 0.0],
                            [0.5, 0.5, -0.5], [0.5, 0.0, -0.5],
                            [0.75, -0.25, -0.25], [0.5, -0.5, 0.0]]]
        return labels, special_kpoints

    def mon_p_a(self):
        labels = [['GP', 'Z', 'C', 'Y', 'GP', 'B', 'D', 'E0', 'A0', 'Y']]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0],
                            [0.5, 0.0, 0.5], [0.0, 0.0, 0.5],
                            [0.0, 0.0, 0.0], [0.0, 0.5, 0.0],
                            [0.5, 0.5, 0.0], [0.5, 0.5, 0.5],
                            [0.0, 0.5, 0.5], [0.0, 0.0, 0.5]]]
        return labels, special_kpoints

    def mon_p_b(self):
        labels = [['GP', 'Z', 'C', 'Y', 'GP', 'B', 'D', 'E0', 'A0', 'Y']]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.0, 0.5, 0.0],
                            [0.5, 0.5, 0.0], [0.5, 0.0, 0.0],
                            [0.0, 0.0, 0.0], [0.0, 0.0, 0.5],
                            [0.0, 0.5, 0.5], [0.5, 0.5, 0.5],
                            [0.5, 0.0, 0.5], [0.5, 0.0, 0.0]]]
        return labels, special_kpoints

    def mon_p_c(self):
        labels = [['GP', 'Z', 'C', 'Y', 'GP', 'B', 'D', 'E0', 'A0', 'Y']]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.0, 0.0, 0.5],
                            [0.0, 0.5, 0.5], [0.0, 0.5, 0.0],
                            [0.0, 0.0, 0.0], [0.5, 0.0, 0.0],
                            [0.5, 0.0, 0.5], [0.5, 0.5, 0.5],
                            [0.5, 0.5, 0.0], [0.0, 0.5, 0.0]]]
        return labels, special_kpoints

    def mon_c_a(self):
        labels = [['GP', 'Y', 'V', 'GP', 'A', 'M', 'L', 'V']]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.5, 0.0, 0.5],
                            [0.0, 0.0, 0.5], [0.0, 0.0, 0.0],
                            [0.0, 0.5, 0.0], [0.5, 0.5, 0.5],
                            [0.0, 0.5, 0.5], [0.0, 0.0, 0.5]]]
        return labels, special_kpoints

    def mon_c_b(self):
        labels = [['GP', 'Y', 'V', 'GP', 'A', 'M', 'L', 'V']]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.5, 0.5, 0.0],
                            [0.5, 0.0, 0.0], [0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.5], [0.5, 0.5, 0.5],
                            [0.5, 0.0, 0.5], [0.5, 0.0, 0.0]]]
        return labels, special_kpoints

    def mon_c_c(self):
        labels = [['GP', 'Y', 'V', 'GP', 'A', 'M', 'L', 'V']]
        special_kpoints = [[[0.0, 0.0, 0.0], [0.0, 0.5, 0.5],
                            [0.0, 0.5, 0.0], [0.0, 0.0, 0.0],
                            [0.5, 0.0, 0.0], [0.5, 0.5, 0.5],
                            [0.5, 0.5, 0.0], [0.0, 0.5, 0.0]]]
        return labels, special_kpoints

    # Returns a list of filenames incase we need to make folders with them
    # hybrid = 0 for dft, 1 for hybrid, 2 for both
    def to_files(self, hybrid=0, split=0, ibzkpt=None):
        filenames = []
        id = 1
        for kpoints in self._kpoints_for_paths:
            kpoints_file = Kpoints(kpoints)

            if hybrid == 0 or hybrid == 2:
                filename = "KPOINTS_%d" % id
                filenames += kpoints_file.to_file(filename)

            if hybrid == 1 or hybrid == 2:
                if split != 0:
                    split_kpoints = kpoints_file.split(split)
                    split_id = 1
                    for split_kpoint in split_kpoints:
                        if ibzkpt is not None:
                            split_kpoint.prepend(ibzkpt)
                        filename = "KPOINTS_%d_SPLIT_%d" % (id, split_id)
                        filenames += split_kpoint.to_file(filename,
                                                          pbe=False, hse=True)
                        split_id += 1
                else:
                    if ibzkpt is not None:
                        kpoints_file.prepend(ibzkpt)
                    filename = "KPOINTS_%d" % id
                    filenames += kpoints_file.to_file(filename,
                                                      pbe=False, hse=True)
            id += 1

        return filenames
