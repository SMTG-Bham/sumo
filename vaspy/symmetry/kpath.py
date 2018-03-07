# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

import math
import spglib
import seekpath

import numpy as np

from vaspy.symmetry.kpoints import get_kpoints

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class Kpath(object):
    """Dummy class providing helper functions for generating k-point paths.

    This class should not be used directly. Instead, one of the PymatgenKpath,
    SeekpathKpath, or BradCrackKpath subclasses should be used. The main use
    of this parent object is for standardisation across the differing k-point
    path generation tools.

    Args:
        structure (Structure): A pymatgen structure object.
        symprec (float): The tolerance for determining the crystal symmetry.

    Attributes:
        kpoints (dict): The high-symmetry k-point labels and their coordinates
            as {label: coords}.
        path (list): The high-symmetry k-point path. Each subpath is provided
            as a list. E.g. [['A', 'B'], ['C', 'D']].
        prim (Structure): The standardised primitive cell structure needed for
            to obtain the correct band structure.
        conv (Structure): The standardised conventional cell structure.
        lattice_type (str): The Bravais lattice system. Hexagonal cells are
            separated into rhombohedral and hexagonal lattices.
        spg_symbol (str): The international space group symbol.
        spg_number (int): The international space group number.
        path_string (str): The high-symmetry k-point path formatted with arrows
            and showing disconnections between subpaths. For example:
            "X -> Gamma | Y -> Z".
    """

    def __init__(self, structure, symprec=1e-3):
        self.structure = structure

        # use sym as a quick way to access the cell data
        sym = SpacegroupAnalyzer(structure, symprec=symprec)
        self._spg_data = sym.get_symmetry_dataset()

        # make primitive and conventional cell from seekpath output
        std = spglib.refine_cell(sym._cell, symprec=symprec)
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

    def correct_structure(self, atol=1e-8):
        """Determine if the structure matches the standard primitive.

        The standard primitive will be different between seekpath and pymatgen
        high-symmetry paths, but this is handled by the specific subclasses.

        Args:
            atol (float): Absolute tolerance used to compare the input
                structure with the one expected as primitive standard.

        Returns:
            True if the structure is the same as the standard primtive, False
            otherwise.
        """
        return np.allclose(self.structure.lattice.matrix,
                           self.prim.lattice.matrix, atol=atol)

    def get_kpoints(self, line_density=20, cart_coords=False, phonopy=False):
        """Calculate a list of k-points along the high-symmetry path.

as       Args:
            line_density (int): The density of k-points along the path.
            cart_coords (bool): Whether the k-points are returned in cartesian
                or reciprocal coordinates.
            phonopy (bool): Format the k-points and labels for use with phonopy.

        Returns:
            A list k-points along the high-symmetry path, together with the
            high symmetry labels for each k-point. Returned as: kpoints, labels.
        """
        return get_kpoints(self.structure, self.kpoints, self.path,
                           line_density=line_density, cart_coords=cart_coords,
                           phonopy=phonopy)

    @property
    def kpoints(self):
        return self._kpath['kpoints']

    @property
    def path(self):
        return self._kpath['path']

    @property
    def lattice_type(self):
        return self.get_lattice_type(self.spg_number)

    @property
    def spg_symbol(self):
        return self._spg_data['international']

    @property
    def spg_number(self):
        return self._spg_data['number']

    @property
    def path_string(self):
        return ' | '.join([' -> '.join(subpath) for subpath in self.path])

    @staticmethod
    def get_lattice_type(number):
        """Obtain the lattice crystal system.

        Hexagonal cells are differentiated into rhombohedral and hexagonal lattices.
        Adapted from pymatgen.symmetry.analyzer.SpacegroupAnalyzer

        Args:
            number (int): The international space group number.

        Returns:
            The lattice crystal system as a string.
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
