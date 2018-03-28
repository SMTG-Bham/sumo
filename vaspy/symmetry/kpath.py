# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Module containing base class for high-symmetry k-point path determination.
"""

import spglib
import seekpath

import numpy as np

from vaspy.symmetry.kpoints import get_kpoints

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class Kpath(object):
    r"""Base class providing helper functions for generating k-point paths.

    This class should not be used directly. Instead, one of the PymatgenKpath,
    SeekpathKpath, or BradCrackKpath subclasses should be used. The main use
    of this parent object is for standardisation across the differing k-point
    path generation classes.

    Args:
        structure (:obj:`~pymatgen.core.structure.Structure`): The structure.
        symprec (:obj:`float`, optional): The tolerance for determining the
            crystal symmetry.

    Attributes:
        prim (:obj:`~pymatgen.core.structure.Structure`): The standardised
            primitive cell structure for the generated k-point path.
        conv (:obj:`~pymatgen.core.structure.Structure`): The standardised
            conventional cell structure.
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
        """Determine if the structure matches the standard primitive structure.

        The standard primitive will be different between seekpath and pymatgen
        high-symmetry paths, but this is handled by the specific subclasses.

        Args:
            atol (:obj:`float`, optional): Absolute tolerance used to compare
                the input structure with the primitive standard structure.

        Returns:
            bool: ``True`` if the structure is the same as the standard
            primtive, otherwise ``False``.
        """
        return np.allclose(self.structure.lattice.matrix,
                           self.prim.lattice.matrix, atol=atol)

    def get_kpoints(self, line_density=20, cart_coords=False, phonopy=False):
        r"""Return a list of k-points and labels along the high-symmetry path.

        Args:
            line_density (:obj:`int`, optional): Density of k-points along the
                path.
            cart_coords (:obj:`bool`, optional): Whether the k-points are
                returned in cartesian or reciprocal coordinates. Defaults to
                ``False`` (fractional coordinates).
            phonopy (:obj:`bool`, optional): Format the k-points and labels for
                use with phonopy. Defaults to ``False``.

        Returns:
            tuple: A :obj:`tuple` of the k-points along the high-symmetry path,
            and k-point labels. Returned as ``(kpoints, labels)``.

            If ``phonopy == False``, then:

                * ``kpoints`` is a :obj:`numpy.ndarray` of the k-point
                  coordinates along the high-symmetry path. For example::

                      [[0, 0, 0], [0.25, 0, 0], [0.5, 0, 0], [0.5, 0, 0.25],
                       [0.5, 0, 0.5]]

                * ``labels`` is a :obj:`list` of the high symmetry labels for
                  each k-point (will be an empty :obj:`str` if the k-point has
                  no label). For example::

                      ['\Gamma', '', 'X', '', 'Y']

            If ``phonopy == True``, then:

                * ``kpoints`` is a :obj:`list` of :obj:`numpy.ndarray`
                  containing the k-points for each branch of the band
                  structure. This means that the first and last k-points of a
                  particular branch may be repeated. For example::

                      [[[0, 0, 0], [0.25, 0, 0], [0.5, 0, 0]],
                       [[0.5, 0, 0], [0.5, 0, 0.25], [0.5, 0, 0.5]]]

                * ``labels`` is a :obj:`list` of the high symmetry labels.
                  For example::

                      ['\Gamma', 'X', 'Y']
        """
        return get_kpoints(self.structure, self.kpoints, self.path,
                           line_density=line_density, cart_coords=cart_coords,
                           phonopy=phonopy)

    @property
    def kpoints(self):
        r"""dict: The high-symmetry k-point labels and their fractional
        coordinates. Formatted as ``{label: coords}``. For example::

            {'\Gamma': [0., 0., 0.], 'X': [0.5, 0. 0.]}
        """
        return self._kpath['kpoints']

    @property
    def path(self):
        r"""list: The high-symmetry k-point path. Each subpath is provided
        as a list. For example, the following covers the path ``\Gamma ->
        X -> C | \Gamma -> Y``::

            [['\Gamma', 'X', 'C'], ['\Gamma', 'Y']].
        """
        return self._kpath['path']

    @property
    def lattice_type(self):
        """str: The Bravais lattice system. Hexagonal cells are separated into
        rhombohedral and hexagonal lattices.
        """
        return self.get_lattice_type(self.spg_number)

    @property
    def spg_symbol(self):
        """str: The international space group symbol."""
        return self._spg_data['international']

    @property
    def spg_number(self):
        """int: The international space group number."""
        return self._spg_data['number']

    @property
    def path_string(self):
        """str: The high-symmetry k-point path formatted with arrows and
        showing disconnections between subpaths. For example::

            "X -> Gamma | Y -> Z"
        """
        return ' | '.join([' -> '.join(subpath) for subpath in self.path])

    @staticmethod
    def get_lattice_type(number):
        """Return the lattice crystal system.

        Hexagonal cells are differentiated into rhombohedral and hexagonal
        lattices.

        Args:
            number (int): The international space group number.

        Returns:
            str: The lattice crystal system.
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
