# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

import string
import numpy as np

from vaspy.symmetry import Kpath


class CustomKpath(Kpath):
    """Class to represent custom k-point paths.

    Args:
        structure (Structure): A pymatgen structure object.
        path_labels (list): A list of labels to use along with kpt_list. These
            should be provided as a list of subpaths, each containing a list of
            labels. For example: [['Gamma', 'Z'], ['X', 'M']], combined with
            the above kpt_list would indicate the path: Gamma -> Z | X -> M.
        symprec (float): The tolerance for determining the crystal symmetry.
        kpt_list (list): Manual list of k-points to use. If kpt_list is set it
            will override the mode selection. Should be formatted as a list of
            subpaths, each containing a list of k-points. For example:
            [[[0., 0., 0.], [0., 0., 0.5]], [[0.5, 0., 0.], [0.5, 0.5, 0.]]]

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

    def __init__(self, structure, kpt_list, path_labels, symprec=1e-3):
        Kpath.__init__(self, structure, symprec=symprec)

        # TODO: add warnings for no labels and incorrect number of labels
        flat_kpts = [x for kpts in kpt_list for x in kpts]
        if path_labels:
            flat_path_labels = [x for labels in path_labels for x in labels]
        else:
            flat_path_labels = [s for s, x in
                                zip(string.ascii_uppercase, flat_kpts)]

        # need this to make sure repeated kpoints have the same labels
        kpt_dict = {}
        for label, kpt in zip(flat_path_labels, flat_kpts):
            if kpt not in kpt_dict.values():
                kpt_dict[label] = kpt

        if not path_labels:
            path_labels = []
            for kpt_sublist in kpt_list:
                labels = []
                for kpt in kpt_sublist:
                    for label, kpt2 in iter(kpt_dict.items()):
                        if np.array_equal(kpt, kpt2):
                            labels.append(label)
                            break
                path_labels.append(labels)

        self._kpath = {'kpoints': kpt_dict, 'path': path_labels}
