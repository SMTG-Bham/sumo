# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Module containing class for generating custom k-point paths.
"""

import string
import numpy as np

from vaspy.symmetry import Kpath


class CustomKpath(Kpath):
    """Class to generate k-points along custom k-point paths.

    Args:
        structure (:obj:`~pymatgen.core.structure.Structure`): The structure.
        symprec (:obj:`float`, optional): The tolerance for determining the
            crystal symmetry.
        kpt_list (list): List of k-points to use, formatted as a list of
            subpaths, each containing a list of fractional k-points. For
            example::

                [ [[0., 0., 0.], [0., 0., 0.5]],
                  [[0.5, 0., 0.], [0.5, 0.5, 0.]] ]

            Will return points along ``0 0 0 -> 0 0 1/2 | 1/2 0 0
            -> 1/2 1/2 0``
        path_labels (:obj:`list`): The k-point labels. These should be provided
            as a :obj:`list` of :obj:`str` for each subpath of the overall
            path. For example::

                [ ['Gamma', 'Z'], ['X', 'M'] ]

            combined with the above example for ``kpt_list`` would indicate the
            path: Gamma -> Z | X -> M. If no labels are provided, letters from
            A -> Z will be used instead.

    Attributes:
        prim (:obj:`~pymatgen.core.structure.Structure`): The standardised
            primitive cell structure for the generated k-point path.
        conv (:obj:`~pymatgen.core.structure.Structure`): The standardised
            conventional cell structure.
    """

    def __init__(self, structure, kpt_list, path_labels=None, symprec=1e-3):
        Kpath.__init__(self, structure, symprec=symprec)

        flat_kpts = [x for kpts in kpt_list for x in kpts]
        if path_labels is None:
            flat_path_labels = [s for s, x in
                                zip(string.ascii_uppercase, flat_kpts)]
        else:
            try:
                assert len(kpt_list) == len(path_labels)
                for kpt_segment, path_segment in zip(kpt_list, path_labels):
                    assert len(kpt_segment) == len(path_segment)
            except AssertionError:
                raise ValueError("kpt_list and path_labels are not consistent")

            flat_path_labels = [x for labels in path_labels for x in labels]

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
