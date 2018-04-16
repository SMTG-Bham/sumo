# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Module containing class for generating custom k-point paths.
"""

from itertools import chain

from sumo.symmetry import Kpath


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

        if path_labels is None:
            path_labels = self._auto_kpath_labels(kpt_list)

        try:
            assert len(kpt_list) == len(path_labels)
            for kpt_segment, path_segment in zip(kpt_list, path_labels):
                assert len(kpt_segment) == len(path_segment)
        except AssertionError:
            raise ValueError("kpt_list and path_labels are not consistent")

        flat_kpts = chain(*kpt_list)
        flat_path_labels = chain(*path_labels)
        kpt_dict = {label: kpt for label, kpt in zip(flat_path_labels,
                                                     flat_kpts)}

        self._kpath = {'kpoints': kpt_dict, 'path': path_labels}

    @staticmethod
    def _auto_kpath_labels(kpt_list):
        """Get a default set of labels (1), (2), (3)... for a k-point path

        Repeated points will be identified and the labels re-used.

        Args:
            kpt_list (list): Nested list representing k-point path segments,
                e.g.::

                  [[[0., 0., 0.], [0., 0., 0.5], [0., 0.5, 0.5]],
                   [[0.5, 0.5, 0.], [0., 0., 0.]]]

        Returns:
            list: Corresponding nested list of labels, e.g.::

              [['(1)', '(2)', '(3)'], ['(4)', '(1)']]

        """

        # Build dict of labels
        label_i = 1
        kpt_labels = {}
        for kpt in chain(*kpt_list):
            if tuple(kpt) in kpt_labels:
                continue
            else:
                kpt_labels.update({tuple(kpt): '({})'.format(label_i)})
                label_i += 1

        # Read out into nested lists
        kpath_labels = [[kpt_labels[tuple(kpt)] for kpt in segment]
                        for segment in kpt_list]

        return kpath_labels
