# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Module providing helper functions for generating k-points along a path.
"""

import string
import logging
import numpy as np


def get_kpoints(structure, kpoints, path, line_density=20, cart_coords=False,
                phonopy=False):
    r"""Calculate a list of k-points along the high-symmetry path.

    Adapted from
    :obj:`pymatgen.symmetry.bandstructure.HighSymmKpath.get_kpoints`.

    Args:
        structure (:obj:`~pymatgen.core.structure.Structure`): The structure.
        kpoints (dict): The high-symmetry k-point labels and their fractional
            coordinates. Formatted as ``{label: coords}``. For example::

                {'\Gamma': [0., 0., 0.], 'X': [0.5, 0. 0.]}

        path (list): The high-symmetry k-point path. Each subpath is provided
            as a list. For example, the following covers the path ``\Gamma ->
            X -> C | \Gamma -> Y``::

                [['\Gamma', 'X', 'C'], ['\Gamma', 'Y']].

        line_density (:obj:`int`, optional): Density of k-points along the
            path.
        cart_coords (:obj:`bool`, optional): Whether the k-points are
            returned in cartesian or reciprocal coordinates. Defaults to
            ``False`` (fractional coordinates).
        phonopy (:obj:`bool`, optional): Format the k-points and labels for
            use with phonopy. Defaults to ``False``.

    Returns:
            list: A list of k-points along the high-symmetry path, together
            with the high symmetry labels for each k-point. Returned as
            ``(kpoints, labels)``.

            If ``phonopy == False``, then:

                * ``kpoints`` is a :obj:`np.ndarray` of the k-point fractional
                  coordinates along the high-symmetry path.
                * ``labels`` is a :obj:`list` of the high symmetry labels for
                  each k-point (will be ``''`` if no label is set).

            If ``phonopy == True``, then:

                * todo: detail phonopy=True output.
    """
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
            nb = int(np.ceil(distance * line_density))
            sym_point_labels.extend([b[i - 1]] + [''] * (nb - 1))

            limit = nb + 1 if phonopy else nb
            kpts = [recip_lattice.get_cartesian_coords(start)
                    + float(i) / float(nb) *
                    (recip_lattice.get_cartesian_coords(end)
                     - recip_lattice.get_cartesian_coords(start))
                    for i in range(0, limit)]

            if phonopy:
                list_k_points.append(kpts)
            else:
                list_k_points.extend(kpts)

        # append last k-point to avoid repition as in pymatgen
        if not phonopy:
            sym_point_labels.append(b[-1])
            list_k_points.append(recip_lattice.get_cartesian_coords(end))

    if phonopy:
        # TODO: fix for multiple band paths
        sym_point_labels = path[0]

    if cart_coords:
        return list_k_points, sym_point_labels
    else:
        if phonopy:
            frac_k_points = [[recip_lattice.get_fractional_coords(k)
                             for k in p] for p in list_k_points]
            # TODO: fix for multiple band paths
            frac_k_points = frac_k_points
        else:
            frac_k_points = [recip_lattice.get_fractional_coords(k)
                             for k in list_k_points]
        return frac_k_points, sym_point_labels


def get_kpoints_from_list(structure, kpt_list, path_labels=None,
                          line_density=60, cart_coords=False, phonopy=False):
    """Generate the k-points along a manually specified path.

    If no labels are provided, letters from A -> Z will be used instead.

    Args:
        structure (:obj:`~pymatgen.core.structure.Structure`): The structure.
        kpt_list (list): List of k-points to use, formatted as a list of
            subpaths, each containing a list of fractional k-points. For
            example::

                [ [[0., 0., 0.], [0., 0., 0.5]],
                  [[0.5, 0., 0.], [0.5, 0.5, 0.]] ]

            Will return points along 0 0 0 -> 0 0 1/2 | 1/2 0 0 -> 1/2 1/2 0
        path_labels (:obj:`list`): The k-point labels. These should be provided
            as a :obj:`list` of :obj:`str` for each subpath of the overall
            path. For example::

                [ ['Gamma', 'Z'], ['X', 'M'] ]

            combined with the above example for ``kpt_list`` would indicate the
            path: Gamma -> Z | X -> M.
        line_density (:obj:`int`, optional): Density of k-points along the
            path.
        cart_coords (:obj:`bool`, optional): Whether the k-points are
            returned in cartesian or reciprocal coordinates. Defaults to
            ``False`` (fractional coordinates).
        phonopy (:obj:`bool`, optional): Format the k-points and labels for
            use with phonopy. Defaults to ``False``.

    Returns:
        tuple: The k-points along the high-symmetry path, together with the
        high symmetry labels for each k-point, a printable string of the
        high-symmetry path, and a dictionary mapping the path labels to the
        k-point coordinates (e.g. ``{label: coords}``). Returned as:
        ``(kpoints, labels, path_string, kpt_dict)``.

        If ``phonopy == False``, then:

            * ``kpoints`` is a :obj:`np.ndarray` of the k-point fractional
                coordinates along the high-symmetry path.
            * ``labels`` is a :obj:`list` of the high symmetry labels for
                each k-point (will be ``''`` if no label is set).

        If ``phonopy == True``, then:

            * todo: detail phonopy=True output.

    """
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

    kpoints, labels = get_kpoints(structure, kpt_dict, path_labels,
                                  line_density=line_density,
                                  cart_coords=cart_coords, phonopy=phonopy)
    path_str = ' | '.join([' -> '.join(subpath) for subpath in path_labels])
    return kpoints, labels, path_str, kpt_dict


def get_path_data(structure, mode='bradcrack', symprec=0.01, spg=None,
                  line_density=60, kpt_list=None, labels=None, phonopy=False):
    """Calculate the high-symmetry k-point path for a structure.

    The format of the returned data will be different if phonopy is True or
    False. This is because phonopy requires the labels and kpoints lists to be
    provided in a different format than kgen.

    Args:
        structure (Structure): A pymatgen structure object.
        mode (str): Sets the method used to calculate the high-symmetry path.
            Choice of 'bradcrack', 'seekpath', and 'pymatgen'.
        symprec (float): The precision used for determining cell symmetry.
        spg (str or int): The space group international number or symbol to
            override the symmetry determined by spglib. This is not recommended
            and only provided for testing purposes.
        line_density (int): The density of k-points along the path.
        kpt_list (list): Manual list of k-points to use. If kpt_list is set it
            will override the mode selection. Should be formatted as a list of
            subpaths, each containing a list of k-points. For example:
            [[[0., 0., 0.], [0., 0., 0.5]], [[0.5, 0., 0.], [0.5, 0.5, 0.]]]
        labels (list): A list of labels to use along with kpt_list. These should
            be provided as a list of subpaths, each containing a list of labels.
            For example: [['Gamma', 'Z'], ['X', 'M']], combined with the above
            kpt_list would indicate the path: Gamma -> Z | X -> M.
            If no labels are provided, letters from A -> Z will be used instead.
        phonopy (bool): Format the k-points and labels for use with phonopy.

    Returns:
        A tuple of the (kpoint path, kpoints, labels).

        If phonopy=False, then:
            `kpoints` is a np.array of the k-point fractional coordinates
                along the high-symmetry path.
            `labels`is the high symmetry labels for each k-point,
                (will be '' if no label set).

        If phonopy is set to True the format will differ:
            TODO: detail phonopy=True output.
    """
    from vaspy.symmetry import BradCrackKpath, SeekpathKpath, PymatgenKpath
    spg = _get_space_group_object(spg, mode)

    if mode == 'bradcrack':
        kpath = BradCrackKpath(structure, symprec=symprec, spg=spg)
    elif mode == 'seekpath':
        kpath = SeekpathKpath(structure, symprec=symprec)
    elif mode == 'pymatgen':
        kpath = PymatgenKpath(structure, symprec=symprec)

    if kpt_list is not None:
        kpoints, labels, path_str, kpt_dict = get_kpoints_from_list(
            structure, kpt_list, path_labels=labels, line_density=density,
            phonopy=phonopy)
    else:
        kpoints, labels = kpath.get_kpoints(line_density=line_density,
                                            phonopy=phonopy)
        path_str = kpath.path_string
        kpt_dict = kpath.kpoints

    logging.info('Structure information:'.format(structure.num_sites))
    logging.info('\tSpace group number: {}'.format(kpath._spg_data['number']))

    logging.info('\tInternational symbol: {}'.format(kpath.spg_symbol))
    logging.info('\tLattice type: {}'.format(kpath.lattice_type))

    logging.info('\nk-point path:\n\t{}'.format(path_str))
    logging.info('\nk-points:')
    for label, kpoint in iter(kpt_dict.items()):
        coord_str = ' '.join(['{}'.format(c) for c in kpoint])
        logging.info('\t{}: {}'.format(label, coord_str))

    return kpath, kpoints, labels

def _get_space_group_object(spg, mode):
    import sys
    if spg and mode != 'bradcrack':
        logging.error("ERROR: Specifying symmetry only supported using "
                      "Bradley and Cracknell path.")
        sys.exit()
    elif spg:
        try:
            if type(spg) is int:
                spg = SpaceGroup.from_int_number(spg)
            else:
                spg = SpaceGroup(spg)
            logging.error("WARNING: Forcing space group not recommended, the"
                          " path is likely\nincorrect. Use at your own risk.\n")
        except ValueError:
            logging.error("ERROR: Space group not recognised.")
            sys.exit()
    return spg
