# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Module providing helper functions for generating k-points along a path.
"""

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
            # for VASP we label every k-point. If a k-point has no
            # high-symmetry label then just use an empty string.
            sym_point_labels.append(b[-1])
            list_k_points.append(recip_lattice.get_cartesian_coords(end))

    if phonopy:
        #sym_point_labels = [y for x in path for y in x]
        # phonopy needs the labels in a funny format. For example if the
        # route is  X -> Y | Z -> R, path will be [['X', 'Y'], ['Z', 'R'],
        # the labels shoudl be ['X', 'Z', 'R']
        sym_point_labels = []
        for i, path_branch in enumerate(path):
            for n, label in enumerate(path_branch):
                if i != 0 and n == 0:
                    sym_point_labels[-1] += " | {}".format(label)
                    #sym_point_labels[-1] = "{}".format(label)
                    #sym_point_labels.append(label)
                    pass
                else:
                    sym_point_labels.append(label)

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
    import logging
    from vaspy.symmetry import (BradCrackKpath, SeekpathKpath, PymatgenKpath,
                                CustomKpath)
    spg = _get_space_group_object(spg, mode)

    if kpt_list:
        kpath = CustomKpath(structure, kpt_list, labels, symprec=symprec)
    elif mode == 'bradcrack':
        kpath = BradCrackKpath(structure, symprec=symprec, spg=spg)
    elif mode == 'seekpath':
        kpath = SeekpathKpath(structure, symprec=symprec)
    elif mode == 'pymatgen':
        kpath = PymatgenKpath(structure, symprec=symprec)

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
