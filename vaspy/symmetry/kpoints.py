# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

import numpy as np


def get_kpoints(structure, kpoints, path, line_density=20, cart_coords=False,
                phonopy=False):
    """Calculate a list of k-points along the high-symmetry path.

    Adapted from pymatgen.symmetry.bandstructure

    Args:
        structure (Structure): Pymatgen structure object.
        kpoints (dict): The high-symmetry k-point labels and their coordinates
            as {label: coords}.
        path (list): The high-symmetry k-point path. Each subpath is provided
            as a list. E.g. [['A', 'B'], ['C', 'D']].
        line_density (int): The density of k-points along the path.
        cart_coords (bool): Whether the k-points are returned in cartesian
            or reciprocal coordinates.
        phonopy (bool): Format the k-points and labels for use with phonopy.

    Returns:
        A list k-points along the high-symmetry path, together with the
        high symmetry labels for each k-point. Returned as: (kpoints, labels).
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
        structure (Structure): A pymatgen structure object.
        kpt_list (list): Manual list of k-points to use. If kpt_list is set it
            will override the mode selection. Should be formatted as a list of
            subpaths, each containing a list of k-points. For example:
            [[[0., 0., 0.], [0., 0., 0.5]], [[0.5, 0., 0.], [0.5, 0.5, 0.]]]
        path_labels (list): A list of labels to use along with kpt_list. These
            should be provided as a list of subpaths, each containing a list of
            labels. For example: [['Gamma', 'Z'], ['X', 'M']], combined with
            the above kpt_list would indicate the path: Gamma -> Z | X -> M.
        line_density (int): The density of k-points along the path.
        cart_coords (bool): Whether the k-points are returned in cartesian
            or reciprocal coordinates.
        phonopy (bool): Format the k-points and labels for use with phonopy.

    Returns:
        A list k-points along the high-symmetry path, together with the
        high symmetry labels for each k-point, a printable string of the
        high-symmetry path, and a dictionary mapping the path labels to the
        k-point coordinates (e.g. {label: coords}). Returned as:
        (kpoints, labels, path_string, kpt_dict).
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
