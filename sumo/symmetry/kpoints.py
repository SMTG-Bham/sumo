# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Module providing helper functions for generating k-points along a path.
"""

import logging
import sys


def get_path_data(
    structure,
    mode="bradcrack",
    symprec=0.01,
    spg=None,
    line_density=60,
    cart_coords=False,
    kpt_list=None,
    labels=None,
    phonopy=False,
):
    r"""Get the k-point path, coordinates and symmetry labels for a structure.

    If a manual :obj:`list` of kpoints is supplied using the ``kpt_list``
    variable, the ``mode`` option will be ignored.

    The format of the returned data will be different if phonopy is ``True`` or
    ``False``. This is because phonopy requires the labels and kpoints to be
    provided in a different format than kgen.

    Args:
        structure (:obj:`~pymatgen.core.structure.Structure`): The structure.
        mode (:obj:`str`, optional): Method used for calculating the
            high-symmetry path. The options are:

            bradcrack
                Use the paths from Bradley and Cracknell. See [brad]_.

            pymatgen
                Use the paths from pymatgen. See [curt]_.

            latimer-munro
                Use the paths from Latimer & Munro. See [lm]_.

            seekpath
                Use the paths from SeeK-path. See [seek]_.

        symprec (:obj:`float`, optional): The tolerance for determining the
            crystal symmetry.
        spg (:obj:`~pymatgen.symmetry.groups.SpaceGroup`, optional): Space
            group used to override the symmetry determined by spglib. This is
            not recommended and only provided for testing purposes.
            This option will only take effect when ``mode = 'bradcrack'``.
        line_density (:obj:`int`, optional): Density of k-points along the
            path.
        cart_coords (:obj:`bool`, optional): Whether the k-points are returned
            in cartesian or reciprocal coordinates. Defaults to ``False``
            (fractional coordinates).
        kpt_list (:obj:`list`, optional): List of k-points to use, formatted as
            a list of subpaths, each containing a list of fractional k-points.
            For example::

                [ [[0., 0., 0.], [0., 0., 0.5]],
                  [[0.5, 0., 0.], [0.5, 0.5, 0.]] ]

            Will return points along ``0 0 0 -> 0 0 1/2 | 1/2 0 0
            -> 1/2 1/2 0``
        path_labels (:obj:`list`, optional): The k-point labels. These should
            be provided as a :obj:`list` of :obj:`str` for each subpath of the
            overall path. For example::

                [ ['Gamma', 'Z'], ['X', 'M'] ]

            combined with the above example for ``kpt_list`` would indicate the
            path: Gamma -> Z | X -> M. If no labels are provided, letters from
            A -> Z will be used instead.
        phonopy (:obj:`bool`, optional): Format the k-points and labels for
            use with phonopy. Defaults to ``False``.

    Returns:
        tuple: A tuple of a :obj:`~sumo.symmetry.kpath` object, the k-points
        along the high-symmetry path, and the k-point labels. Returned as
        ``(kpath, kpoints, labels)``.

        The type of ``kpath`` object will depend on the value of ``mode`` and
        whether ``kpt_list`` is set.

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
    from sumo.symmetry import (
        BradCrackKpath,
        CustomKpath,
        LatimerKpath,
        PymatgenKpath,
        SeekpathKpath,
    )

    spg = _get_space_group_object(spg, mode)

    if kpt_list:
        kpath = CustomKpath(structure, kpt_list, labels, symprec=symprec)
    elif mode == "bradcrack":
        kpath = BradCrackKpath(structure, symprec=symprec, spg=spg)
    elif mode == "seekpath":
        kpath = SeekpathKpath(structure, symprec=symprec)
    elif mode == "pymatgen":
        kpath = PymatgenKpath(structure, symprec=symprec)
    elif mode == "latimer-munro":
        kpath = LatimerKpath(structure, symprec=symprec)

    kpoints, labels = kpath.get_kpoints(
        line_density=line_density, phonopy=phonopy, cart_coords=cart_coords
    )
    path_str = kpath.path_string
    kpt_dict = kpath.kpoints

    logging.info("Structure information:")
    logging.info(f"\tSpace group number: {kpath._spg_data['number']}")

    logging.info(f"\tInternational symbol: {kpath.spg_symbol}")
    logging.info(f"\tLattice type: {kpath.lattice_type}")

    logging.info(f"\nk-point path:\n\t{path_str}")
    logging.info("\nk-points:")

    for label, kpoint in iter(kpt_dict.items()):
        coord_str = " ".join([f"{c}" for c in kpoint])
        logging.info(f"\t{label}: {coord_str}")

    return kpath, kpoints, labels


def _get_space_group_object(spg, mode):
    from pymatgen.symmetry.groups import SpaceGroup

    if spg and mode != "bradcrack":
        logging.error(
            "ERROR: Specifying symmetry only supported using "
            "Bradley and Cracknell path."
        )
        sys.exit()
    elif spg:
        try:
            if isinstance(spg, int):
                spg = SpaceGroup.from_int_number(spg)
            else:
                spg = SpaceGroup(spg)
            logging.error(
                "WARNING: Forcing space group not recommended, the "
                "path is likely\nincorrect. Use at your own risk.\n"
            )
        except ValueError:
            logging.error("ERROR: Space group not recognised.")
            sys.exit()
    return spg
