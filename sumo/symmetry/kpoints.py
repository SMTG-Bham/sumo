# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Module providing helper functions for generating k-points along a path.
"""

import os
import sys
import math
import errno
import shutil
import logging

from pymatgen.io.vasp.inputs import Kpoints


def get_path_data(structure, mode='bradcrack', symprec=0.01, spg=None,
                  line_density=60, cart_coords=False, kpt_list=None,
                  labels=None, phonopy=False):
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
    from sumo.symmetry import (BradCrackKpath, SeekpathKpath, PymatgenKpath,
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
                                        phonopy=phonopy,
                                        cart_coords=cart_coords)
    path_str = kpath.path_string
    kpt_dict = kpath.kpoints

    logging.info('Structure information:')
    logging.info('\tSpace group number: {}'.format(kpath._spg_data['number']))

    logging.info('\tInternational symbol: {}'.format(kpath.spg_symbol))
    logging.info('\tLattice type: {}'.format(kpath.lattice_type))

    logging.info('\nk-point path:\n\t{}'.format(path_str))
    logging.info('\nk-points:')

    for label, kpoint in iter(kpt_dict.items()):
        coord_str = ' '.join(['{}'.format(c) for c in kpoint])
        logging.info('\t{}: {}'.format(label, coord_str))

    return kpath, kpoints, labels


def write_kpoint_files(filename, kpoints, labels, make_folders=False,
                       ibzkpt=None, kpts_per_split=None, directory=None,
                       cart_coords=False):
    r"""Write the k-points data to files.

    Folders are named as 'split-01', 'split-02', etc ...
    KPOINTS files are named KPOINTS_band_split_01 etc ...

    Args:
        filename (:obj:`str`): Path to VASP structure file.
        kpoints (:obj:`numpy.ndarray`): The k-point coordinates along the
            high-symmetry path. For example::

                [[0, 0, 0], [0.25, 0, 0], [0.5, 0, 0], [0.5, 0, 0.25],
                [0.5, 0, 0.5]]

        labels (:obj:`list`) The high symmetry labels for each k-point (will be
            an empty :obj:`str` if the k-point has no label). For example::

                ['\Gamma', '', 'X', '', 'Y']

        make_folders (:obj:`bool`, optional): Generate folders and copy in
            required files (INCAR, POTCAR, POSCAR, and possibly CHGCAR) from
            the current directory.
        ibzkpt (:obj:`str`, optional): Path to IBZKPT file. If set, the
            generated k-points will be appended to the k-points in this file
            and given a weight of 0. This is necessary for hybrid band
            structure calculations.
        kpts_per_split (:obj:`int`, optional): If set, the k-points are split
            into separate k-point files (or folders) each containing the number
            of k-points specified. This is useful for hybrid band structure
            calculations where it is often intractable to calculate all
            k-points in the same calculation.
        directory (:obj:`str`, optional): The output file directory.
        cart_coords (:obj:`bool`, optional): Whether the k-points are returned
            in cartesian or reciprocal coordinates. Defaults to ``False``
            (fractional coordinates).
    """
    if kpts_per_split:
        kpt_splits = [kpoints[i:i+kpts_per_split] for
                      i in range(0, len(kpoints), kpts_per_split)]
        label_splits = [labels[i:i+kpts_per_split] for
                        i in range(0, len(labels), kpts_per_split)]
    else:
        kpt_splits = [kpoints]
        label_splits = [labels]

    if cart_coords:
        coord_type = 'cartesian'
        style = Kpoints.supported_modes.Cartesian
    else:
        coord_type = 'reciprocal'
        style = Kpoints.supported_modes.Reciprocal

    kpt_files = []
    for kpt_split, label_split in zip(kpt_splits, label_splits):
        if ibzkpt:
            # hybrid calculation so set k-point weights to 0
            kpt_weights = ibzkpt.kpts_weights + [0] * len(kpt_split)
            kpt_split = ibzkpt.kpts + kpt_split
            label_split = [''] * len(ibzkpt.labels) + label_split
        else:
            # non-SCF calculation so set k-point weights to 1
            kpt_weights = [1] * len(kpt_split)

        segment = ' -> '.join([label for label in label_split if label])
        kpt_file = Kpoints(comment=segment, num_kpts=len(kpt_split),
                           kpts=kpt_split, kpts_weights=kpt_weights,
                           style=style, coord_type=coord_type,
                           labels=label_split)
        kpt_files.append(kpt_file)

    pad = int(math.floor(math.log10(len(kpt_files)))) + 2
    if make_folders:
        for i, kpt_file in enumerate(kpt_files):
            folder = 'split-{}'.format(str(i+1).zfill(pad))
            if directory:
                folder = os.path.join(directory, folder)

            try:
                os.makedirs(folder)
            except OSError as e:
                if e.errno == errno.EEXIST:
                    logging.error("\nERROR: Folders already exist, won't "
                                  "overwrite.")
                    sys.exit()
                else:
                    raise

            kpt_file.write_file(os.path.join(folder, 'KPOINTS'))
            vasp_files = [filename, "INCAR", "POTCAR", "job"]
            vasp_files += [] if ibzkpt else ['CHGCAR']
            for vasp_file in vasp_files:
                if os.path.isfile(vasp_file):
                    shutil.copyfile(vasp_file, os.path.join(folder, vasp_file))
    else:
        for i, kpt_file in enumerate(kpt_files):
            if len(kpt_files) > 1:
                kpt_filename = 'KPOINTS_band_split_{:0d}'.format(i + 1)
            else:
                kpt_filename = 'KPOINTS_band'
            if directory:
                kpt_filename = os.path.join(directory, kpt_filename)
            kpt_file.write_file(kpt_filename)


def _get_space_group_object(spg, mode):
    from pymatgen.symmetry.groups import SpaceGroup
    if spg and mode != 'bradcrack':
        logging.error("ERROR: Specifying symmetry only supported using "
                      "Bradley and Cracknell path.")
        sys.exit()
    elif spg:
        try:
            if isinstance(spg, int):
                spg = SpaceGroup.from_int_number(spg)
            else:
                spg = SpaceGroup(spg)
            logging.error("WARNING: Forcing space group not recommended, the "
                          "path is likely\nincorrect. Use at your own risk.\n")
        except ValueError:
            logging.error("ERROR: Space group not recognised.")
            sys.exit()
    return spg
