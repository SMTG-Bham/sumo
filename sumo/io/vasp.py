import errno
import logging
import math
import os
import shutil
import sys

from pymatgen.io.vasp.inputs import Kpoints


def write_kpoint_files(
    filename,
    kpoints,
    labels,
    make_folders=False,
    ibzkpt=None,
    kpts_per_split=None,
    directory=None,
    cart_coords=False,
):
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
        kpt_splits = [
            kpoints[i : i + kpts_per_split]
            for i in range(0, len(kpoints), kpts_per_split)
        ]
        label_splits = [
            labels[i : i + kpts_per_split]
            for i in range(0, len(labels), kpts_per_split)
        ]
    else:
        kpt_splits = [kpoints]
        label_splits = [labels]

    if cart_coords:
        coord_type = "cartesian"
        style = Kpoints.supported_modes.Cartesian
    else:
        coord_type = "reciprocal"
        style = Kpoints.supported_modes.Reciprocal

    kpt_files = []
    for kpt_split, label_split in zip(kpt_splits, label_splits):
        if ibzkpt:
            # hybrid calculation so set k-point weights to 0
            kpt_weights = ibzkpt.kpts_weights + [0] * len(kpt_split)
            kpt_split = ibzkpt.kpts + kpt_split
            label_split = [""] * len(ibzkpt.labels) + label_split
        else:
            # non-SCF calculation so set k-point weights to 1
            kpt_weights = [1] * len(kpt_split)

        segment = " -> ".join([label for label in label_split if label])
        kpt_file = Kpoints(
            comment=segment,
            num_kpts=len(kpt_split),
            kpts=kpt_split,
            kpts_weights=kpt_weights,
            style=style,
            coord_type=coord_type,
            labels=label_split,
        )
        kpt_files.append(kpt_file)

    pad = int(math.floor(math.log10(len(kpt_files)))) + 2
    if make_folders:
        for i, kpt_file in enumerate(kpt_files):
            folder = f"split-{str(i + 1).zfill(pad)}"
            if directory:
                folder = os.path.join(directory, folder)

            try:
                os.makedirs(folder)
            except OSError as e:
                if e.errno == errno.EEXIST:
                    logging.error("\nERROR: Folders already exist, won't overwrite.")
                    sys.exit()
                else:
                    raise

            kpt_file.write_file(os.path.join(folder, "KPOINTS"))
            vasp_files = [filename, "INCAR", "POTCAR", "job"]
            vasp_files += [] if ibzkpt else ["CHGCAR"]
            for vasp_file in vasp_files:
                if os.path.isfile(vasp_file):
                    shutil.copyfile(vasp_file, os.path.join(folder, vasp_file))
    else:
        for i, kpt_file in enumerate(kpt_files):
            if len(kpt_files) > 1:
                kpt_filename = f"KPOINTS_band_split_{i + 1:0d}"
            else:
                kpt_filename = "KPOINTS_band"
            if directory:
                kpt_filename = os.path.join(directory, kpt_filename)
            kpt_file.write_file(kpt_filename)
