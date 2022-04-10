# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
This module provides helper functions for dealing with phonopy.
"""

import logging
import sys

from phonopy import Phonopy, file_IO
from phonopy.structure.cells import determinant
from phonopy.units import Bohr, Hartree, VaspToTHz
from pymatgen.io.phonopy import get_phonopy_structure


def load_phonopy(
    filename,
    structure,
    dim,
    symprec=1e-5,
    primitive_matrix=None,
    factor=VaspToTHz,
    symmetrise=True,
    born=None,
    write_fc=False,
):
    """Load phonopy output and return an ``phonopy.Phonopy`` object.

    Args:
        filename (str): Path to phonopy output. Can be any of ``FORCE_SETS``,
            ``FORCE_CONSTANTS``, or ``force_constants.hdf5``.
        structure (:obj:`~pymatgen.core.structure.Structure`): The unitcell
            structure.
        dim (list): The supercell size, as a :obj:`list` of :obj:`float`.
        symprec (:obj:`float`, optional): The tolerance for determining the
            crystal symmetry.
        primitive_matrix (:obj:`list` or :obj:`str`, optional): The
            transformation matrix from the conventional to primitive cell. Only
            required when the conventional cell was used as the starting
            structure. Should be provided as a 3x3 :obj:`list` of :obj:`float`.
            Alternatively the str 'auto' may be provided.
        factor (:obj:`float`, optional): The conversion factor for phonon
            frequency. Defaults to :obj:`phonopy.units.VaspToTHz`.
        symmetrise (:obj:`bool`, optional): Symmetrise the force constants.
            Defaults to ``True``.
        born (:obj:`str`, optional): Path to file containing Born effective
            charges. Should be in the same format as the file produced by the
            ``phonopy-vasp-born`` script provided by phonopy.
        write_fc (:obj:`bool` or :obj:`str`,  optional): Write the force
            constants to disk. If ``True``, a ``FORCE_CONSTANTS`` file will be
            written. Alternatively, if set to ``"hdf5"``, a
            ``force_constants.hdf5`` file will be written. Defaults to
            ``False`` (force constants not written).
    """
    unitcell = get_phonopy_structure(structure)
    num_atom = unitcell.get_number_of_atoms()
    num_satom = determinant(dim) * num_atom

    phonon = Phonopy(
        unitcell, dim, primitive_matrix=primitive_matrix, factor=factor, symprec=symprec
    )

    if "FORCE_CONSTANTS" in filename or ".hdf5" in filename:
        # if force constants exist, use these to avoid recalculating them
        if ".hdf5" in filename:
            fc = file_IO.read_force_constants_hdf5(filename)

        elif "FORCE_CONSTANTS" in filename:
            fc = file_IO.parse_FORCE_CONSTANTS(filename=filename)

        if fc.shape[0] != num_satom:
            msg = (
                "\nNumber of atoms in supercell is not consistent with the "
                "matrix shape of\nforce constants read from {}.\nPlease"
                "carefully check --dim."
            )
            logging.error(msg.format(filename))
            sys.exit()
        phonon.force_constants = fc

    elif "FORCE_SETS" in filename:
        # load the force sets from file and calculate force constants
        fs = file_IO.parse_FORCE_SETS(filename=filename)

        if fs["natom"] != num_satom:
            msg = (
                "\nNumber of atoms in supercell is not consistent with the "
                "the data in FORCE_SETS\nPlease carefully check --dim."
            )
            logging.error(msg.format(filename))
            sys.exit()

        phonon.dataset = fs

        logging.info("Calculating force constants...")
        phonon.produce_force_constants()

    if symmetrise:
        phonon.symmetrize_force_constants()

    if born:
        # load born parameters from a file
        nac_params = file_IO.parse_BORN(
            phonon._primitive, symprec=symprec, filename=born
        )
        # set the nac unit conversion factor manual,  specific to VASP
        nac_params["factor"] = Hartree * Bohr
        phonon.nac_params = nac_params

    if write_fc == "hdf5":
        file_IO.write_force_constants_to_hdf5(phonon.force_constants)
        logging.info("Force constants written to force_constants.hdf5.")

    elif write_fc:
        file_IO.write_FORCE_CONSTANTS(phonon.force_constants)
        logging.info("Force constants written to FORCE_CONSTANTS.")

    return phonon
