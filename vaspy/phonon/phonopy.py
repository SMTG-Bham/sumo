# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

import logging

from phonopy import Phonopy, file_IO
from phonopy.structure.cells import print_cell, determinant
from phonopy.units import VaspToTHz, Hartree, Bohr

from pymatgen.io.phonopy import get_phonopy_structure

def load_phonopy(filename, structure, dim, symprec=0.01, primitive_matrix=None, 
                 factor=VaspToTHz, symmetrise=True, born=None, write_fc=False):
    unitcell = get_phonopy_structure(structure)
    num_atom = unitcell.get_number_of_atoms()
    num_satom = determinant(dim) * num_atom

    phonon = Phonopy(unitcell, dim, primitive_matrix=primitive_matrix,
                     factor=factor, symprec=symprec)

    if 'FORCE_CONSTANTS' == filename or '.hdf5' in filename:
        if '.hdf5' in filename:
            fc = file_IO.read_force_constants_hdf5(filename)
        elif 'FORCE_CONSTANTS' == filename:
            fc = file_IO.parse_FORCE_CONSTANTS(filename=filename)

        if fc.shape[0] != num_satom:
            msg = ("\nNumber of atoms in supercell is not consistent with the "
                   "matrix shape of\nforce constants read from {}.\nPlease"
                   "carefully check --dim.")
            logging.error(msg.format(filename))
            sys.exit()

        phonon.set_force_constants(fc)

    elif 'FORCE_SETS' == filename:
        fs = file_IO.parse_FORCE_SETS()
        if fs['natom'] != num_satom:
            msg = ("\nNumber of atoms in supercell is not consistent with the "
                   "the data in FORCE_SETS\nPlease carefully check --dim.")
            logging.error(msg.format(filename))
            sys.exit()

        phonon.set_displacement_dataset(fs)
        logging.info("Calculating force constants...")
        phonon.produce_force_constants()

    if born:
        nac_params = file_IO.parse_BORN(atoms, filename=born)
        # Have to set the nac unit conversion factor manually
        # this is specific to the VASP code (I got these values from looking at
        # the phonopy source code in file 'scripts/phonopy.py')
        nac_params['factor'] = Hartree * Bohr
        phonon.set_nac_params(nac_params)

    if symmetrise:
        phonon.symmetrize_force_constants()

    if write_fc == 'hdf5':
        file_IO.write_force_constants_to_hdf5(phonon.get_force_constants())
        logging.info("Force constants written to force_constants.hdf5.")
    elif write_fc:
        file_IO.write_FORCE_CONSTANTS(phonon.get_force_constants())
        logging.info("Force constants written to force_constants.hdf5.")

    return phonon
