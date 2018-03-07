# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

import numpy as np


def get_reconstructed_band_structure(list_bs, efermi=None):
    """
    This method takes a list of band structures and reconstructs
    one band structure object from all of them.

    This is typically very useful when you split non self consistent
    band structure runs in several independent jobs and want to merge back
    the results.

    This script will also ensure that any BandStructure objects will contain
    branches.

    Args:
        list_bs: A list of BandStructure or BandStructureSymmLine objects.
        efermi: The Fermi energy of the reconstructed band structure. If
            None is assigned an average of all the Fermi energy in each
            object in the list_bs is used.
    Returns:
        A BandStructure or BandStructureSymmLine object (depending on
        the type of the list_bs objects)
    """
    from pymatgen.electronic_structure.core import Spin
    from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine

    if efermi is None:
        efermi = sum([b.efermi for b in list_bs]) / len(list_bs)

    kpoints = []
    labels_dict = {}
    rec_lattice = list_bs[0].lattice_rec
    nb_bands = min([list_bs[i].nb_bands for i in range(len(list_bs))])

    kpoints = np.concatenate([[k.frac_coords for k in bs.kpoints]
                              for bs in list_bs])


    dicts = [bs.labels_dict for bs in list_bs]
    labels_dict = {k: v.frac_coords for d in dicts for k, v in d.items()}

    # pymatgen band structure objects support branches. These are formed when
    # two kpoints with the same label are next to each other. This bit of code
    # will ensure that the band structure will contain branches, if it doesn't
    # already.
    dup_ids = []
    for i, k in enumerate(kpoints):
        dup_ids.append(i)
        if (tuple(k) in tuple(map(tuple, labels_dict.values()))
                and i != 0 and i != len(kpoints) - 1
                and (not np.array_equal(kpoints[i+1], k)
                or not np.array_equal(kpoints[i-1], k))):
            dup_ids.append(i)

    kpoints = kpoints[dup_ids]

    eigenvals = {}
    eigenvals[Spin.up] = np.concatenate([bs.bands[Spin.up][:nb_bands]
                                         for bs in list_bs], axis=1)
    eigenvals[Spin.up] = eigenvals[Spin.up][:, dup_ids]

    if list_bs[0].is_spin_polarized:
        eigenvals[Spin.down] = np.concatenate([bs.bands[Spin.down][:nb_bands]
                                               for bs in list_bs], axis=1)
        eigenvals[Spin.down] = eigenvals[Spin.up][:, dup_ids]

    projections = {}
    if len(list_bs[0].projections) != 0:
        projs = [bs.projections[Spin.up][:nb_bands][dup_ids] for bs in list_bs]
        projections[Spin.up] = np.concatenate(projs, axis=1)[:, dup_ids]

        if list_bs[0].is_spin_polarized:
            projs = [bs.projections[Spin.down][:nb_bands][dup_ids] for bs in list_bs]
            projections[Spin.down] = np.concatenate(projs, axis=1)[:, dup_ids]

    if isinstance(list_bs[0], BandStructureSymmLine):
        return BandStructureSymmLine(kpoints, eigenvals, rec_lattice,
                                     efermi, labels_dict,
                                     structure=list_bs[0].structure,
                                     projections=projections)
    else:
        return BandStructure(kpoints, eigenvals, rec_lattice, efermi,
                             labels_dict, structure=list_bs[0].structure,
                             projections=projections)
