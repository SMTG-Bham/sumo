# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Module containing helper functions for dealing with band structures.

todo:
  * Extend get projections to allow specifying lm orbitals and atomic sites.
"""

import numpy as np
import itertools as it
from copy import deepcopy

from collections import defaultdict

from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.bandstructure import (BandStructure,
                                                         BandStructureSymmLine)


def get_projections_by_branches(bs, selection, normalise=None):
    """Returns orbital projections for each branch in a band structure.

    Args:
        bs (:obj:`~pymatgen.electronic_structure.bandstructure.BandStructureSymmLine`):
            The band structure.
        selection (list): A list of :obj:`tuple` or :obj:`string`
            identifying which projections to return. Projections can be
            specified by both element and orbital, for example::

                [('Sn', 's'), ('Bi', 'p'), ('S', 'p')]

            If just the element is specified then all the orbitals of that
            element are combined. For example, the following will combine
            all the S orbitals into a single projection::

                [('Bi', 's'), ('Bi', 'p'), 'S']

            Particular orbitals can also be combined, for example::

                [('Bi', 's'), ('Bi', 'p'), ('S', ('s', 'p', 'd'))]

        normalise (:obj:`str`, optional): Normalisation the projections.
            Options are:

              * ``'all'``: Projections normalised against the sum of all
                   other projections.
              * ``'select'``: Projections normalised against the sum of the
                   selected projections.
              * ``None``: No normalisation performed.

            Defaults to ``None``.

    Returns:
        list: A ``list`` of orbital projections for each branch of the band
        structure, in the same order as specified in ``selection``, with
        the format::

            [ [ {spin: projections} ], [ {spin: projections} ], ... ]

        Where spin is a :obj:`pymatgen.electronic_structure.core.Spin`
        object and projections is a :obj:`numpy.array` of::

            projections[band_index][kpoint_index]

        If there are no projections in the band structure, then an array of
        zeros is returned for each spin.
    """
    spins = bs.bands.keys()
    projections = get_projections(bs, selection, normalise=normalise)

    branches = []
    for b in bs.branches:
        s = b['start_index']
        e = b['end_index'] + 1

        branch_proj = deepcopy(projections)
        for spin, i in it.product(spins, range(len(projections))):
            branch_proj[i][spin] = projections[i][spin][:, s:e]

        branches.append(branch_proj)
    return branches


def get_projections(bs, selection, normalise=None):
    """Returns orbital projections from a band structure.

    Args:
        bs (:obj:`~pymatgen.electronic_structure.bandstructure.BandStructureSymmLine`):
            The band structure.
        selection (list): A list of :obj:`tuple` or :obj:`string`
            identifying which projections to return. Projections can be
            specified by both element and orbital, for example::

                [('Bi', 's'), ('Bi', 'p'), ('S', 'p')]

            If just the element is specified then all the orbitals of
            that element are combined. For example, the following will combine
            all the S orbitals into a single projection::

                [('Bi', 's'), ('Bi', 'p'), 'S']

            Particular orbitals can also be combined, for example::

                [('Bi', 's'), ('Bi', 'p'), ('S', ('s', 'p', 'd'))]

        normalise (:obj:`str`, optional): Normalisation the projections.
            Options are:

              * ``'all'``: Projections normalised against the sum of all
                   other projections.
              * ``'select'``: Projections normalised against the sum of the
                   selected projections.
              * ``None``: No normalisation performed.

            Defaults to ``None``.

    Returns:
        list: A ``list`` of orbital projections, in the same order as specified
        in ``selection``, with the format::

            [ {spin: projections}, {spin: projections} ... ]

        Where spin is a :obj:`pymatgen.electronic_structure.core.Spin`
        object and projections is a :obj:`numpy.array` of::

            projections[band_index][kpoint_index]

        If there are no projections in the band structure, then an array of
        zeros is returned for each spin.
    """
    spins = bs.bands.keys()
    nbands = bs.nb_bands
    nkpts = len(bs.kpoints)

    # if we are to normalise the data later we need access to all projections
    elements = bs.structure.symbol_set
    all_orbitals = ['s', 'p', 'd', 'f']

    # dictio has the form: {'el1': [s, p, d, f], 'el2': [s, p, d, f]...}
    dictio = dict(zip(elements, [all_orbitals] * len(elements)))

    # bs.get_projection_on_elements_and_orbitals() returns the data in a
    # really frustrating format, namely:
    #     {spin: [band_index][kpoint_index]{element: {orbital: projection}}}
    all_proj = bs.get_projections_on_elements_and_orbitals(dictio)

    # Make a defaultdict of defaultdicts
    dict_proj = defaultdict(lambda: defaultdict(dict))
    sum_proj = dict(zip(spins, [np.zeros((nbands, nkpts))] * len(spins)))

    # store the projections for all elements and orbitals in a useable format
    for spin, element, orbital in it.product(spins, elements, all_orbitals):

        # convert data to [nb][nk][projection]
        el_orb_proj = [[all_proj[spin][nb][nk][element][orbital]
                        for nk in range(nkpts)] for nb in range(nbands)]

        dict_proj[element][orbital][spin] = np.array(el_orb_proj)

        if normalise == 'all':
            sum_proj[spin] += el_orb_proj

    # now go through the selected orbitals and extract what's needed
    spec_proj = []
    for spec in selection:

        if isinstance(spec, str):
            # spec is just an element type, therefore sum all orbitals
            element = spec
            orbitals = all_orbitals
        else:
            element, orbitals = spec
            # even if there is only one orbital, make sure we can loop over it
            orbitals = tuple(orbitals)

        proj = dict(zip(spins, [np.zeros((nbands, nkpts))] * len(spins)))
        for spin, orbital in it.product(spins, orbitals):
            proj[spin] += dict_proj[element][orbital][spin]

            if normalise == 'select':
                sum_proj[spin] += dict_proj[element][orbital][spin]

        spec_proj.append(proj)

    if normalise:
        # to prevent warnings/errors relating to divide by zero,
        # catch warnings and surround divide with np.nan_to_num
        with np.errstate(divide='ignore', invalid='ignore'):
            for spin, i in it.product(spins, range(len(spec_proj))):
                spec_proj[i][spin] = np.nan_to_num(spec_proj[i][spin] /
                                                   sum_proj[spin])
    return spec_proj


def get_reconstructed_band_structure(list_bs, efermi=None):
    """Combine a list of band structures into a single band structure.

    This is typically very useful when you split non self consistent
    band structure runs in several independent jobs and want to merge back
    the results.

    This method will also ensure that any BandStructure objects will contain
    branches.

    Args:
        list_bs (:obj:`list` of \
        :obj:`~pymatgen.electronic_structure.bandstructure.BandStructure` \
        or :obj:`~pymatgen.electronic_structure.bandstructure.BandStructureSymmLine`):
            The band structures.
        efermi (:obj:`float`, optional): The Fermi energy of the reconstructed
            band structure. If `None`, an average of all the Fermi energies
            across all band structures is used.

    Returns:
        :obj:`pymatgen.electronic_structure.bandstructure.BandStructure` or \
        :obj:`pymatgen.electronic_structure.bandstructureBandStructureSymmLine`:
        A band structure object. The type depends on the type of the band
        structures in ``list_bs``.
    """
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
                and (not np.array_equal(kpoints[i + 1], k)
                     or not np.array_equal(kpoints[i - 1], k))):
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
            projs = [bs.projections[Spin.down][:nb_bands][dup_ids]
                     for bs in list_bs]
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
