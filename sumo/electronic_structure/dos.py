# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Module containing helper functions for dealing with
:obj:`~pymatgen.electronic_structure.dos.Dos` and
:obj:`~pymatgen.electronic_structure.dos.CompleteDos` objects.
"""

from __future__ import unicode_literals

import os
import logging
import numpy as np

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.electronic_structure.core import Orbital, Spin


def load_dos(vasprun, elements=None, lm_orbitals=None, atoms=None,
             gaussian=None, total_only=False, log=False,
             adjust_fermi=True):
    """Load a vasprun and extract the total and projected density of states.

    Args:
        vasprun (str): Path to a vasprun.xml or vasprun.xml.gz file or
            a :obj:`pymatgen.io.vasp.outputs.Vasprun` object.
        elements (:obj:`dict`, optional): The elements and orbitals to extract
            from the projected density of states. Should be provided as a
            :obj:`dict` with the keys as the element names and corresponding
            values as a :obj:`tuple` of orbitals. For example, the following
            would extract the Bi s, px, py and d orbitals::

                {'Bi': ('s', 'px', 'py', 'd')}

            If an element is included with an empty :obj:`tuple`, all orbitals
            for that species will be extracted. If ``elements`` is not set or
            set to ``None``, all elements for all species will be extracted.
        lm_orbitals (:obj:`dict`, optional): The orbitals to decompose into
            their lm contributions (e.g. p -> px, py, pz). Should be provided
            as a :obj:`dict`, with the elements names as keys and a
            :obj:`tuple` of orbitals as the corresponding values. For example,
            the following would be used to decompose the oxygen p and d
            orbitals::

                {'O': ('p', 'd')}

        atoms (:obj:`dict`, optional): Which atomic sites to use when
            calculating the projected density of states. Should be provided as
            a :obj:`dict`, with the element names as keys and a :obj:`tuple` of
            :obj:`int` specifying the atomic indices as the corresponding
            values. The elemental projected density of states will be summed
            only over the atom indices specified. If an element is included
            with an empty :obj:`tuple`, then all sites for that element will
            be included. The indices are 0 based for each element specified in
            the POSCAR. For example, the following will calculate the density
            of states for the first 4 Sn atoms and all O atoms in the
            structure::

                {'Sn': (1, 2, 3, 4), 'O': (, )}

            If ``atoms`` is not set or set to ``None`` then all atomic sites
            for all elements will be considered.
        gaussian (:obj:`float`, optional): Broaden the density of states using
            convolution with a gaussian function. This parameter controls the
            sigma or standard deviation of the gaussian distribution.
        total_only (:obj:`bool`, optional): Only extract the total density of
            states. Defaults to ``False``.
        log (:obj:`bool`): Print logging messages. Defaults to ``False``.
        adjust_fermi (:obj:`bool`, optional): Shift the Fermi level to sit at
            the valence band maximum (does not affect metals).

    Returns:
        dict: The total and projected density of states. Formatted as a
        :obj:`tuple` of ``(dos, pdos)``, where ``dos`` is a
        :obj:`~pymatgen.electronic_structure.dos.Dos` object containing the
        total density of states and ``pdos`` is a :obj:`dict` of
        :obj:`dict` mapping the elements and their orbitals to
        :obj:`~pymatgen.electronic_structure.dos.Dos` objects. For example::

            {
                'Bi': {'s': Dos, 'p': Dos ... },
                'S': {'s': Dos}
            }
    """
    if isinstance(vasprun, str):
        vr = Vasprun(vasprun)
    else:
        vr = vasprun

    band = vr.get_band_structure()
    dos = vr.complete_dos

    if band.is_metal():
        if log:
            logging.info('System is metallic')
        zero_point = vr.efermi
    else:
        if log:
            logging.info('Band gap: {:.3f}'.
                         format(band.get_band_gap()['energy']))
            logging.info('DOS band gap: {:.3f}'.format(dos.get_gap()))
        zero_point = band.get_vbm()['energy']

    if adjust_fermi:
        dos.efermi -= dos.efermi - zero_point

    if vr.parameters['ISMEAR'] in [-1, 0, 1]:
        dos.energies -= vr.parameters['SIGMA']

    if gaussian:
        dos.densities = dos.get_smeared_densities(gaussian)
        for site in dos.pdos:
            for orbital in dos.pdos[site]:
                dos.pdos[site][orbital] = dos.get_site_orbital_dos(
                    site, orbital).get_smeared_densities(gaussian)

    if vr.parameters['LSORBIT']:
        # pymatgen includes the spin down channel for SOC calculations, even
        # though there is no density here. We remove this channel so the
        # plotting is easier later on.
        del dos.densities[Spin.down]
        for site in dos.pdos:
            for orbital in dos.pdos[site]:
                del dos.pdos[site][orbital][Spin.down]

    pdos = {}
    if not total_only:
        pdos = get_pdos(dos, lm_orbitals=lm_orbitals, atoms=atoms,
                        elements=elements)
    return dos, pdos


def get_pdos(dos, lm_orbitals=None, atoms=None, elements=None):
    """Extract the projected density of states from a CompleteDos object.

    Args:
        dos (:obj:`~pymatgen.electronic_structure.dos.CompleteDos`): The
            density of states.
        elements (:obj:`dict`, optional): The elements and orbitals to extract
            from the projected density of states. Should be provided as a
            :obj:`dict` with the keys as the element names and corresponding
            values as a :obj:`tuple` of orbitals. For example, the following
            would extract the Bi s, px, py and d orbitals::

                {'Bi': ('s', 'px', 'py', 'd')}

            If an element is included with an empty :obj:`tuple`, all orbitals
            for that species will be extracted. If ``elements`` is not set or
            set to ``None``, all elements for all species will be extracted.
        lm_orbitals (:obj:`dict`, optional): The orbitals to decompose into
            their lm contributions (e.g. p -> px, py, pz). Should be provided
            as a :obj:`dict`, with the elements names as keys and a
            :obj:`tuple` of orbitals as the corresponding values. For example,
            the following would be used to decompose the oxygen p and d
            orbitals::

                {'O': ('p', 'd')}

        atoms (:obj:`dict`, optional): Which atomic sites to use when
            calculating the projected density of states. Should be provided as
            a :obj:`dict`, with the element names as keys and a :obj:`tuple` of
            :obj:`int` specifying the atomic indices as the corresponding
            values. The elemental projected density of states will be summed
            only over the atom indices specified. If an element is included
            with an empty :obj:`tuple`, then all sites for that element will
            be included. The indices are 0 based for each element specified in
            the POSCAR. For example, the following will calculate the density
            of states for the first 4 Sn atoms and all O atoms in the
            structure::

                {'Sn': (1, 2, 3, 4), 'O': (, )}

            If ``atoms`` is not set or set to ``None`` then all atomic sites
            for all elements will be considered.

    Returns:
        dict: The projected density of states. Formatted as a :obj:`dict` of
        :obj:`dict` mapping the elements and their orbitals to
        :obj:`~pymatgen.electronic_structure.dos.Dos` objects. For example::

            {
                'Bi': {'s': Dos, 'p': Dos ... },
                'S': {'s': Dos}
            }
    """
    if not elements:
        symbols = dos.structure.symbol_set
        elements = dict(zip(symbols, [None] * len(symbols)))
    pdos = {}
    for el in elements:
        if atoms and el not in atoms:
            continue

        # select which sites to consider, if no sites were specified then
        # select all. Make a list of the sites of particular elements first
        # due to the dosplot atoms list specification (e.g. starts at 0 for
        # each element
        element_sites = [site for site in dos.structure.sites
                         if site.specie == get_el_sp(el)]
        sites = [site for i, site in enumerate(element_sites)
                 if not atoms or (el in atoms and i in atoms[el])]
        lm = lm_orbitals[el] if (lm_orbitals and el in lm_orbitals) else None
        orbitals = elements[el] if elements and el in elements else None

        pdos[el] = get_element_pdos(dos, el, sites, lm, orbitals)
    return pdos


def get_element_pdos(dos, element, sites, lm_orbitals=None, orbitals=None):
    """Get the projected density of states for an element.

    Args:
        dos (:obj:`~pymatgen.electronic_structure.dos.CompleteDos`): The
            density of states.
        element (str): Element symbol. E.g. 'Zn'.
        sites (tuple): The atomic indices over which to sum the density of
            states, as a :obj:`tuple`. Indices are zero based for each
            element. For example, ``(0, 1, 2)`` will sum the density of states
            for the 1st, 2nd and 3rd sites of the element specified.
        lm_orbitals (:obj:`tuple`, optional): The orbitals to decompose into
            their lm contributions (e.g. p -> px, py, pz). Should be provided
            as a :obj:`tuple` of :obj:`str`. For example, ``('p')``, will
            extract the projected density of states for the px, py, and pz
            orbitals. Defaults to ``None``.
        orbitals (:obj:`tuple`, optional): The orbitals to extract from the
            projected density of states. Should be provided as a :obj:`tuple`
            of :obj:`str`. For example, ``('s', 'px', 'dx2')`` will extract the
            s, px, and dx2 orbitals, only. If ``None``, all orbitals will be
            extracted. Defaults to ``None``.

    Returns:
        dict: The projected density of states. Formatted as a :obj:`dict`
        mapping the orbitals to :obj:`~pymatgen.electronic_structure.dos.Dos`
        objects. For example::

            {
                's': Dos,
                'p': Dos
            }
    """
    el_dos = {}
    for site in sites:
        # build a list of which orbitals we are after
        # start with s, p, and d orbitals only
        spd = [orb for orb in dos.get_element_spd_dos(element).keys() if
               ((orbitals and orb.name in orbitals) or not orbitals) and
               ((lm_orbitals and orb.name not in lm_orbitals) or
                not lm_orbitals)]

        # now add any lm decomposed orbitals
        lm = [orb for orb in Orbital
              if lm_orbitals and orb.name[0] in lm_orbitals]

        # extract the data
        for orb in spd:
            pdos = dos.get_site_spd_dos(site)[orb]
            el_dos[orb.name] = (el_dos[orb.name] + pdos if orb.name in el_dos
                                else pdos)

        for orb in lm:
            pdos = dos.get_site_orbital_dos(site, orb)
            el_dos[orb.name] = (el_dos[orb.name] + pdos if orb.name in el_dos
                                else pdos)
    return el_dos


def write_files(dos, pdos, prefix=None, directory=None, zero_to_efermi=True):
    """Write the density of states data to disk.

    Args:
        dos (:obj:`~pymatgen.electronic_structure.dos.Dos` or \
             :obj:`~pymatgen.electronic_structure.dos.CompleteDos`): The total
            density of states.
        pdos (dict): The projected density of states. Formatted as a
            :obj:`dict` of :obj:`dict` mapping the elements and their orbitals
            to :obj:`~pymatgen.electronic_structure.dos.Dos` objects. For
            example::

                {
                    'Bi': {'s': Dos, 'p': Dos},
                    'S': {'s': Dos}
                }

        prefix (:obj:`str`, optional): A prefix for file names.
        directory (:obj:`str`, optional): The directory in which to save files.
        zero_to_efermi (:obj:`bool`, optional): Normalise the energy such
             that the Fermi level is set as 0 eV.
    """
    # defining these cryptic lists makes formatting the data much easier later
    if len(dos.densities) == 1:
        sdata = [[Spin.up, 1, '']]
    else:
        sdata = [[Spin.up, 1, '(up)'], [Spin.down, -1, '(down)']]

    header = ['energy']
    eners = dos.energies - dos.efermi if zero_to_efermi else dos.energies
    tdos_data = [eners]
    for spin, sign, label in sdata:
        header.append('dos{}'.format(label))
        tdos_data.append(dos.densities[spin] * sign)
    tdos_data = np.stack(tdos_data, axis=1)

    filename = "{}_total_dos.dat".format(prefix) if prefix else 'total_dos.dat'
    if directory:
        filename = os.path.join(directory, filename)
    np.savetxt(filename, tdos_data, header=" ".join(header))

    spin = len(dos.densities)
    for el, el_pdos in pdos.items():
        header = ['energy']
        pdos_data = [eners]
        for orb in sort_orbitals(el_pdos):
            for spin, sign, label in sdata:
                header.append('{}{}'.format(orb, label))
                pdos_data.append(el_pdos[orb].densities[spin] * sign)
        pdos_data = np.stack(pdos_data, axis=1)

        if prefix:
            filename = '{}_{}_dos.dat'.format(prefix, el)
        else:
            filename = '{}_dos.dat'.format(el)
        if directory:
            filename = os.path.join(directory, filename)
        np.savetxt(filename, pdos_data, header=" ".join(header))


def sort_orbitals(element_pdos):
    """Sort the orbitals of an element's projected density of states.

    Sorts the orbitals based on a standard format. E.g. s < p < d.
    Will also sort lm decomposed orbitals. This is useful for plotting/saving.

    Args:
        element_pdos (dict): An element's pdos. Should be formatted as a
            :obj:`dict` of ``{orbital: dos}``. Where dos is a
            :obj:`~pymatgen.electronic_structure.dos.Dos` object. For example::

                {'s': dos, 'px': dos}

    Returns:
        list: The sorted orbitals.
    """
    sorted_orbitals = ['s', 'p', 'py', 'pz', 'px',
                       'd', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2',
                       'f', 'f_3', 'f_2', 'f_1', 'f0', 'f1', 'f2', 'f3']
    unsorted_keys = element_pdos.keys()

    sorted_keys = []
    for key in sorted_orbitals:
        if key in unsorted_keys:
            sorted_keys.append(key)

    return sorted_keys
