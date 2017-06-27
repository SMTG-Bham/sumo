# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

from itertools import product

from pymatgen.core.periodic_table import get_el_sp
from pymatgen.electronic_structure.core import Orbital

"""
This module contains helper functions for dealing with pymatgen Dos objects
"""

__author__ = "Alex Ganose"
__copyright__ = "Copyright 2017, Scanlon Materials Theory Group"
__version__ = "0.1"
__date__ = "Jun 23, 2017"


def get_pdos(dos, lm_orbitals=None, atoms=None, elements=None):
    """Gets the entire partial dos.

    Args:
        dos: A complete Dos object from a Vasprun
        lm_orbitals: A list of orbitals for which the lm decomposed
            contributions should be calculated, in the form {Element: [orbs]}
        atoms: A list of atoms of over which to sum the DOS. 0 indexed.
            If only the atom symbol is specified all the
            atoms are considered. Provided in the form {Element: [atoms]}
        elements: A dict of the elements containing which orbitals to plot

    Returns:
        The pdos in the form {Element: {Orbital: dos}}
    """
    pdos = {}
    elements = elements if elements else dos.structure.symbol_set
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
                 if (el in atoms and i in atoms[el]) or not atoms]
        lm = lm_orbitals[el] if (lm_orbitals and el in lm_orbitals) else None
        orbitals = elements[el] if elements and el in elements else None

        pdos[el] = get_element_pdos(dos, el, lm, orbitals)
    return pdos


def get_element_pdos(dos, element, lm_orbitals=None, orbitals=None):
    """Get the projected DOS for an element.

    Args:
        dos: A complete Dos object from a Vasprun
        element: The element as a string
        lm_orbitals: A list of orbitals for which the lm decomposed
            contributions should be calculated, in the form {Element: [orbs]}
        orbitals: A list of orbitals to include

    Returns:
        dict of {Orbital: Dos}
    """
    el_dos = {}
    for site in sites:
        # No we build up a list of exactly which elements we are after
        # First consider only the spd orbitals
        spd = [orb.name for orb in dos.get_element_spd_dos(element).keys() if
               ((orbitals and orb.name in orbitals) or not orbitals) and
               ((lm_orbitals and orb.name not in lm_orbitals) or not lm_orbitals)]
        # Now add the lm decomposed orbitals
        lm += [orb for orb in Orbital
               if lm_orbitals and orb.name[0] in lm_orbitals]
        for orb in spd:
            pdos = dos.get_site_spd_dos(site, orb)
            el_dos[orb] = el_dos[orb] + pdos if orb in el_dos else pdos
        for orb in lm:
            pdos = dos.get_site_orbital_dos(site, orb)
            el_dos[orb] = el_dos[orb] + pdos if orb in el_dos else pdos
    return el_dos


def write_files(dos, pdos, prefix=None, directory=None):
    """Write the VASP DOS to a series of files containing the formatted data

    Args:
        dos: A Dos or complete Dos object
        prefix: Prefix for file names
        directory: The directory in which to save files
    """
    # defining these cryptic lists makes formatting the data much easier later
    if len(dos.densities) == 1:
        sdata = [[Spin.up, 1, '']]
    else:
        sdata = [[Spin.up, 1, '(up)'], [Spin.down, -1, '(down)']]

    header = ['#energy']
    tdos_data = [dos.energies]
    for spin, sign, label in sdata:
        header.append('dos{}'.format(label))
        tdos_data.append(dos.densities[spin] * sign)
    tdos_data = np.stack(tdos_data, axis=1)

    filename = "{}_total_dos.dat".format(prefix) if prefix else 'total_dos.dat'
    if directory:
        filename = os.path.join(directory, filename)
    np.savefile(filename, tdos_data, header=" ".join(header))

    spin = len(dos.densities)
    for el, el_pdos in pdos.iteritems():
        header = ['#energy']
        pdos_data = [dos.energies]
        for orb, spin, sign, label in product(sort_orbitals(el_pdos), sdata):
            header.append('{}{}'.format(orb, label))
            tdos_data.append(dos.densities[spin] * sign)
        pdos_data = np.stack(pdos_data, axis=1)

        if prefix:
            filename = '{}_{}_dos.dat'.format(prefix, el)
        else:
            filename = '{}_dos.dat'.format(el)
        if directory:
            filename = os.path.join(directory, filename)
        np.savefile(filename, pdos_data, header=" ".join(header))


def sort_orbitals(element_pdos):
    """Sort the orbitals of an element's partial DOS

    Sorts the orbitals based on a standard format. E.g. s -> p -> d.
    Will also sort lm decomposed orbitals. This is useful for plotting/saving.

    Args:
        element_pdos: An element pdos in the form {Orbital: Dos}

    Returns:
        A list of the sorted orbitals
    """
    sorted_orbitals = ['s', 'p', 'py', 'pz', 'px',
                       'd', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2',
                       'f', 'f_3', 'f_2', 'f_1', 'f_0', 'f1', 'f2', 'f3']
    sorted_keys = []
    unsorted_keys = element_pdos.keys()
    for key in sorted_orbitals:
        if key in unsorted_keys:
            sorted_keys.append(key)
    return sorted_keys
