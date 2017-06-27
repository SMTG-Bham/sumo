"""
This module defines classes to represent the density of states, etc.
"""
__author__ = "Alex Ganose"
__copyright__ = "Copyright 2012, David Scanlon Materials Theory"
__version__ = "0.1"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "Sept 1, 2015"

import numpy as np


class Dos(object):
    """Basic DoS object containing just the energies and densities

    Args:
        efermi: The Fermi level
        energies: A list of energies
        densities ({Spin: np.array}): The densities for each spin

    Attributes:
        efermi: The Fermi level
        energies: The sequence of energies
        densities: A dict of spin densitites, e.g {0: [...], 1: [...]}
    """

    lm_orbitals = ['s',
                   'py', 'pz', 'px',
                   'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2']
    lm_f_orbitals = ['s',
                     'py', 'pz', 'px',
                     'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2',
                     'f-3', 'f-2', 'f-1', 'f-0', 'f1', 'f2', 'f3']

    orbitals = ['s', 'p', 'd']
    f_orbitals = ['s', 'p', 'd', 'f']

    def __init__(self, efermi, energies, densities):
        self.efermi = efermi
        self.energies = np.array(energies)
        self.densities = {k: np.array(d) for k, d in densities.items()}

    def get_densities(self, spin=None):
        """Returns the density of states for a particular spin.

        Args:
            spin: The spin where 0 = Up and 1 is Down.

        Returns:
            Returns the density of states for a particular spin. If the spin is
            None, the sum of all spins is returned.
        """
        if spin is None:
            if self.spin == 2:
                result = self.densities[0] + self.densities[1]
            else:
                result = self.densities[0]
        else:
            result = self.densities[spin]
        return result

    def get_smeared_dos(self, sigma):
        """Gaussian smear the Dos.

        Args:
            sigma: Std dev of Gaussian smearing function.

        Returns:
            Dict of Gaussian-smeared densities
        """
        smeared_dens = smear_dens(self.energies, self.densities, sigma)
        return Dos(self.efermi, self.energies, smeared_dens)

    def get_cbm_vbm(self, tol=0.001, abs_tol=False, spin=None):
        """Find the CBM and VBM based on the Fermi level and a occ tolerance.

        If the system is a metal then the VBM and the CBM both equal the fermi
        level

        Args:
            tol: Tolerance used to decide if the states are occupied
            abs_tol: Use an absolute tolerance (true) or relative (False)
            spin: Find the gap both up and down (None), Spin up only (0) or
                spin down only (1)

        Returns:
            (cbm, vbm): Floats in eV
        """
        tdos = self.get_densities(spin)
        if not abs_tol:
            tol = tol * tdos.sum() / tdos.shape[0]

        i_fermi = 0
        while self.energies[i_fermi] <= self.efermi:
            i_fermi += 1

        # work backwards until tolerance is reached
        i_gap_start = i_fermi
        while i_gap_start - 1 >= 0 and tdos[i_gap_start - 1] <= tol:
            i_gap_start -= 1

        # work forwards until tolerance is reached
        i_gap_end = i_gap_start
        while i_gap_end < tdos.shape[0] and tdos[i_gap_end] <= tol:
            i_gap_end += 1
        i_gap_end -= 1
        return self.energies[i_gap_end], self.energies[i_gap_start]

    def get_band_gap(self, tol=0.001, abs_tol=False, spin=None):
        """Find the fundamental band gap

        Args:
            tol: Tolerance used to decide if the states are occupied
            abs_tol: Use an absolute tolerance (true) or relative (False)
            spin: Find the gap both up and down (None), Spin up only (0) or
                spin down only (1)

        Returns:
            The band gap in eV
        """
        (cbm, vbm) = self.get_cbm_vbm(tol, abs_tol, spin)
        return max(cbm - vbm, 0.0)

    def to_file(self, prefix=""):
        """Write the DOS to a file containing the formatted data.

        Args:
            prefix: Prefix for the file names
        """
        filename = prefix + "_total_dos.dat"
        head = '#Energy  ' + fmt_head(self.spin, ['total_DOS'], '(down)')
        form = "%.3f  " + fmt_head(self.spin, ['%G'], '-', prefix=True)
        with open(filename, 'w') as f:
            f.write(head + "\n")
            for i, energy in enumerate(self.energies):
                dos_line = [energy]
                for spin in range(self.spin):
                    dos_line.append(self.densities[spin][i])
                f.write((form + '\n') % tuple(dos_line))

    @property
    def spin(self):
        """The number of spins in the DOS

        Returns:
            spin: 2 if spin polarised and 1 if not
        """
        return 2 if len(self.densities) > 1 else 1


class VaspDos(Dos):
    """All encompasing DOS object, with access to the pdos and smearing

    Args:
        efermi: The Fermi level
        energies: A list of energies
        tdos ({Spin: np.array}): The total DOS for each spin
        pdoss: The partial density of states supplied as {Element:[{Orbital:{
            Spin:Densities}]}
        structure: The structure, in the form of {El: natoms} TODO: Update this


    Attributes:
        efermi: The Fermi energy
        energies: The sequence of energies
        densities: A dict of spin densitites, e.g {0: [...], 1: [...]}
        pdos: The partial density of states in the form {Element:[{Orbital:{
            Spin:Densities}]}
    """

    def __init__(self, efermi, energies, tdos, pdoss, structure):
        super(VaspDos, self).__init__(
            efermi=efermi, energies=energies,
            densities={k: np.array(d) for k, d in tdos.iteritems()})
        self.pdos = pdoss
        self.structure = structure
        self.remake_total()

    def remake_total(self):
        tdos = {}
        for el, el_pdos in self.pdos.iteritems():
            for site in el_pdos:
                for orb, pdos in site.iteritems():
                    if not tdos:
                        tdos = pdos
                    else:
                        tdos = add_dens(pdos, tdos)
        self.densities = tdos

    def get_shifted_dos(self, tol=0.2):
        """Get the DOS but shifted so the VBM is at 0.

        If there is no bandgap then the DOS is shifted by the fermi energy
        (which is equal to the VBM)

        Args:
            tol: The tolerance to use when finding the VBM

        Returns:
            VaspDos with all the energies shifted so that the VBM is at 0
        """
        cbm, vbm = self.get_cbm_vbm(tol=tol)
        efermi = self.efermi - vbm
        energies = self.energies - vbm
        return VaspDos(efermi, energies, self.densities, self.pdos,
                       self.structure)

    def get_smeared_vaspdos(self, sigma):
        """Smears the total DOS and all the partial DOSs

        Args:
            sigma: Std dev of Gaussian smearing function.

        Returns:
            Dict of Gaussian-smeared densities
            smear:
        """
        s_tdos = smear_dens(self.energies, self.densities, sigma)
        s_pdos = {}
        for el, el_pdos in self.pdos.iteritems():
            s_el_pdos = []
            for site in el_pdos:
                s_el_orbs = {}
                for orb, pdos in site.iteritems():
                    smeared_dens = smear_dens(self.energies, pdos, sigma)
                    s_el_orbs[orb] = smeared_dens
                s_el_pdos.append(s_el_orbs)
            s_pdos[el] = s_el_pdos
        return VaspDos(self.efermi, self.energies, s_tdos, s_pdos,
                       self.structure)

    def get_pdos(self, split_orbs=None, atoms=None, elements=None):
        """Gets the entire partial dos.

        Args:
            split_orbs: A list of the orbitals to split for each elment, in the
                form {Element: [orbs]}
            atoms: A list of atoms of over which to sum the DOS. The index
                starts at 1. If only the atom symbol is specified all the
                atoms are considered. Provided in the form {Element: [atoms]}
            elements: A dict of the elements containing which orbitals to plot

        Returns:
            The pdos in the form {Element: {Orbital: dos}}
        """
        pdoss = {}
        for el in self.pdos.keys():
            if (elements and el not in elements) or (atoms and el not in atoms):
                continue
            so = split_orbs[el] if (split_orbs and el in split_orbs) else None
            if atoms and el in atoms:
                sites = atoms[el] if (atoms[el] is not []) else\
                    range(self.structure[el])  # if empty list add all sites
            else:
                sites = None
            if elements and el in elements:
                orbs = elements[el]
            else:
                orbs = None
            pdoss[el] = self.get_element_dos(el, so, sites, orbs)
        return pdoss

    def get_element_dos(self, el, split_orbs=None, sites=None, orbs=None):
        """Get the projected DOS for an element.

        Args:
            el: Element from Poscar
            split_orbs: A list of the orbitals to split
            atoms: A list of atoms of over which to sum the DOS. The index
                starts at 1. If nothing is specified all the atoms are
                considered
            orbs: A list of orbitals to plot

        Returns:
            dict of {Orbital: dos}
        """
        el_dos = {}
        for i, site in enumerate(self.pdos[el]):
            if sites and i+1 not in sites:  # i+1 as VASP indices start at 1
                continue
            for orb, pdos in site.iteritems():
                if orbs and orb[0] not in orbs:
                    continue
                # Use the combined orbital if we don't want its split orbitals
                corb = orb if (split_orbs and orb[0] in split_orbs) else orb[0]
                if corb not in el_dos:
                    el_dos[corb] = pdos
                else:
                    el_dos[corb] = add_dens(el_dos[corb], pdos)
        return {orb: Dos(self.efermi, self.energies, densities)
                for orb, densities in el_dos.items()}

    def to_files(self, prefix="", split_orbs=None, atoms=None, elements=None):
        """Write the VASP DOS to a series of files containing the formatted data

        Args:
            prefix: Prefix for the file names
            split_elements: Write each element out as its own file
            split_orbitals: Which orbitals to split, in the form {El:[Orb]}
        """
        self.to_file(prefix)  # write tdos
        pdoss = self.get_pdos(split_orbs, atoms, elements)
        for el, el_pdos in pdoss.iteritems():
            data = []
            head = '#Energy'
            form = '%.3f'

            # here we add all the data requested into a big list
            # and generate the headers and formatting for when we want to write
            # the orbitals are sorted so that the DOS is always in the right
            # order
            for orb in sort_orbitals(el_pdos):
                pdos = el_pdos[orb]
                data.append(pdos.densities)
                orb_list = [orb]
                form_list = ['%G']
                head += "  " + fmt_head(self.spin, orb_list, '(down)')
                form += "  " + fmt_head(self.spin, form_list, '-', prefix=True)

            # used enumerate because we don't know how many pdos data sets
            # there will be in `data`. But we look over each energy and take the
            # densitiy of each dataset at that energy and add it to our string
            with open(prefix + "_" + el + '_dos.dat', 'w') as f:
                f.write(head + "\n")
                for i, energy in enumerate(self.energies):
                    dos_line = [energy]
                    for pdos in data:
                        for spin in range(self.spin):
                            dos_line.append(pdos[spin][i])
                    f.write((form + '\n') % tuple(dos_line))


def sort_orbitals(el_pdoss):
    """Get the sorted orbitals of a pdoss

    Sort the orbitals based in a standard format. E.g. s -> p -> d.
    Will also sort individual orbitals. This is useful for plotting/saving.

    Args:
        pdoss: An element pdoss in the form {Orbital: dens}

    Returns:
        A list of the sorted orbitals
    """
    sorted_orbitals = ['s', 'p', 'py', 'pz', 'px', 'd', 'dxy', 'dyz', 'dz2',
                       'dxz', 'dx2-y2', 'f', 'f-3', 'f-2', 'f-1', 'f-0', 'f1', 'f2',
                       'f3']
    sorted_keys = []
    unsorted_keys = el_pdoss.keys()
    for key in sorted_orbitals:
        if key in unsorted_keys:
            sorted_keys.append(key)
    return sorted_keys


def smear_dens(energies, dens, sigma):
    """Gaussian smear the densities.

    Returns the Dict represenation of the total densities, the
    {Spin: densities} but with a Gaussian smearing of std dev sigma applied
    about the Fermi level.
    Currently using Arons group's smearning method but could also try the
    one on pymatgen:

        from scipy.ndimage.filters import gaussian_filter1d
        smeared_dens = {}
        diff = [self.energies[i + 1] - self.energies[i]
                for i in range(len(self.energies) - 1)]
        avgdiff = sum(diff) / len(diff)
        for spin, dens in self.densities.items():
            smeared_dens[spin] = gaussian_filter1d(dens, sigma / avgdiff)

    Args:
        energies: The energies
        dens: The densities in the form {Spin: dens}
        sigma: Std dev of Gaussian smearing function.

    Returns:
        Dict of Gaussian-smeared densities
    """
    smeared_dens = {}
    kernx_p = np.arange(0.0, 10.0*sigma, energies[1]-energies[0])
    kernx = [-1.0 * v for v in kernx_p[::-1]]
    kernx = kernx + [v for v in kernx_p[1:]]
    kernx = np.array(kernx, dtype=np.float64)
    kern = np.exp(-1.0 * (kernx**2 / (2 * sigma**2)))
    kern = kern / kern.sum()
    for spin, dens in dens.iteritems():
        smeared_dens[spin] = np.convolve(dens, kern, mode='same')
    return smeared_dens


def add_dens(dens1, dens2):
    """Sum two sets of densities

    Doesn't check to make sure they have the same number of spins.

    Args:
        dens1: First density
        dens2: Second density

    Returns:
        Dict of {Spin: dens}
    """
    return {spin: np.array(dens1[spin]) + np.array(dens2[spin])
            for spin in dens1.keys()}


def fmt_head(spin, header, soc_addition, prefix=False):
    """Format the the headers for file writing.

    Formats a header based on whether or not it is for a SOC DOS calculation.
    E.g. With a header of [a, b], an addition of "-" and soc set to True, the
    formatted header will be:
        "a  -a  b  -b"

    Args:
        spin: 1 for a non-SOC header, 2 for a SOC header
        header: A list of items to go in the header
        soc_addition: The addition to add to the header items
        prefix: Whether or not the addition should go before or after the item.

    Returns:
        String containing the formatted header
    """
    fmt_head = []
    if spin == 2:
        for item in header:
            if prefix:
                fmt_head.append("%s  %s%s" % (item, soc_addition, item))
            else:
                fmt_head.append("%s  %s%s" % (item, item, soc_addition))
    else:
        fmt_head = header
    return "  ".join(fmt_head)
