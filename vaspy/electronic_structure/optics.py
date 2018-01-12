# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

import numpy as np

from scipy.ndimage.filters import gaussian_filter1d


def broaden_eps(dielectric, sigma):
    """Apply Gaussian broadening to the dielectric response.

    Args:
        dielectric_data (tuple): The dielectric data formatted
            as given by pymatgen.io.vasp.Vasprun.dielectric. E.g. a tuple
            of 3 values containing the energy, the real part tensor, and
            the imaginary part tensor ([energies],[[real_partxx,real_partyy,
            real_partzz,real_partxy,real_partyz,real_partxz]],[[imag_partxx,
            imag_partyy,imag_partzz,imag_partxy, imag_partyz, imag_partxz]])
        sigma (float): Standard deviation for Gaussian kernel

    Returns:
        The broadened dielectric response, formatted as a tuple
        of 3 values containing the energy, the real part tensor, and
        the imaginary part tensor ([energies],[[real_partxx,real_partyy,
        real_partzz,real_partxy,real_partyz,real_partxz]],[[imag_partxx,
        imag_partyy,imag_partzz,imag_partxy, imag_partyz, imag_partxz]])
    """
    e = dielectric[0]
    diff = [e[i + 1] - e[i] for i in range(len(e) - 1)]
    diff_avg = sum(diff) / len(diff)
    real = [gaussian_filter1d(np.array(dielectric[1])[:, x], sigma / diff_avg)
            for x in range(6)]
    imag = [gaussian_filter1d(np.array(dielectric[2])[:, x], sigma / diff_avg)
            for x in range(6)]

    return (e, np.array(real).T, np.array(imag).T)


def calculate_alpha(dielectric, average=True):
    """Calculate the absorption coefficient based on the dielectric data.
    The unit of alpha is cm^-1.

    Based on http://course.ee.ust.hk/elec342/notes/Lecture%207_absorption%20and%20dispersion.pdf

    Refactive index n has real and imaginary parts:

        n = [(e' + ie"/e_0]^(1/2) = n' + in"

    Relationship between a and imaginary n":

        a = 4*pi*n"/wavelength

    Where wavelength is = hc/E

    Args:
        dielectric_data (tuple): The dielectric data formatted
            as given by pymatgen.io.vasp.Vasprun.dielectric. E.g. a tuple
            of 3 values containing the energy, the real part tensor, and
            the imaginary part tensor ([energies],[[real_partxx,real_partyy,
            real_partzz,real_partxy,real_partyz,real_partxz]],[[imag_partxx,
            imag_partyy,imag_partzz,imag_partxy, imag_partyz, imag_partxz]])
        average (bool, optional): Average the dielectric response across
            all lattice directions.

    Returns:
        The optical absorption, formatted as a tuple of ([energies], [alpha]).
        If average is set to false, the data is returned as ([energies], [alphaxx,
        alphayy, alphazz]).
    """
    real_eps = np.array(dielectric[1])[:, :3]
    imag_eps = np.array(dielectric[2])[:, :3]
    energies = np.array(dielectric[0])

    if average:
        real_eps = np.average(real_eps, axis=1)
        imag_eps = np.average(imag_eps, axis=1)

    eps = real_eps + 1j * imag_eps
    imag_ref_index = np.sqrt(eps).imag

    if average:
        alpha = imag_ref_index * energies * 4 * np.pi / 1.23984212E-4
    else:
        alpha = imag_ref_index * energies[:, None] * 4 * np.pi / 1.23984212E-4

    return (energies, alpha)


def write_files(abs_data, prefix=None, directory=None):
    """Write the absorption spectra to a (series of) file(s).

    Args:
        abs_data (list): A list of absorption spectra. Each spectra should
            be formatted as a tuple of ([energies], [alpha]) or if the data
            has not been averaged: ([energies], [alphaxx, alphayy, alphazz]).
        prefix (str): Prefix for file names.
        directory (str): The directory in which to save files.
    """
    for i, absorption in enumerate(abs_data):
        num = '_{}'.format(i) if len(abs_data) > 1 else ''
        basename = 'absorption{}.dat'.format(num)
        filename = '{}_{}'.format(prefix, basename) if prefix else basename
        if directory:
            filename = os.path.join(directory, filename)

        header = 'energy(eV)'
        if len(absorption[1].shape) == 2:
            header += ' alpha_xx alpha_yy alpha_zz'
            data = np.concatenate((absorption[0][:, None], absorption[1]), axis=1)
        else:
            header += ' alpha'
            data = np.stack((absorption[0], absorption[1]), axis=1)

        np.savetxt(filename, data, header=header)
