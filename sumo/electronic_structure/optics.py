# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Module containing functions to process dielectric and optical absorption data.

TODO:
    * Remove magic values
"""

import os
import numpy as np

from scipy.ndimage.filters import gaussian_filter1d


def broaden_eps(dielectric, sigma):
    """Apply gaussian broadening to the dielectric response.

    Args:
        dielectric_data (tuple): The high-frequency dielectric data, following
            the same format as
            :attr:`pymatgen.io.vasp.outputs.Vasprun.dielectric`.
            This is a :obj:`tuple` containing the energy, the real part of the
            dielectric tensor, and the imaginary part of the tensor, as a
            :obj:`list` of :obj:`floats`. E.g.::

                (
                    [energies],
                    [[real_xx, real_yy, real_zz, real_xy, real_yz, real_xz]],
                    [[imag_xx, imag_yy, imag_zz, imag_xy, imag_yz, imag_xz]]
                )

        sigma (float): Standard deviation for gaussian broadening.

    Returns:
        :obj:`tuple` of :obj:`list` of :obj:`list` of :obj:`float`: The
        broadened dielectric response. Returned as a tuple containing the
        energy, the real part of the dielectric tensor, and the imaginary
        part of the tensor. E.g.::

            (
                [energies],
                [[real_xx, real_yy, real_zz, real_xy, real_yz, real_xz]],
                [[imag_xx, imag_yy, imag_zz, imag_xy, imag_yz, imag_xz]]
            )
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
    r"""Calculate the optical absorption from the high-frequency dielectric.

    The unit of alpha is :math:`\mathrm{cm}^{-1}`.

    Refactive index :math:`n` has real and imaginary parts:

    .. math::

        n = [(e^\prime + ie^{\prime\prime} / e_0]^{1/2}
          = n^\prime + in^{\prime\prime}

    Relationship between :math:`a` and imaginary :math:`n^{\prime\prime}`:

    .. math::

        a = 4 \pi n^{\prime\prime} / \lambda

    Where:

    .. math:: \lambda = hc/E

    Args:
        dielectric_data (tuple): The high-frequency dielectric data, following
            the same format as :obj:`pymatgen.io.vasp.Vasprun.dielectric`.
            This is a :obj:`tuple` containing the energy, the real part of the
            dielectric tensor, and the imaginary part of the tensor, as a
            :obj:`list` of :obj:`floats`. E.g.::

                (
                    [energies],
                    [[real_xx, real_yy, real_zz, real_xy, real_yz, real_xz]],
                    [[imag_xx, imag_yy, imag_zz, imag_xy, imag_yz, imag_xz]]
                )

        average (:obj:`bool`, optional): Average the dielectric response across
            all lattice directions. Defaults to ``True``.

    Returns:
        :obj:`tuple` of :obj:`list` of :obj:`float`: The optical absorption in
        :math:`\mathrm{cm}^{-1}`. If ``average`` is ``True``, the data will be
        returned as::

            ([energies], [alpha]).

        If ``average`` is ``False``, the data will be returned as::

            ([energies], [alpha_xx, alpha_yy, alpha_zz]).
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
    """Write the absorption spectra to a file.

    Note that this function expects to recieve an iterable series of spectra.

    Args:
        abs_data (tuple): Series (either :obj:`list` or :obj:`tuple`) of
            optical absorption spectra. Each spectra should be formatted as a
            :obj:`tuple` of :obj:`list` of :obj:`float`. If the absorption data
            has been averaged, each spectra should be::

                ([energies], [alpha])

            Else, if the data has not been averaged, each spectra should be::

                ([energies], [alpha_xx, alpha_yy, alpha_zz]).

        prefix (:obj:`str`, optional): Prefix for file names.
        directory (:obj:`str`, optional): The directory in which to save files.
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
            data = np.concatenate((absorption[0][:, None], absorption[1]),
                                  axis=1)
        else:
            header += ' alpha'
            data = np.stack((absorption[0], absorption[1]), axis=1)

        np.savetxt(filename, data, header=header)
