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


def calculate_dielectric_properties(dielectric, properties,
                                    average=True):
    r"""Calculate optical properties from the dielectric function

    Supported properties:

    Absorption
    ~~~~~~~~~~

    The unit of alpha is :math:`\mathrm{cm}^{-1}`.

    Refractive index :math:`n` has real and imaginary parts:

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

        properties (set):
            The set of properties to return. Intermediate properties will be
            calculated as needed. Accepted values: 'eps_real', 'eps_im',
            'absorption', 'loss', 'n_real', 'n_imag'

        average (:obj:`bool`, optional): Average the dielectric response across
            the xx, yy, zz directions and calculate properties with scalar
            maths. Defaults to ``True``. If False, solve dielectric matrix to
            obtain directional properties, returning xx, yy, zz components.
            This may be significantly slower!

    Returns:
        dict of :obj:`tuple` of :obj:`list` of :obj:`float`: The optical
        absorption in :math:`\mathrm{cm}^{-1}`. If ``average`` is ``True``, the
        data will be returned as::

            ([energies], [property]).

        If ``average`` is ``False``, the data will be returned as::

            ([energies], [property_xx, property_yy, property_zz]).

        In both cases these are collected in a results dictionary with keys
        corresponding to the selected properties, e.g.::

            {'absorption': ([energies], [absorption]),
             'eps_real': ([energies], [eps_real])}
    """

    results = {}

    def _update_results(keys_vals):
        """Update results dict with selected properties only"""
        results.update({prop: (energies, data)
                        for prop, data in keys_vals.items()
                        if (prop in properties)})
        return results

    energies = np.array(dielectric[0])
    real_eps = np.array(dielectric[1])
    imag_eps = np.array(dielectric[2])

    if average:
        real_eps = np.average(real_eps[:, :3], axis=1)
        imag_eps = np.average(imag_eps[:, :3], axis=1)

        results = _update_results({'eps_real': real_eps,
                                   'eps_imag': imag_eps})

        eps = real_eps + 1j * imag_eps

        if 'loss' in properties:
            loss = -np.imag(1/eps)
            _update_results({'loss': loss})

        if properties.intersection({'n_real', 'n_imag', 'absorption'}):
            n = np.sqrt(eps)
            _update_results({'n_real': n.real,
                             'n_imag': n.imag})

        if 'absorption' in properties:
            alpha = n.imag * energies * 4 * np.pi / 1.23984212E-4
            _update_results({'absorption': alpha})

    else:
        # Work with eps as complex numbers in 9-column 'flattened' matrix
        # First interpret 6-column data as symmetric matrix
        # Input form xx yy zz xy yz xz
        # Indices     0  1  2  3  4  5
        n_rows = real_eps.shape[0]
        eps = real_eps + 1j * imag_eps
        eps = np.array([eps[:, 0], eps[:, 3], eps[:, 5],
                        eps[:, 3], eps[:, 1], eps[:, 4],
                        eps[:, 5], eps[:, 4], eps[:, 2]]).T

        _update_results(
            {'eps_real': eps.real[:, [0, 4, 8]],
             'eps_imag': eps.imag[:, [0, 4, 8]]})
        # Invert epsilon to obtain energy-loss function
        if 'loss' in properties:
            def matrix_loss_func(eps_row):
                eps_matrix = eps_row.reshape(3, 3)
                return -np.linalg.inv(eps_matrix).imag.flatten()

            loss = np.array([matrix_loss_func(row) for row in eps])

            _update_results({'loss': loss[:, [0, 4, 8]]})

        if properties.intersection({'n_real', 'n_imag', 'absorption'}):
            def matrix_n(eps_row):
                eps_matrix = eps_row.reshape(3, 3)
                eigenvals, v = np.linalg.eig(eps_matrix)
                d = np.diag(eigenvals)
                n = v @ np.sqrt(d) @ np.linalg.inv(v)  # Py3.5 matrix mult @ =D
                return n.flatten()

            n = np.array([matrix_n(row) for row in eps])

            _update_results({'n_real': n.real[:, [0, 4, 8]],
                             'n_imag': n.imag[:, [0, 4, 8]]})

        if 'absorption' in properties:
            alpha = (n.imag * energies.reshape(n_rows, 1) *
                     4 * np.pi / 1.23984212E-4)
            _update_results({'absorption': alpha[:, [0, 4, 8]]})

    return results


def write_files(abs_data, basename='absorption', prefix=None, directory=None):
    """Write the absorption or loss spectra to a file.

    Note that this function expects to receive an iterable series of spectra.

    Args:
        abs_data (tuple): Series (either :obj:`list` or :obj:`tuple`) of
            optical absorption or loss spectra. Each spectrum should be
            formatted as a :obj:`tuple` of :obj:`list` of :obj:`float`. If the
            data has been averaged, each spectrum should be::

                ([energies], [alpha])

            Else, if the data has not been averaged, each spectrum should be::

                ([energies], [alpha_xx, alpha_yy, alpha_zz]).

        prefix (:obj:`str`, optional): Prefix for file names.
        directory (:obj:`str`, optional): The directory in which to save files.
    """

    for i, absorption in enumerate(abs_data):
        num_txt = '_{}'.format(i + 1) if len(abs_data) > 1 else ''
        prefix_txt = '{}_'.format(prefix) if prefix else ''
        filename = prefix_txt + basename + num_txt + '.dat'

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


def kkr(de, eps_imag, cshift=1e-6):
    """Kramers Kronig transformation of imaginary dielectric function

    Args:
        de (:obj:`float`): Energy difference between evenly-spaced energy
            values corresponding to dielectric data
        eps_imag (:obj:`list` or :obj:`np.array`): Evenly-spaced sequence of
            frequency-dependent 3x3 dielectric matrices (imaginary component
            only)
        cshift (:obj:`float`, optional): imaginary finite shift used in
            integration; this should be small (and results should not be very
            sensitive)

        returns:
            :obj:`numpy.array`
                Real part of frequency-dependent dielectric function
                corresponding to eps_imag. Array shape (NEDOS, 3, 3)
    """
    eps_imag = np.array(eps_imag)
    nedos = eps_imag.shape[0]
    cshift = complex(0, cshift)
    w_i = np.arange(0, (nedos - 0.5)*de, de, dtype=np.complex_)
    w_i = np.reshape(w_i, (nedos, 1, 1))

    def integration_element(w_r):
        factor = w_i / (w_i**2 - w_r**2 + cshift)
        total = np.sum(eps_imag * factor, axis=0)
        return total * (2/np.pi) * de + np.diag([1, 1, 1])

    return np.real([integration_element(w_r) for w_r in w_i[:, 0, 0]])

