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
        dielectric (tuple): The high-frequency dielectric data, following
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
    real = [
        gaussian_filter1d(np.array(dielectric[1])[:, x], sigma / diff_avg)
        for x in range(6)
    ]
    imag = [
        gaussian_filter1d(np.array(dielectric[2])[:, x], sigma / diff_avg)
        for x in range(6)
    ]

    return e, np.array(real).T, np.array(imag).T


def calculate_dielectric_properties(dielectric, properties, mode="average"):
    r"""Calculate optical properties from the dielectric function

    Supported properties:

    *Absorption*

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
        dielectric (tuple): The high-frequency dielectric data, following
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
            calculated as needed. Accepted values: 'eps_real', 'eps_imag',
            'absorption', 'loss', 'n_real', 'n_imag'

        mode (str): The format of the results: Options are:

            - "average": Average the dielectric response for all directions.
            - "trace": The trace of the dielectric tensor (i.e., along xx, yy, zz).
            - "full": The full dielectric tensor (i.e., xx, xy, xz, yx, yy, yz, zx, zy,
              zz).
            - "eigs": Calculate the eigenvalues of the dielectric response. The
              eigenvalues are sorted from smallest to largest.

    Returns:
        tuple of ``(energies, {property_name: property_value})``: The optical absorption is given in
        :math:`\mathrm{cm}^{-1}`. If ``mode`` is ``"average"``, the
        property_value will be returned as::

            [property]

        If ``mode`` is ``"trace"``, the data will be returned as::

            [property_xx, property_yy, property_zz]

        If ``mode`` is ``"full"``, the data will be returned as::

            [xx, xy, xz, yx, yy, yz, zx, zy, zz]

        If ``mode`` is ``"eigs"``, the data will be returned as::

            [eig_1, eig_2, eig_3]

        In all cases these are collected in a results dictionary with keys
        corresponding to the selected properties, e.g.::

            {'absorption': [absorption], 'eps_real': [eps_real]}
    """
    energies = np.array(dielectric[0])

    # Work with eps as complex numbers in Nx3x3 matrix
    # First interpret 6-column data as symmetric matrix
    # Input form xx yy zz xy yz xz
    # Indices     0  1  2  3  4  5
    real_eps = np.array(dielectric[1])[:, [[0, 3, 5], [3, 1, 4], [5, 4, 2]]]
    imag_eps = np.array(dielectric[2])[:, [[0, 3, 5], [3, 1, 4], [5, 4, 2]]]
    eps_full = real_eps + 1j * imag_eps

    # take sqrt of eps matrix; if eps = V S V^-1; then eps^1/2 = V S^{1/2} V^-1;
    eigvals, eigvecs = np.linalg.eig(eps_full)

    # fancy einsum to calculate V S^{1/2} V^-1 at every energy
    n = np.einsum("ijk,ik,ikl->ijl", eigvecs, np.sqrt(eigvals), np.linalg.inv(eigvecs))

    # calculate optical absorption
    alpha = n.imag * energies[:, None, None] * 4 * np.pi / 1.23984212e-4

    # Invert epsilon to obtain energy-loss function
    loss = -np.linalg.inv(eps_full).imag

    if mode == "average":
        eps = np.linalg.eigvals(eps_full).mean(axis=1)
        n = np.linalg.eigvals(n).mean(axis=1)
        loss = np.linalg.eigvalsh(loss).mean(axis=1)
        alpha = np.linalg.eigvalsh(alpha).mean(axis=1)
    elif mode == "eigs":
        eps = eigvals
        n = np.linalg.eigvals(n)
        loss = np.linalg.eigvalsh(loss)
        alpha = np.linalg.eigvalsh(alpha)
    elif mode == "full":
        eps = eps_full.reshape(-1, 9)
        n = n.reshape(-1, 9)
        loss = loss.reshape(-1, 9)
        alpha = alpha.reshape(-1, 9)
    elif mode == "trace":
        eps = np.diagonal(eps_full.T)
        n = np.diagonal(n.T)
        loss = np.diagonal(loss.T)
        alpha = np.diagonal(alpha.T)
    else:
        raise ValueError(
            f"Unsupported mode: {mode}. Options are 'average', 'eigs', 'full', 'trace'"
        )

    data = {
        "eps_real": eps.real,
        "eps_imag": eps.imag,
        "n_real": n.real,
        "n_imag": n.imag,
        "loss": loss,
        "absorption": alpha,
    }

    return energies, {k: data[k] for k in properties}


def write_files(abs_data, basename="absorption", prefix=None, directory=None):
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
        num_txt = f"_{i + 1}" if len(abs_data) > 1 else ""
        prefix_txt = f"{prefix}_" if prefix else ""
        filename = prefix_txt + basename + num_txt + ".dat"

        if directory:
            filename = os.path.join(directory, filename)

        header = "energy(eV)"
        if len(absorption[1].shape) == 2:
            header += " alpha_xx alpha_yy alpha_zz"
            data = np.concatenate((absorption[0][:, None], absorption[1]), axis=1)
        else:
            header += " alpha"
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

    Returns:
        (:obj:`numpy.array`) Real part of frequency-dependent dielectric function
        corresponding to eps_imag. Array shape (NEDOS, 3, 3)
    """
    eps_imag = np.array(eps_imag)
    nedos = eps_imag.shape[0]
    cshift = complex(0, cshift)
    w_i = np.arange(0, (nedos - 0.5) * de, de, dtype=np.complex_)
    w_i = np.reshape(w_i, (nedos, 1, 1))

    def integration_element(w_r):
        factor = w_i / (w_i**2 - w_r**2 + cshift)
        total = np.sum(eps_imag * factor, axis=0)
        return total * (2 / np.pi) * de + np.diag([1, 1, 1])

    return np.real([integration_element(w_r) for w_r in w_i[:, 0, 0]])
