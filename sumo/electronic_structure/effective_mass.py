# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Module containing helper functions for calculating band effective masses.
"""

import numpy as np

from scipy.optimize import curve_fit
from scipy.constants import physical_constants

eV_to_hartree = physical_constants['electron volt-hartree relationship'][0]
bohr_to_m = physical_constants['Bohr radius'][0]
angstrom_to_bohr = bohr_to_m / 1e-10


def get_fitting_data(bs, spin, band_id, kpoint_id, num_sample_points=3):
    """Extract fitting data for band extrema based on spin, kpoint and band.

    Searches forward and backward from the extrema point, but will only sample
    there data if there are enough points in that direction.

    Args:
        bs (:obj:`~pymatgen.electronic_structure.bandstructure.BandStructureSymmLine`):
            The band structure.
        spin (:obj:`~pymatgen.electronic_structure.core.Spin`): Which spin
            channel to sample.
        band_id (int): Index of the band to sample.
        kpoint_id (int): Index of the kpoint to sample.

    Returns:
        list: The data necessary to calculate the effective mass, along with
        some metadata. Formatted as a :obj:`list` of :obj:`dict`, each with the
        keys:

        'energies' (:obj:`numpy.ndarray`)
            Band eigenvalues in eV.

        'distances' (:obj:`numpy.ndarray`)
            Distances of the k-points in reciprocal space.

        'band_id' (:obj:`int`)
            The index of the band,

        'spin' (:obj:`~pymatgen.electronic_structure.core.Spin`)
            The spin channel

        'start_kpoint' (:obj:`int`)
            The index of the k-point at which the band extrema occurs

        'end_kpoint' (:obj:`int`)
            The k-point towards which the data has been sampled.
    """
    # branch data provides data about the start and end points
    # of specific band paths
    branch_data = [b for b in bs.get_branch(kpoint_id)
                   if b['index'] == kpoint_id][0]
    start_kpoint = bs.kpoints[kpoint_id]

    fitting_data = []

    # check to see if there are enough points to sample from first
    # check in the forward direction
    if kpoint_id + num_sample_points <= branch_data['end_index']:

        # calculate sampling limits
        start_id = kpoint_id
        end_id = kpoint_id + num_sample_points + 1

        energies = np.array(bs.bands[spin][band_id][start_id:end_id].copy())
        dists = np.array(bs.distance[start_id:end_id].copy())

        # normalise eigenvalues and distances to starting point
        energies -= bs.bands[spin][band_id][kpoint_id]
        dists -= bs.distance[kpoint_id]

        # symmetrise the data to make fitting more reliable
        energies = np.concatenate([energies[::-1], energies[1:]])
        dists = np.concatenate([-dists[::-1], dists[1:]])

        end_kpoint = bs.kpoints[branch_data['end_index']]
        data = {'energies': energies, 'distances': dists, 'band_id': band_id,
                'spin': spin, 'start_kpoint': start_kpoint,
                'end_kpoint': end_kpoint}
        fitting_data.append(data)

    # check in the backward direction
    if kpoint_id - num_sample_points >= branch_data['start_index']:

        # calculate sampling limits
        start_id = kpoint_id - num_sample_points
        end_id = kpoint_id + 1

        energies = bs.bands[spin][band_id][start_id:end_id].copy()
        dists = bs.distance[start_id:end_id].copy()

        # normalise eigenvalues and distances to starting point
        energies -= bs.bands[spin][band_id][kpoint_id]
        dists -= bs.distance[kpoint_id]

        # symmetrise the data to make fitting more reliable
        energies = np.concatenate([energies[:-1], energies[::-1]])
        dists = np.concatenate([dists[:-1], -dists[::-1]])

        end_kpoint = bs.kpoints[branch_data['start_index']]
        data = {'energies': energies, 'distances': dists, 'band_id': band_id,
                'spin': spin, 'start_kpoint': start_kpoint,
                'end_kpoint': end_kpoint}
        fitting_data.append(data)

    return fitting_data


def fit_effective_mass(distances, energies, parabolic=True):
    """Fit the effective masses using either a parabolic or nonparabolic fit.

    Args:
        distances (:obj:`numpy.ndarray`): The x-distances between k-points in
            reciprocal Angstroms, normalised to the band extrema.
        energies (:obj:`numpy.ndarray`): The band eigenvalues normalised to the
            eigenvalue of the band extrema.
        parabolic (:obj:`bool`, optional): Use a parabolic fit of the band
            edges. If ``False`` then nonparabolic fitting will be attempted.
            Defaults to ``True``.

    Returns:
        float: The effective mass in units of electron rest mass, :math:`m_0`.
    """
    if parabolic:
        fit = np.polyfit(distances, energies, 2)
        c = 2 * fit[0]  # curvature therefore 2 * the exponent on the ^2 term

    else:
        # Use non parabolic description of the bands
        def f(x, alpha, d):
            top = np.sqrt(4 * alpha * d * x**2 + 1) - 1
            bot = 2 * alpha
            return top / bot

        # set boundaries for curve fitting: alpha > 1e-8
        # as alpha = 0 causes an error
        bounds = ((1e-8, -np.inf), (np.inf, np.inf))
        popt, _ = curve_fit(f, distances, energies, p0=[1., 1.],
                            bounds=bounds)
        c = 2 * popt[1]

    # coefficient is currently in eV/Angstrom^2/h_bar^2
    # want it in atomic units so Hartree/bohr^2/h_bar^2
    eff_mass = (angstrom_to_bohr**2 / eV_to_hartree) / c
    return eff_mass
