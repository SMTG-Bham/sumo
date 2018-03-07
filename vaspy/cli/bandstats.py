# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

import os
import sys
import logging
import glob
import argparse

import numpy as np
from scipy.optimize import curve_fit

from vaspy.electronic_structure.bandstructure import \
    get_reconstructed_band_structure

from pymatgen.io.vasp.outputs import BSVasprun
from pymatgen.electronic_structure.core import Spin

"""
A script to calculate effective masses and band degeneracies
"""

__author__ = "Alex Ganose"
__version__ = "1.0"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "March 4, 2018"


kpt_str = '[{k[0]:.2f}, {k[1]:.2f}, {k[2]:.2f}]'

# TODO:
#  - Would be good to reimplement get_vbm and get_cbm methods, with the ability to
#    give degeneracies + other band edges within a certain tolerance.
#  - Sample custom k-point, band, spin combinations


def bandstats(filenames=None, num_sample_points=3, temperature=None,
              degeneracy_tol=1e-4, parabolic=True):
    """Calculate the effective masses of the bands of a semiconductor.

    Args:
        filenames (list, None): List of vasprun.xml files from which to extract
            the band structure. If None, the script will look for folders named
            split-*/ and the current directory for vaspruns.
        num_sample_points (int): Number of points to sample when fitting the
            effective masses.
        temperature (int, None): If not None, the script will attempt to find
            band edges within kB * T of VBM and CBM.
        degeneracy_tol (float): Tolerance for determining the degeneracy of the
            VBM/CBM.
        parabolic (bool): Use a parabolic fit of the band edges. If False then
            nonparabolic fitting will be attempted.

    Returns:
        A dictionary of the hole and electron effective masses, with keys:
        'hole_data' and 'electron_data'. The data is a list of dictionaries
        with the format:
            {'effective_mass': The effective mass in m_0.
             'energies': np.array of band eigenvalues in eV used for fitting,
             'distances': np.array of distances in reciprocal space used for
                 fitting,
             'band_id': The index of the sampled band,
             'spin': The spin direction of the sampled band,
             'start_kpoint': The k-point of the band extrema
             'end_kpoint': The k-point towards which the data has been
                 sampled.}
    """
    if not filenames:
        folders = glob.glob('split-*')
        folders = sorted(folders) if folders else ['.']
        filenames = []
        for fol in folders:
            vr_file = os.path.join(fol, 'vasprun.xml')
            if os.path.exists(vr_file):
                filenames.append(vr_file)
            else:
                logging.error('ERROR: No vasprun.xml found in {}!'.format(fol))
                sys.exit()

    bandstructures = []
    for vr_file in filenames:
        vr = BSVasprun(vr_file, parse_projected_eigen=False)
        bs = vr.get_band_structure(line_mode=True)
        bandstructures.append(bs)
    bs = get_reconstructed_band_structure(bandstructures)

    if bs.is_metal():
        logging.error('ERROR: System is metallic!')
        sys.exit()

    log_band_gap_information(bs)

    vbm_data = bs.get_vbm()
    cbm_data = bs.get_cbm()

    logging.info('\nValence band maximum:')
    log_band_edge_information(bs, vbm_data)

    logging.info('\nConduction band minimum:')
    log_band_edge_information(bs, cbm_data)

    if parabolic:
        logging.info('\nUsing parabolic fitting of the band edges')
    else:
        logging.info('\nUsing nonparabolic fitting of the band edges')

    if temperature:
        logging.error('ERROR: This feature is not yet supported!')

    else:
        # Work out where the hole and electron band edges are.
        # Fortunately, pymatgen does this for us. Points at which to calculate
        # the effective mass are identified as a tuple of:
        # (spin, band_index, kpoint_index)
        hole_extrema = []
        for spin, bands in vbm_data['band_index'].items():
            hole_extrema.extend([(spin, band, kpoint) for band in bands
                                 for kpoint in vbm_data['kpoint_index']])

        elec_extrema = []
        for spin, bands in cbm_data['band_index'].items():
            elec_extrema.extend([(spin, band, kpoint) for band in bands
                                 for kpoint in cbm_data['kpoint_index']])

        # extract the data we need for fitting from the band structure
        hole_data = []
        for extrema in hole_extrema:
            hole_data.extend(get_fitting_data(bs, *extrema,
                             num_sample_points=num_sample_points))

        elec_data = []
        for extrema in elec_extrema:
            elec_data.extend(get_fitting_data(bs, *extrema,
                             num_sample_points=num_sample_points))

    # calculate the effective masses and log the information
    logging.info('\nHole effective masses:')
    for data in hole_data:
        eff_mass = fit_effective_mass(data['distances'], data['energies'],
                                      parabolic=parabolic)
        data['effective_mass'] = eff_mass
        log_effective_mass_data(data, bs.is_spin_polarized)

    logging.info('\nElectron effective masses:')
    for data in elec_data:
        eff_mass = fit_effective_mass(data['distances'], data['energies'],
                                     parabolic=parabolic)
        data['effective_mass'] = eff_mass
        log_effective_mass_data(data, bs.is_spin_polarized)

    return {'hole_data': hole_data, 'electron_data': elec_data}


def get_fitting_data(bs, spin, band_id, kpoint_id, num_sample_points=3):
    """Extract fitting data for band extrema based on spin, kpoint and band.

    Searches forward and backward from the extrema point, to check there is
    enough data to sample.

    Args:
        bs (BandStructureSymmLine): Band structure object from which to extract
            fitting data.
        spin (Spin): Spin of the band to sample.
        band_id (int): Index of the band to sample.
        kpoint_id (int): Index of the kpoint to sample.

    Returns:
        A list of dicts containing the data necessary to calculate the
        effective mass, along with some metadata. Formatted as:
            {'energies': np.array of band eigenvalues in eV,
             'distances': np.array of distances in reciprocal space,
             'band_id': The index of the sampled band,
             'spin': The spin direction of the sampled band,
             'start_kpoint': The k-point of the band extrema
             'end_kpoint': The k-point towards which the data has been
             sampled.}
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
        distances (np.array): The x-distances between k-points in reciprocal
            Angstroms, normalised to the band extrema.
        energies (np.array): The band eigenvalues normalised to the eigenvalue
            of the band extrema.
        parabolic (bool): Use a parabolic fit. If False then nonparabolic
            fitting will be attempted.

    Returns:
        The effective mass in units of electron rest mass m_0.
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
        popt, pconv = curve_fit(f, distances, energies, p0=[1., 1.],
                                bounds=bounds)
        c = 2 * popt[1]

    # coefficient is currently in eV/Angstrom^2/h_bar^2
    # want it in units of m_e
    eff_mass = 7.61996348863/c
    return eff_mass


def log_band_gap_information(bs):
    """Log data about the direct and indirect band gaps.

    Args:
        bs (BandStructureSymmLine): The band structure object.
    """
    bg_data = bs.get_band_gap()
    if not bg_data['direct']:
        logging.info('Indirect band gap: {:.3f} eV'.format(bg_data['energy']))

    direct_data = bs.get_direct_band_gap_dict()

    if bs.is_spin_polarized:
        spins = direct_data.keys()
        direct_bg = direct_data[spins[0]]['value']
        logging.info('Direct band gap: {:.3f} eV'.format(direct_bg))

        for spin in direct_bg.keys():
            direct_kindex = direct_data[spin]['kpoint_index']
            direct_kpoint = bs.kpoints[direct_kindex].frac_coords
            direct_kpoint = kpt_str.format(k=direct_kpoint)
            eq_kpoints = bs.get_equivalent_kpoints(direct_kindex)
            k_indexes = ', '.join(map(str, eq_kpoints))
            b_indexes = ', '.join(map(str, direct_data[spin]['band_indices']))

            logging.info('  {}:'.format(spin.name.capitalize()))
            logging.info('    k-point: {}'.format(direct_kpoint))
            logging.info('    k-point indexes: {}'.format(k_indexes))
            logging.info('    Band indexes: {}'.format(b_indexes))

    else:
        direct_bg = direct_data[Spin.up]['value']
        logging.info('Direct band gap: {:.3f} eV'.format(direct_bg))

        direct_kindex = direct_data[Spin.up]['kpoint_index']
        direct_kpoint = kpt_str.format(k=bs.kpoints[direct_kindex].frac_coords)
        k_indexes = ', '.join(map(str, bs.get_equivalent_kpoints(direct_kindex)))
        b_indexes = ', '.join(map(str, direct_data[Spin.up]['band_indices']))

        logging.info('  k-point: {}'.format(direct_kpoint))
        logging.info('  k-point indexes: {}'.format(k_indexes))
        logging.info('  Band indexes: {}'.format(b_indexes))


def log_band_edge_information(bs, edge_data):
    """Log data about the valence band maximum or conduction band minimum.

    Args:
        bs (BandStructureSymmLine): The band structure object.
        edge_data (dict): The dict from bs.get_vbm() or bs.get_cbm()
    """
    if bs.is_spin_polarized:
        spins = edge_data['band_index'].keys()
        b_indexes = [', '.join(map(str, edge_data['band_index'][spin]))
                     + '({})'.format(spin.name.capitalise()) for spin in spins]
        b_indexes = ', '.join(b_indexes)
    else:
        b_indexes = ', '.join(map(str, edge_data['band_index'][Spin.up]))

    kpoint = edge_data['kpoint']
    kpoint_str = kpt_str.format(k=kpoint.frac_coords)
    k_indexes = ', '.join(map(str, edge_data['kpoint_index']))

    if kpoint.label:
        k_loc = kpoint.label
    else:
        branch = bs.get_branch(edge_data['kpoint_index'][0])[0]
        k_loc = 'between {}'.format(branch['name'])

    logging.info('  Energy: {:.3f} eV'.format(edge_data['energy']))
    logging.info('  k-point: {}'.format(kpoint_str))
    logging.info('  k-point location: {}'.format(k_loc))
    logging.info('  k-point indexes: {}'.format(k_indexes))
    logging.info('  Band indexes: {}'.format(b_indexes))


def log_effective_mass_data(data, is_spin_polarized):
    """Log data about the effective masses and their directions.

    Args:
        data (dict): A dictionary containing the data about the effective
            masses. This is mainly the data from get_fitting_data() but
            with the addition of the 'effective_mass' key. Formatted as:
                {'effective_mass': The effective mass in m_0.
                 'band_id': The index of the sampled band,
                 'spin': The spin direction of the sampled band,
                 'start_kpoint': The k-point of the band extrema
                 'end_kpoint': The k-point towards which the data has
                    been sampled.}
        is_spin_polarized (bool): Whether the system is spin polarized.
    """
    s = ' ({})'.format(data['spin'].name) if is_spin_polarized else ''
    band_str = 'band {}{}'.format(data['band_id'], s)

    start_kpoint = data['start_kpoint']
    end_kpoint = data['end_kpoint']
    eff_mass = data['effective_mass']

    kpoint_str = kpt_str.format(k=start_kpoint.frac_coords)
    if start_kpoint.label:
        kpoint_str += ' ({})'.format(start_kpoint.label)
    kpoint_str += ' -> '
    kpoint_str += kpt_str.format(k=end_kpoint.frac_coords)
    if end_kpoint.label:
        kpoint_str += ' ({})'.format(end_kpoint.label)

    logging.info('  m_e: {:.3f} | {} | {}'.format(eff_mass, band_str,
                                                  kpoint_str))

def main():
    parser = argparse.ArgumentParser(description="""
    bandstats is provides information on the band gap and effective
    masses of semiconductors.""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('-f', '--filenames', default=None, nargs='+',
                        help="one or more vasprun.xml files to plot")
    parser.add_argument('-n', '--nonparabolic', default=True,
                        action='store_false',
                        help="""Use a nonparabolic model to fit the
                                effective masses""")
    parser.add_argument('-s', '--sample-points', default=3, type=int,
                        dest='sample_points',
                        help="number of k-points to sample in fitting")
    parser.add_argument('-t', '--temperature', default=None, type=float,
                        help="Find band edges within kB * T of VBM and CBM.")

    args = parser.parse_args()
    logging.basicConfig(filename='vaspy-bandstats.log', level=logging.DEBUG,
                        filemode='w', format='%(message)s')
    console = logging.StreamHandler()
    logging.info(" ".join(sys.argv[:]))
    logging.getLogger('').addHandler(console)

    bandstats(filenames=args.filenames, num_sample_points=args.sample_points,
              temperature=args.temperature, parabolic=args.nonparabolic)

if __name__ == "__main__":
    main()
