#!/usr/bin/env python

"""
A script to extra effective masses and band gaps from the a bandstructure
"""

__author__ = "Alex Ganose"
__version__ = "0.1"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "Oct 9, 2015"

import logging
import argparse

from itertools import izip, tee

from vaspy.vasp_output import Procar
from vaspy.vasp_input import Poscar

import numpy as np
from numpy import sqrt, diag, pi, polyfit, array, dot, append
from numpy.linalg import inv, norm, eig, det

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
import scipy.constants

const = scipy.constants.physical_constants


def bandstats(poscar, procar, steps=2, tol=0.0002, labels=None, vkpts=None,
              ckpts=None):
    ang_to_bohr = const['Angstrom star'][0] / const['Bohr radius'][0]
    #recip_lattice = inv(poscar.lattice).T * 2e0 * pi
    recip_lattice = inv(np.array(poscar.lattice) * ang_to_bohr).T * 2 * pi
    vb = array(procar.eigenvals[procar.vbm['band']-1])
    cb = array(procar.eigenvals[procar.cbm['band']-1])

    kpoints = array(procar.kpoints)
    logging.info("")
    vbms, cbms, eg = get_band_edges(kpoints, vb, cb, procar.vbm, procar.cbm,
                                    tol, vkpts, ckpts)
    vb_data = get_data(kpoints, vb, vbms, steps, labels)
    cb_data = get_data(kpoints, cb, cbms, steps, labels)

    for vb_dir in vb_data:
        em = calc_em(vb_dir, recip_lattice, eg)
        if vb_dir['label'] != "":
            label = vb_dir['label']
        else:
            label = "kpoints %d to %d" % (vb_dir['start']+1, vb_dir['end']+1)

        logging.info("m_h from %s = %.4f" % (label, em))

    logging.info("")
    for cb_dir in cb_data:
        em = calc_em(cb_dir, recip_lattice, eg)
        if cb_dir['label'] != "":
            label = cb_dir['label']
        else:
            label = "kpoints %d to %d" % (cb_dir['start']+1, cb_dir['end']+1)
        logging.info("m_e from %s = %.4f" % (label, em))


def calc_em(data, recip_lattice, eg):
    eV_to_hartree = const['electron volt-hartree relationship'][0]
    x = [norm(dot(kpt, recip_lattice)) for kpt in data['kpoints']]
    e = data['vals'] * eV_to_hartree

    norm_x = np.append(-np.flipud(x), x)
    norm_e = np.append(np.flipud(e), e)

    #norm_x = x
    #norm_e = e
    pfit = polyfit(x, e, 2)
    co = pfit[0] * 2  # curvature therefore 2 * the exponent on the ^2 term

    def f(x, alpha, d):
        top = sqrt(4 * alpha * d * x**2 + 1) - 1
        bot = 2 * alpha
        return top/bot
    popt, pconv = curve_fit(f, norm_x, norm_e, p0=[0.5, 0.5])
    fxs = [f(x_i, popt[0], popt[1]) for x_i in norm_x]

    if False:
        plt.figure(figsize=(6, 6))
        plt.plot(norm_x, norm_e)
        plt.plot(norm_x, fxs)
        plt.plot(norm_x, fxs)
        pfit = polyfit(norm_x, norm_e, 2)
        poly = np.poly1d(pfit)
        plt.plot(norm_x, poly(norm_x))
        plt.show()
    #return 1/co
    return 1 / (2 * popt[1])


def get_data(kpoints, band, band_edges, steps, band_labels):
    data = []
    for edge in band_edges:
        for x in [-steps, steps]:
            start_tmp = edge['num'] - 1
            if start_tmp + x > 0 and start_tmp + x < len(band):
                end_tmp = start_tmp + x
                start = min([start_tmp, end_tmp])  # start from the lowest
                end = max([start_tmp, end_tmp])  # go to the highest
                direction = kpoints[start] - kpoints[end]
                direction /= norm(direction)

                vals_tmp = band[start:end + 1]
                kpts_tmp = kpoints[start:end + 1]
                vals_rev = vals_tmp[::-1]  # copy the list and reverse it

                if x < 0:
                    vals = vals_rev - vals_tmp[-1]
                    kpts = kpts_tmp[::-1] - kpts_tmp[-1]
                else:
                    vals = vals_tmp - vals_tmp[0]
                    kpts = kpts_tmp - kpts_tmp[0]

                # if going forward then everything is ok, if going backwards
                # then we need to reverse the kpoints

                # now append the lists together in the correct order depending
                # on whether we are moving forwards or backwards. E.g. if we
                # have gone forwards two steps, then we need to put the reversed
                # list behind the current list. We also split the list so that
                # the minima/maxima is not repeated. Axis=0 is used so that the
                # kpoint lists are not merged into one big list when using
                # numpy's append function. E.g. the kpts list stays as a list of
                # lists.
                #if x < 0:
                #    vals = append(vals_tmp, vals_rev[1:])
                #    # This is equal to VBM_kpt - (VBM_kpt-1) + VBM_kpt for each
                #    # step. And then we need to reverse the list
                #    # Note the VBM/CBM appears at the end of the kpts_tmp list
                #    # when we are going backwards
                #    kpts_repeated = 2 * np.array(kpts_tmp[steps]) - kpts_tmp
                #    print kpts_tmp
                #    print kpts_repeated
                #    kpts_rev = kpts_repeated[::-1]
                #    kpts = append(kpts_tmp, kpts_rev[1:], axis=0)
                #else:
                #    vals = append(vals_rev, vals_tmp[1:])
                #    # notice that now the VBM/CBM is at the start of our
                #    # kpts_tmp list.
                #    kpts_repeated = 2 * kpts_tmp[0] - kpts_tmp
                #    kpts_rev = kpts_repeated[::-1]
                #    kpts = append(kpts_rev, kpts_tmp[1:], axis=0)

                # TODO sort label generation
                label = ""
                data.append({'start': start_tmp, 'end': end_tmp, 'vals': vals,
                             'label': label, 'dir': direction,
                             'kpoints': kpts})
    return data


def get_band_edges(kpoints, vb, cb, n_vbm, n_cbm, tol, vkpts, ckpts):
    # check the vbm and cbm given by the Procar
    vbm = {'energy': float("-inf")}
    cbm = {'energy': float("inf")}
    for v_e, c_e in zip(vb, cb):
        if v_e > vbm['energy']:
            vbm['energy'] = v_e
        if c_e < cbm['energy']:
            cbm['energy'] = c_e

    if n_vbm['energy'] != vbm['energy']:
        logging.info("A lower VBM was found outside of the high symmetry band" +
                     " path\ni.e. in the weighted kpoints\n!!!!! CHECK YOUR" +
                     " BAND PATH !!!!!\nenergy: %.5f\tkpoint: %d"
                     % (n_vbm['energy'], n_vbm['kpoint']))
    if n_cbm['energy'] != cbm['energy']:
        logging.info("A lower CBM was found outside of the high symmetry band" +
                     " path\ni.e. in the weighted kpoints\n!!!!! CHECK YOUR" +
                     " BAND PATH !!!!!\nenergy: %.5f\tkpoint: %d"
                     % (n_cbm['energy'], n_cbm['kpoint']))

    # loop through the bands and find other points which are within a tolerance
    # of the vbm and cbm
    vbms = []
    cbms = []
    bg_direct = [{'e': float("inf")}]
    k_num = 1
    for kpoint, v_e, c_e in zip(kpoints, vb, cb):
        if k_num in vkpts or abs(v_e - vbm['energy']) < tol:
            vbms.append({'num': k_num, 'loc': kpoint, 'e': v_e})
        if k_num in ckpts or abs(c_e - cbm['energy']) < tol:
            cbms.append({'num': k_num, 'loc': kpoint, 'e': c_e})

        # calculate the direct bandgaps, will add to the list of direct bandgaps
        # if it is within the tolerance, otherwise if the new band gap is
        # smaller than the current band gap by more than the tolerance it will
        # clear the list and start again.
        bg = c_e - v_e
        if bg_direct[0]['e'] - bg > tol:
            bg_direct = [{'num': k_num, 'loc': kpoint, 'e': bg}]
        elif abs(bg_direct[0]['e'] - bg) < tol:
            bg_direct.append({'num': k_num, 'loc': kpoint, 'e': bg})

        k_num += 1

    logging.info("\nFound the VBM at %d points:" % (len(vbms)))
    for p in vbms:
        logging.info("\t%.3f %.3f %.3f  (kpoint #%d)  with E = %.8f" %
                     (p['loc'][0], p['loc'][1], p['loc'][2], p['num'], p['e']))

    logging.info("Found the CBM at %d points:" % (len(cbms)))
    for p in cbms:
        logging.info("\t%.3f %.3f %.3f  (kpoint #%d)  with E = %.8f" %
                     (p['loc'][0], p['loc'][1], p['loc'][2], p['num'], p['e']))

    bg_indirect = cbms[0]['e'] - vbms[0]['e']
    logging.info("\nIndirect band gap: %.4f" % (bg_indirect))
    logging.info("Direct band gap: %.4f found at %d points:" %
                 (bg_direct[0]['e'], len(bg_direct)))
    for p in bg_direct:
        logging.info("\t%.3f %.3f %.3f  (kpoint #%d)" %
                     (p['loc'][0], p['loc'][1], p['loc'][2], p['num']))
    return vbms, cbms, bg_indirect


def main():
    parser = argparse.ArgumentParser(description="""This script gets information
    such as the nature and size of the band gaps, and the effective masses from
    band structure data.""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('-t', '--tolerance', type=float, default=0.0002,
                        help="""Tolerance for finding the individual VBMs/CBMs
                        in eV. Default is 0.0002""")
    parser.add_argument('-s', '--steps', type=int, default=2,
                        help="""The number of points from the VBM/CBM to sample.
                        Default is 2""")
    parser.add_argument('-v', '--vkpts', type=str,
                        help="""Specific kpoints to include in the calculation
                        of the hole effective masses. List as comma seperated
                        numbers.""")
    parser.add_argument('-c', '--ckpts', type=str,
                        help="""Specific kpoints to include in the calculation
                        of the electron effective masses. List as comma
                        seperated numbers.""")

    args = parser.parse_args()

    logging.basicConfig(filename='bandstats.log', level=logging.DEBUG,
                        filemode='w', format='%(message)s')
    console = logging.StreamHandler()
    logging.getLogger('').addHandler(console)

    poscar = Poscar.load("CONTCAR")
    procar = Procar.find_procars()[0]
    vkpts = map(int, args.vkpts.split(',')) if args.vkpts else []
    ckpts = map(int, args.ckpts.split(',')) if args.ckpts else []
    bandstats(poscar, procar, steps=args.steps, tol=args.tolerance,
              vkpts=vkpts, ckpts=ckpts)

if __name__ == "__main__":
    main()
