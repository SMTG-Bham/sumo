#!/usr/bin/env python

"""
A script to extract the optical absorption data from either Graeme's optics or
standard VASP optics.
"""

__author__ = "Alex Ganose"
__version__ = "0.1"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "Nov 3, 2015"

import sys
import math
import cmath
import argparse

from numpy import array, diag
from scipy.ndimage.filters import gaussian_filter1d
from xml.etree.cElementTree import iterparse


def parse_sarray(data):
    """Helper function for parsing the EPS file.

    Args:
        The lines of the EPS file for a specfic component of the EPS in the form

    Returns:
        The data converted into floats
    """
    return [[float(i) for i in l.split()] for l in data]


def optics(prefix, furthmueller, smear=False, sigma=0.06/2.3548, isymm=6):
    if furthmueller:
        eps = furthmueller_optics(isymm)
    else:
        eps = standard_optics()

    if smear:
        eps = smear_eps(eps)

    abs_data = calc_absorption(eps)

    to_file(prefix + '_optics.dat', abs_data)


def standard_optics():
    """Parse the vasprun.xml file to read the optics data

    Returns:
        A dict containing: {'energies':e, 'real':real, 'imag':imag}
    """
    with open("vasprun.xml", "rt") as f:
        for event, elem in iterparse(f):
            if elem.tag == "dielectricfunction":
                diel = elem
                break

    try:  # TODO: Proper checking for this
        diel
    except:
        print "Couldn't find the absorption data in vasprun.xml"
        print "If you want Forthmueller's optics run with -f"
        sys.exit(0)

    imag = [[float(l) for l in r.text.split()]
            for r in diel.find("imag").find("array").find("set").findall("r")]
    real = [[float(l) for l in r.text.split()]
            for r in diel.find("real").find("array").find("set").findall("r")]

    energies = [e[0] for e in imag]
    real_zip = [e[1:] for e in real]
    imag_zip = [e[1:] for e in imag]

    real_avg = [sum(real_val[:3])/3 for real_val in real_zip]
    imag_avg = [sum(imag_val[:3])/3 for imag_val in imag_zip]
    return {'energies': energies, 'real': real_avg, 'imag': imag_avg}


def furthmueller_optics(isymm):
    """Parse the EPS file produced by the Furthmueller optics binary.

    Args:
        isymm: The value of ISYMM used to produce the EPS output

    Returns:
        A dict containing: {'energies':e, 'real':real, 'imag':imag}
    """
    with open('EPS', 'r') as f:
        eps = f.readlines()
    nedos = len(eps)/6  # this actually equals nedos+1 but it makes things easy

    xx = array(parse_sarray(eps[1:nedos]))
    yy = array(parse_sarray(eps[nedos*1+1:nedos*2]))
    zz = array(parse_sarray(eps[nedos*2+1:nedos*3]))
    xy = array(parse_sarray(eps[nedos*3+1:nedos*4]))
    yz = array(parse_sarray(eps[nedos*4+1:nedos*5]))
    zx = array(parse_sarray(eps[nedos*5+1:nedos*6]))
    energies = array(xx[:, 0])
    real_zip = array(zip(xx[:, 1], yy[:, 1], zz[:, 1],
                         xy[:, 1], yz[:, 1], zx[:, 1]))
    imag_zip = array(zip(xx[:, 2], yy[:, 2], zz[:, 2],
                         xy[:, 2], yz[:, 2], zx[:, 2]))

    # VASP seems to fail pretty badly at diagonalising the EPS matrix using this
    # method, so lets do it ourselves (thanks Alexey Sokol and John Buckeridge)
    real_avg = []
    imag_avg = []
    for energy, real, imag in zip(energies, real_zip, imag_zip):
        real_mat = [[real[0], real[3], real[5]],
                    [real[3], real[1], real[4]],
                    [real[5], real[4], real[2]]]
        imag_mat = [[imag[0], imag[3], imag[5]],
                    [imag[3], imag[1], imag[4]],
                    [imag[5], imag[4], imag[2]]]
        real_avg.append(sum(diag(real_mat))/3)
        imag_avg.append(sum(diag(imag_mat))/3)

    return {'energies': energies, 'real': real_avg, 'imag': imag_avg}


def smear_eps(eps, sigma=0.06/2.3548):
    """Gaussian smear the EPS data by the amount sigma.

    Args:
        eps: EPS data in the form {'energies':e, 'real':real, 'imag':imag}
        sigma: The amount of smearing

    Returns:
        A dict containing: {'energies':e, 'real':real, 'imag':imag}
    """
    e = eps['energies']
    diff = [e[i + 1] - e[i] for i in range(len(e) - 1)]
    diff_avg = sum(diff) / len(diff)
    real = gaussian_filter1d(eps['real'], sigma / diff_avg)
    imag = gaussian_filter1d(eps['imag'], sigma / diff_avg)
    return {'energies': e, 'real': real, 'imag': imag}


def calc_absorption(eps):
    """Calculate the absorption data based on the EPS data. The units of alpha
    is cm^-1.

    Based on http://course.ee.ust.hk/elec342/notes/Lecture%207_absorption%20and%20dispersion.pdf

    Refactive index n has real and imaginary parts:

        n = [(e' + ie"/e_0]^(1/2) = n' + in"

    Relationship between a and imaginary n":

        a = 4*pi*n"/wavelength

    Where wavelength is = hc/E

    Args:
        eps: EPS data in the form {'energies':e, 'real':real, 'imag':imag}

    Returns:
        A list of [energy, alpha, (alpha*hv)^2]
    """

    abs_data = []
    for energy, real, imag in zip(eps['energies'], eps['real'], eps['imag']):
        total = complex(real, imag)
        imag_ref_index = cmath.sqrt(total).imag
        alpha = ((imag_ref_index * energy * 4 * math.pi) / 1.23984212E-6)/100
        ahvsq = (alpha * energy)**2
        abs_data.append([energy, alpha, ahvsq])
    return abs_data


def to_file(filename, data):
    """Write the absorption data to disk.

    Args:
        filename: The name of the file to write the data to
        data: The data to write, in the form [energy, alpha, (alpha*hv)^2]
    """
    with open(filename, 'w') as f:
        f.write("#energy  alpha  (alpha*hv)^2\n")
        for energy, alpha, ahvsq in data:
            f.write("%g %g %g\n" % (energy, alpha, ahvsq))


def main():
    parser = argparse.ArgumentParser(description="""
    optics is a program to parse VASP optics from either the standard method or
    using Furthmueller's method (Graeme's VASP addition). Alpha is given in
    cm^-1.""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('prefix', help='The prefix for the generated files')
    parser.add_argument('-f', '--furthmueller', action='store_true',
                        help="""Enable Furthmueller's optics and look for an
                        EPS file""")
    parser.add_argument('-s', '--smear', action='store_true', help="""Smear the
                        EPS data.""")
    parser.add_argument('--sigma', type=float, default=0.06/2.3548,
                        help="""The smearing value to use, default = 0.025.""")
    parser.add_argument('-i', '--isymm', type=int, default=6,
                        help="""The ISYMM value used in the OPTCTR file (should
                        really be 6). Currently only 6 is supported.""")

    args = parser.parse_args()
    optics(args.prefix, args.furthmueller, args.smear, args.sigma, args.isymm)


if __name__ == "__main__":
    main()
