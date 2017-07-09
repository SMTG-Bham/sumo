#!/usr/bin/env python
import sys
import os
import shutil
import argparse
import logging
import math

from vaspy.io.vasp_input import Poscar, Kpoints
from vaspy.electronic_structure.band import Symmetry, PathGen


def main():

    parser = argparse.ArgumentParser(description=(
        "Generate high symmetry k-point paths for a crystal "
        "structure. Non-self-consistent (preferred for GGA) and"
        " self-consistent (preferred for hybrid DFT) modes are"
        " supported."))

    sym_group = parser.add_mutually_exclusive_group()
    sym_group.add_argument(
        "-tolerance",
        type=int,
        choices=[1, 2, 3, 4, 5],
        default=2,
        help=("Number of decimal places for symmetry "
              "analysis tolerance - default: 2"))
    sym_group.add_argument(
        "-symmetry",
        help="Override kgen symmetry determination."
        "This should be in the form of \"centering,type\", "
        "e.g. \"P,orthorhombic\" ")
    parser.add_argument(
        "-density",
        type=int,
        default=400,
        help="k-point line density - default: 400")
    parser.add_argument(
        "-breakup",
        type=int,
        default=0,
        help="Number of kpoints per file for hybrid calculations")
    parser.add_argument(
        "-ybrid",
        action="store_true",
        help="Generate k-points file for self-consistent band structure (zero-"
        "weighting along paths); this is the preferred approach for DFT "
        "with hybrid exchange/correlation functionals.")
    parser.add_argument(
        "-all",
        action="store_true",
        help="Make KPOINTS files for non-self-consistent (typical GGA) and "
        "self-consistent (needed for hybrid) calculations.")
    parser.add_argument(
        "-folders",
        action="store_true",
        help=("Generate folders and copy files into them "
              "(POSCAR, INCAR, POTCAR, job)."))
    parser.add_argument(
        "-old",
        action="store_true",
        help="Don't use the new path for Orthorhombic p. "
        "(This exists for compatibility with old projects).")

    args = parser.parse_args()

    poscar = Poscar.load("CONTCAR")
    logging.basicConfig(
        filename='kgen.log',
        level=logging.DEBUG,
        format='%(message)s',
        filemode='w')
    console = logging.StreamHandler()
    logging.getLogger('').addHandler(console)

    if args.symmetry is None:
        sym = Symmetry.from_poscar(poscar, tolerance=args.tolerance)
        tol = 1.0 / (10 ** args.tolerance)
        logging.info("Tolerance = %g, %s" % (tol, sym))
        logging.info("Conventional cell %s" % sym.str_conv_cell())
    else:
        split_sym = args.symmetry.split(",")
        if len(split_sym) != 2:
            sys.exit("symmetry format invalid. "
                     "Should be e.g. \"P,orthorhombic\"")
        centre, system = split_sym
        logging.info("specified symmetry is {0}-centred"
                     " {0} lattice".format(centre, system))
        logging.info("this method probably isn't very good")

        # Here we trust people to know what they are doing and
        # use the poscar as a conventional cell.
        # What could possibly go wrong?
        conv_cell = poscar.lattice_abc + poscar.lattice_angles
        sym = Symmetry({'spacegroup_symbol': centre,
                        'crystal_system': system}, conv_cell)

    hybrid = 0
    if args.ybrid:
        hybrid = 1
        try:
            ibzkpt = Kpoints.load_ibzkpt()
            logging.info("Hybrid calculation specified and found IBZKPT, "
                         "will add to KPOINTS file")
        except:
            logging.info("No IBZKPT file was found; please remember to "
                         "append the k-points manually.")
            ibzkpt = None

    else:
        ibzkpt = None

    if args.all:
        hybrid = 2

    try:
        path = PathGen(poscar, sym, density=args.density, old=args.old)
    except Exception as e:
        print e
        sys.exit("Exiting")

    logging.info("%s\n" % path.message)
    logging.info("Determined Brillouin zone path is: " + path.str_labels())
    logging.info("There are {0} k-points, the coordinates of "
                 "which are: ".format(path.str_nspecial_kpoints()))
    logging.info(path.str_special_kpoints())
    logging.info("\n" + path.str_nkpoints_per_step())

    breakup = 0
    if not args.breakup and args.ybrid:
        logging.info("\nThere are a total of {0} kpoints, do you want to "
                     "break them up? (y/n)".format(path.str_nkpoints()))
        if raw_input()[0].lower() == 'y':
            logging.info("How many kpoints per file?")
            breakup = input()
            try:
                int(breakup)
            except:
                print "not a valid number"
                sys.exit()
    elif args.breakup:
        breakup = args.breakup

    filenames = path.to_files(hybrid=hybrid, split=breakup, ibzkpt=ibzkpt)

    if args.folders:
        vasp_files = []
        for vasp_file in [poscar.filename, "INCAR", "POTCAR", "job", "CHGCAR"]:
            if os.path.isfile(vasp_file):
                vasp_files.append(vasp_file)

        logging.info("Checking for VASP files in current directory to copy...")
        logging.info("Found: " + ", ".join(vasp_files))

        for kpoint_file in filenames:
            new_folder = kpoint_file.replace("KPOINTS", "band").lower()
            os.makedirs(new_folder)
            shutil.move(kpoint_file, new_folder + '/KPOINTS')

            for vasp_file in vasp_files:
                if "CHGCAR" in vasp_file and "PBE" not in kpoint_file:
                    continue
                elif vasp_file == poscar.filename:
                    filename = "POSCAR"
                else:
                    filename = vasp_file
                shutil.copy(vasp_file, new_folder + '/' + filename)


if __name__ == '__main__':
    main()
