#!/usr/bin/env python
import logging
import argparse
import sys
import os
import re

from vaspy.vasp_output import Procar


class BandgenOptions(object):

    def __init__(self, options={}):
        if options:
            self.options = options
        else:
            parser = argparse.ArgumentParser(
                description="generate band structure data")
            parser.add_argument("-prefix", help="prefix for generated files")
            parser.add_argument("-file_type", type=int,
                                choices=[0, 1, 2], default=0,
                                help="possible type for PROCAR either both "
                                "(0), pbe (1), hybrid (2) - default 0")

            comb = parser.add_mutually_exclusive_group()
            comb.add_argument("-combine", action="store_true",
                              help="look for PROCARS in the current directory "
                              "as PROCARXX and combine")
            comb.add_argument("-search", action="store_true",
                              help="look for band_X(_SPLIT_Y)_PBE/HSE "
                              "folders and combine")

            args = parser.parse_args()

            prefix = args.prefix+"_" if args.prefix else ""

            if args.file_type == 0:
                suffix = "PBE|HSE"
            elif args.file_type == 1:
                suffix = "PBE"
            else:
                suffix = "HSE"

            self.options = {'prefix': prefix, 'suffix': suffix,
                            'combine': args.combine, 'search': args.search}


class Bandgen(object):

    def __init__(self, bandgen_options):
        self.prefix = bandgen_options.options['prefix']
        self.suffix = bandgen_options.options['suffix']
        self.combine = bandgen_options.options['combine']
        self.search = bandgen_options.options['search']

    def _find_procars(self):
        if self.search:
            logging.info("searching for matching directories")

            dirs = next(os.walk('.'))[1]

            r = re.compile('band_(\d*)_(%s)' % self.suffix, re.I)
            r_split = re.compile('band_(\d*)_split_(\d*)_(%s)' %
                                 self.suffix, re.I)
            nosplit_band_dirs = filter(r.match, dirs)
            split_band_dirs = filter(r_split.match, dirs)

            if nosplit_band_dirs and split_band_dirs:
                sys.exit("both split and non split band paths found, "
                         "too confusing")
            elif not nosplit_band_dirs and not split_band_dirs:
                sys.exit("no matching folders found")

            band_dirs = split_band_dirs + nosplit_band_dirs
            band_dirs.sort()

            band_dirs = sorted(band_dirs, key=lambda x:
                    int(x.split('_')[3]) if len(x.split('_')) == 5 else
                    int(x.split('_')[1]))
            print band_dirs

            filenames = []
            for folder in band_dirs:
                filename = folder + "/PROCAR"
                if os.path.isfile(filename):
                    logging.info("found PROCAR in %s", folder)
                else:
                    sys.exit("no PROCAR in %s" % folder)

                # Organise the filenames so we can use them more easily later
                num = int(folder.split("_")[1])
                if len(filenames) >= num:
                    filenames[num-1].append(filename)
                else:
                    filenames.append([filename])

            logging.info("\n%d band structure paths found from a total of %d "
                         "PROCARs", len(filenames), sum(map(len, filenames)))

        elif self.combine:
            logging.info("looking for PROCARS to combine")
            files = [f for f in os.listdir('.') if os.path.isfile(f)]
            r = re.compile('PROCAR\d+')
            filenames = [filter(r.match, files)]
            filenames[0].sort()

            if filenames[0]:
                logging.info("found %d PROCARs", len(filenames[0]))
            else:
                sys.exit("no matching folders found")

        else:
            if os.path.isfile("PROCAR"):
                filenames = [["PROCAR"]]
                logging.info("loaded PROCAR")
            else:
                sys.exit("no PROCAR in current directory")

        return filenames

    def generate_bandstructure(self):
        filenames = self._find_procars()

        self.vbm = {'band': 0, 'energy': float('-inf'), 'string': ""}
        self.cbm = {'band': 0, 'energy': float('inf'), 'string': ""}

        generated_files = []
        for path in filenames:
            index = filenames.index(path) + 1
            logging.info("\nfor path %d", index)

            if len(path) > 1:
                procar = Procar.combine(path)
            else:
                procar = Procar.load(filename=path[0])

            logging.info("   num k-points: %d\tnum bands: %d\tnum ions: %d",
                         procar.num_kpoints, procar.num_bands, procar.num_ions)

            if procar.vbm['energy'] > self.vbm['energy']:
                self.vbm = {'band': procar.vbm['band'],
                            'energy': procar.vbm['energy'],
                            'string': procar.str_VBM()}
            if procar.cbm['energy'] < self.cbm['energy']:
                self.cbm = {'band': procar.cbm['band'],
                            'energy': procar.cbm['energy'],
                            'string': procar.str_CBM()}

            filename = self.prefix + "band_" + str(index) + ".dat"
            procar.to_file(filename)
            logging.info("   saved data as %s",filename)
            generated_files.append(filename)

        logging.info("\nhighest occupied band is: %d", self.vbm['band'])
        logging.info("lowest unoccupied band is: %d", self.cbm['band'])
        logging.info("\n" + self.vbm['string'])
        logging.info(self.cbm['string'])
        logging.info("band gap is: %.3f eV",
                     (self.cbm['energy'] - self.vbm['energy']))

        return generated_files

def main():
    logging.basicConfig(filename='bandgen.log', level=logging.DEBUG,
                        format='%(message)s', filemode='w')
    console = logging.StreamHandler()
    logging.getLogger('').addHandler(console)
    options = BandgenOptions()
    bandgen = Bandgen(options)
    bandgen.generate_bandstructure()
    sys.exit(0)


if __name__ == '__main__':
    main()
