import os
import re
import sys
import logging
import itertools

from vaspy.misc.math_func import multl_const, mag, angle


class Incar(dict):
    """
    This is more or less copied from pymatgen because it does the job
    """

    def __init__(self, params=None):
        super(Incar, self).__init__()
        if params:
            self.update(params)

    def __setitem__(self, key, val):
        # strip and process values before adding them
        super(Incar, self).__setitem__(
            key.strip(), Incar.proc_val(key.strip(), val.strip())
            if isinstance(val, str) else val)

    @staticmethod
    def from_file(filename):
        with open(filename, "rt") as f:
            return Incar.from_string(f.read())

    @staticmethod
    def from_string(string):
        lines_a = list(strip_comments(re.split(r'[\r|\n|\n\r]+', string)))
        lines = list(strip_comments(re.split(r'[\r|\n|\n\r|;]+', "\n".join(lines_a))))
        params = {}
        for line in lines:
            m = re.match("(\w+)\s*=\s*(.*)", line)
            if m:
                key = m.group(1).strip()
                val = m.group(2).strip()
                val = Incar.proc_val(key, val)
                params[key] = val
        return Incar(params)

    @staticmethod
    def proc_val(key, val):
        """
        Static helper method to convert INCAR parameters to proper types, e.g.,
        integers, floats, lists, etc.
        Args:
            key: INCAR parameter key
            val: Actual value of INCAR parameter.
        """
        list_keys = ("LDAUU", "LDAUL", "LDAUJ", "MAGMOM")
        bool_keys = ("LDAU", "LWAVE", "LSCALU", "LCHARG", "LPLANE",
                     "LHFCALC", "ADDGRID")
        float_keys = ("EDIFF", "SIGMA", "TIME", "ENCUTFOCK", "HFSCREEN",
                      "POTIM", "EDIFFG")
        int_keys = ("NSW", "NBANDS", "NELMIN", "ISIF", "IBRION", "ISPIN",
                    "ICHARG", "NELM", "ISMEAR", "NPAR", "LDAUPRINT", "LMAXMIX",
                    "ENCUT", "NSIM", "NKRED", "NUPDOWN", "ISPIND", "LDAUTYPE")
        caps_keys = ("GGA")

        def smart_int_or_float(numstr):
            if numstr.find(".") != -1 or numstr.lower().find("e") != -1:
                return float(numstr)
            else:
                return int(numstr)

        try:
            if key in list_keys:
                output = []
                toks = re.findall(r"(-?\d+\.?\d*)\*?(-?\d+\.?\d*)?", val)
                for tok in toks:
                    if tok[1]:
                        output.extend([smart_int_or_float(tok[1])]
                                      * int(tok[0]))
                    else:
                        output.append(smart_int_or_float(tok[0]))
                return output

            if key in "ROPT":
                return([float(i) for i in val.split()])

            if key in bool_keys:
                m = re.match(r"^\.?([T|F|t|f])[A-Za-z]*\.?", val)
                if m:
                    if m.group(1) == "T" or m.group(1) == "t":
                        return True
                    else:
                        return False
                raise ValueError(key + " should be a boolean type!")

            if key in float_keys:
                return float(re.search(r"^-?\d*\.?\d*[e|E]?-?\d*", val)
                             .group(0))

            if key in int_keys:
                return int(re.match(r"^-?[0-9]+", val).group(0))

            if key in caps_keys:
                return val.upper()

        except ValueError:
            pass

        # Not in standard keys. We will try a hirerachy of conversions.
        try:
            val = int(val)
            return val
        except ValueError:
            pass
        try:
            val = float(val)
            return val
        except ValueError:
            pass

        if "true" in val.lower():
            return True

        if "false" in val.lower():
            return False
        try:
            if key not in ("TITEL", "SYSTEM"):
                return re.search(r"^-?[0-9]+", val.capitalize()).group(0)
            else:
                return val.capitalize()
        except:
            return val.capitalize()

    def to_file(self, filename):
        """
        Write Incar to a file.
        Args:
            filename (str): filename to write to.
        """
        with open(filename, "wt") as f:
            f.write(self.get_string())


    def write_file(self, filename):
        """
        Write Incar to a file. This is for compatability reasons with PMG
        Args:
            filename (str): filename to write to.
        """
        self.to_file(filename)


    def get_string(self, sort_keys=False, formatted=True):
        """
        Returns a string representation of the INCAR

        Args:
            sort_keys (bool): Set to True to sort the INCAR parameters
                alphabetically. Defaults to False.
            formatted (bool): Set to True for pretty aligned output. Defaults
                to False.
        """
        keys = self.keys()
        if sort_keys:
            keys = sorted(keys)

        lines = []
        lines_dict = {}
        for k in keys:
            if k == "MAGMOM" and isinstance(self[k], list):
                value = []
                for m, g in itertools.groupby(self[k]):
                    n = len(tuple(g))
                    if n is 1:
                        value.append("{}".format(m))
                    else:
                        value.append("{}*{}".format(n, m))
                lines.append([k, " ".join(value)])
                lines_dict[k] = " ".join(value)
            elif isinstance(self[k], list):
                lines.append([k, " ".join([str(i) for i in self[k]])])
                lines_dict[k] = " ".join([str(i) for i in self[k]])
            else:
                lines.append([k, self[k]])
                lines_dict[k] = self[k]

        if formatted:
            return format_incar(lines_dict) + "\n"
        else:
            return "\n".join([" = ".join([str(m) for m in line])
                              for line in lines])

    @property
    def hse06(self):
        if 'GGA' in self and self['GGA'] == 'PS':
            return False
        if 'LHFCALC' in self and self['LHFCALC'] is True:
            return True
        else:
            return False

    @property
    def pbesol(self):
        if 'GGA' in self and self['GGA'] == 'PS' and not self.hse06:
            return True
        else:
            return False

    def diff(self, other):
        """
        Diff function for Incar.  Compares two Incars and indicates which
        parameters are the same and which are not. Useful for checking whether
        two runs were done using the same parameters.
        Args:
            other (Incar): The other Incar object to compare to.
        Returns:
            Dict of the following format:
            {"Same" : parameters_that_are_the_same,
            "Different": parameters_that_are_different}
            Note that the parameters are return as full dictionaries of values.
            E.g. {"ISIF":3}
        """
        similar_param = {}
        different_param = {}
        for k1, v1 in self.items():
            if k1 not in other:
                different_param[k1] = {"INCAR1": v1, "INCAR2": None}
            elif v1 != other[k1]:
                different_param[k1] = {"INCAR1": v1, "INCAR2": other[k1]}
            else:
                similar_param[k1] = v1
        for k2, v2 in other.items():
            if k2 not in similar_param and k2 not in different_param:
                if k2 not in self:
                    different_param[k2] = {"INCAR1": None, "INCAR2": v2}
        return {"Same": similar_param, "Different": different_param}


def format_incar(lines_dict, show_hidden=True, comments=True):
    """
    Returns the formatted INCAR exactly how I like it to look.

    This is a bit of a mess but it works pretty nicely
    Args:
        lines_dict (dict): All the tags with pre formatted values.
        show_hidden (bool): Show the tags that are missing from my INCAR
        show_comments (bool): Show relevant comments for the tags
    """

    titles = ("start parameters", "parallelisation", "electronic relaxation",
              "ionic relaxation", "misc", "hybrid-dft", "magnetic", "dft+u",
              "decomposed charge density")

    start = (('NWRITE', '2', '! Medium-level information output'),
             ('ISTART', '1', '! read existing wavefunction; if there'),
             ('INIWAV', '1', '! Random initial wavefunction; otherwise'),
             ('ICORELEVEL', '1', '! Print core levels'),
             ('ICHARG', '11', '! Non-selfconsistent: GGA/LDA band structures'),
             ('NBANDS', '35', '! No. bands'),
             ('NELECT', '35', '! No. bands'))

    para = (('NCORE', '12', '! No. cores per orbital'),
            ('LPLANE', '', '! Real space distribution; supercells'),
            ('KPAR', '4', "! k-point parallelisation"))

    elec = (('PREC', 'Accurate', '! Precision level'),
            ('ALGO', 'Fast', '! SCF minimisation algorithm; 38/48 combo'),
            ('ENMAX', '500', '! Plane-wave cutoff'),
            ('NELM', '200', '! Max SCF steps'),
            ('NELMIN', '2', '! Min SCF steps'),
            ('EDIFF', '1E-05', '! SCF energy convergence'),
            ('GGA', 'PS', '! PBEsol exchange-correlation'),
            ('LASPH', 'True', '! Non-spherical elements; d/f convergence'),
            ('LREAL', 'Auto', '! Projection operators: automatic'),
            ('ADDGRID', 'True', '! Increase grid; helps GGA convergence'))

    ioni = (('EDIFFG', '-0.01', '! Ionic convergence; eV/AA^3'),
            ('NSW', '200', '! Max ionic steps'),
            ('IBRION', '1', '! Algorithm: 0-MD; 1-Quasi-New; 2-CG'),
            ('ISIF', '3', '! Stress/relaxation: 2-Ions, 3-Shape/Ions/V, 7-Vol'),
            ('ISYM', '2', '! Symmetry: 0-none; 2=GGA; 3=hybrids'),
            ('NBLOCK', '1', '! Update XDATCAR every X steps'),
            ('KBLOCK', '40', '! Update PCDAT and DOSCAR every X*NBLOCK steps'),
            ('ISMEAR', '0', '! Gaussian smearing; metals:1'),
            ('SIGMA', '0.05', '! Smearing value in eV; metals:0.2'),
            ('IWAVPR', '1', '! charge density extrapolation: '
             '0-non 1-charg 2-wave 3-comb'),
            ('POTIM', '0.04', '! Timestep in fs'))

    misc = (('LORBIT', '11', '! PAW radii for projected DOS'),
            ('NEDOS', '2000', '! DOSCAR points'),
            ('LVHAR', 'True', '! Ionic and Hartree potential'),
            ('LORBIT', '1', '! Supply radii for projected DOS'),
            ('RWIGS', '1.5 1.5', '! Radii for each atomic species'),
            ('LOPTICS', 'True', '! Output OPTIC file'),
            ('LVTOT', 'True', '! Electrostatic potential'),
            ('LELF', 'True', '! Localisation function'))

    hybr = (('LHFCALC', 'True', '! Activate HF'),
            ('PRECFOCK', 'Fast', '! HF FFT grid'),
            ('NKRED', '2', '! Reduce k-grid-even only'),
            ('ALGO', 'All', '! SCF Combo; ALGO=58; Damped for MD'),
            ('TIME', '0.30', '! Timestep for IALGO5X'),
            ('HFLMAX', '4', '! HF cut-off: 4d,6f'),
            ('HFSCREEN', '0.207', '! Switch to screened exchange; e.g. HSE06'),
            ('AEXX', '0.25', '! 25% HF exchange; e.g. PBE0'))

    magn = (('ISPIN', '2', '! Enable spin polarisation'),
            ('MAGMOM', '5 0', '! Initial magnetic momoment on each ion'),
            ('NUPDOWN', '-1', '! Enforce spin multiplet'),
            ('LSORBIT', 'True', '! Spin-orbit coupling'))

    dftu = (('LDAU', 'True', '! Activate DFT+U'),
            ('LDATYPE', '2', '! Dudarev; only U-J matters'),
            ('LDAUL', '2 -1', '! Orbitals for each species'),
            ('LDAUU', '2  0', '! U for each species'),
            ('LDAUJ', '0  0', '! J for each species'),
            ('LMAXMIX', '4', '! Mixing cut-off; 4-d, 6-f'))

    part = (('LPARD', 'True', '! Generate PARCHG'),
            ('EINT', '-10 0', '! Energy range'),
            ('NBMOD', '-3', '! With reference to Ef'),
            ('KPUSE', '1', '! Over k-points'),
            ('IBAND', '20', '! Over bands'))

    def add_section(params, length=23):
        sec_string = []
        hidden_string = []
        max_len = max(map(len, [k[0] for k in params]))
        for key, val, comment in params:
            # ALGO is in two sections, want to show comments if it is in the
            # other section
            print_algo = False
            if (key == 'ALGO' and val == 'All' and
                    key in lines_dict and lines_dict[key] in 'Damped All'):
                print_algo = True
            elif (key == 'ALGO' and val != 'All' and
                    key in lines_dict and lines_dict[key] not in 'Damped All'):
                print_algo = True

            if key in lines_dict and ((key == 'ALGO' and print_algo) or
                                      key != 'ALGO'):
                spaces = " " * (max_len - len(key))
                s = "  %s%s = %s" % (key, spaces, str(lines_dict[key]))
                # these lines are usually too long for comments
                if key + ' ' not in 'MAGMOM LDAUL LDAUU LDAUJ ' and comments:
                    s += ' ' * (length - len(s))
                    s += comment

                del lines_dict[key]
                sec_string.append(s)
            elif show_hidden:
                spaces = " " * (max_len - len(key))
                s = "  !%s%s = %s" % (key, spaces, val)
                # we want to make sure you can see when lines are commented out
                s += " " * (1 + length - len(s))
                if comments:
                    s += comment
                hidden_string.append(s)

        return sec_string + hidden_string

    if 'SYSTEM' in lines_dict:
        s = ["SYSTEM = " + lines_dict["SYSTEM"]]
        del lines_dict["SYSTEM"]

    for title, params in zip(titles, [start, para, elec, ioni, misc,
                                      hybr, magn, dftu, part]):
        sec = add_section(params)
        if sec:
            s.append("\n" + title)
            s += sec

    if lines_dict:  # we missed some tags
        s.append("\nothers")
        max_len = max(itertools.imap(len, lines_dict))
        for key in lines_dict:
            spaces = " " * (max_len - len(key))
            string = "  %s%s = %s" % (key, spaces, str(lines_dict[key]))
            s.append(string)

    return "\n".join(s)
