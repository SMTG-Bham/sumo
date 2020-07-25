"""Read dielectric data from Vasp RPA calculations"""

from xml.etree import ElementTree
import numpy as np

def dielectric_from_vasprun(filename, theory='auto'):
    """Read a vasprun.xml file and return dielectric function from RPA

    This function can interpret calculations performed with ALGO=CHI. For
    regular LOPTICS Vasp runs, use the standard (Pymatgen-based) importer

    Args:
        filename (:obj:`float`):
            Path to ``vasprun.xml`` output of ALGO=CHI calculation.
        theory (:obj:`str`):
            Calculation type to use. Depending on LRPA tag, 'rpa' or 'lfe'
            (local field effects including XC field) is available. If 'auto',
            use whichever of these methods is available. In either case 'ipa'
            (independent particle approximation) is available.

    Returns:
        :obj:`tuple`
            The format imitates the ``dielectric`` attribute of
            :obj:`pymatgen.io.vasp.outputs.Vasprun`: a tuple of the form::

                ([energy1, energy2, ...],
                 [[real_xx_1, real_yy_1, real_zz_1,
                   real_xy_1, real_yz_1, real_xz_1],
                  [real_xx_2, real_yy_2, real_zz_2,
                   real_xy_2, real_yz_2, real_xz_2], ...],
                 [[imag_xx_1, imag_yy_1, imag_zz_1,
                   imag_xy_1, imag_yz_1, imag_xz_1],
                  [imag_xx_2, imag_yy_2, imag_zz_2,
                   imag_xy_2, imag_yz_2, imag_xz_2], ...])
    """

    vrxml = ElementTree.parse(filename)
    root = vrxml.getroot()
    dielectric_functions = root.findall('dielectricfunction')    

    calcs_tags = {'ipa': 'INDEPENDENT PARTICLE',
                  'rpa': 'RPA',
                  'lfe': 'local field effects'}

    data = {}
    for calc, tag in calcs_tags.items():
        try:
            diel_xml = next((diel for diel in dielectric_functions
                            if tag in diel.attrib['comment']))
            data[calc] = {}
            for component in ('real', 'imag'):
                text_rows = ((el.text for el in
                              diel_xml.find(component).find('array').find('set')))
                data[calc][component] = np.genfromtxt(text_rows)
        except StopIteration:
            continue

    theory = theory.lower()
    if theory == 'auto':
        if 'rpa' in data:
            theory = 'rpa'
        elif 'lfe' in data:
            theory = 'lfe'
        else:
            raise Exception("No RPA or LFE data in vasprun")

    return (data[theory]['real'][:, 0],
            data[theory]['real'][:, 1:],
            data[theory]['imag'][:, 1:])
