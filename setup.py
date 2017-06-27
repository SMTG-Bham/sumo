"""
Vaspy: SMTG utils for working with Vasp
"""

from os.path import abspath, dirname
from setuptools import setup, find_packages

project_dir = abspath(dirname(__file__))

setup(
    name='vaspy',
    version='1.0.0',
    description='Convenience tools for working with VASP',
    long_description="""
Some handy tools developed by and for the Scanlon Materials Theory Group
""",
    url="https://github.com/smtg-ucl/vaspy",
    author="Alex Ganose, Adam J. Jackson",
    author_email="d.scanlon@ucl.ac.uk",
    license='GPL v2',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics'
        ],
    keywords='chemistry ase dft vasp',
    packages=find_packages(),
    install_requires=['ase', 'spglib', 'numpy', 'scipy'],
    entry_points={'console_scripts': [
                      'vaspy-bandgen = vaspy.utils.bandgen:main',
                      'vaspy-bandplot = vaspy.utils.bandplot:main',
                      'vaspy-bandstats = vaspy.utils.bandstats:main',
                      'vaspy-dosplot = vaspy.utils.dosplot:main',
                      'vaspy-kgen = vaspy.utils.kgen:main',
                      'vaspy-optics = vaspy.utils.optics:main']}
    )
