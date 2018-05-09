"""
sumo: Heavy weight plotting tools.
"""


from setuptools import setup, find_packages

import unittest


def load_test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='test*.py')
    return test_suite

with open('README.rst', 'r') as file:
    long_description = file.read()

setup(
    name='sumo',
    version='1.0.8',
    description='Heavy weight plotting tools for ab initio calculations',
    url='https://github.com/smtg-ucl/sumo',
    author='Alex Ganose, Adam J. Jackson',
    author_email='d.scanlon@ucl.ac.uk',
    long_description=long_description,
    license='MIT',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics'
        ],
    keywords='chemistry pymatgen dft vasp dos band',
    test_suite='setup.load_test_suite',
    packages=find_packages(),
    install_requires=['spglib', 'numpy', 'scipy', 'pymatgen>=2017.12.30', 'h5py',
                      'phonopy', 'matplotlib', 'seekpath'],
    extras_require={'docs': ['sphinx', 'sphinx-argparse']},
    package_data={'sumo': ['symmetry/bradcrack.json',
                           'plotting/orbital_colours.conf']},
    data_files=['examples/orbital_colours.conf', 'LICENSE', 'requirements.txt'],
    entry_points={'console_scripts': [
                      'sumo-bandplot = sumo.cli.bandplot:main',
                      'sumo-bandstats = sumo.cli.bandstats:main',
                      'sumo-dosplot = sumo.cli.dosplot:main',
                      'sumo-kgen = sumo.cli.kgen:main',
                      'sumo-phonon-bandplot = sumo.cli.phonon_bandplot:main',
                      'sumo-optplot = sumo.cli.optplot:main']}
    )
