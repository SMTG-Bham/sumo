sumo-kgen
==========

``sumo-kgen`` is a program for generating input files containing k-points along high-symmetry
k-point paths, for electronic band structure calculations (e.g. ``KPOINTS`` file for VASP).
The path is unique for each of the 14 Bravais lattice types and, as such, will depend
on the symmetry of the unitcell. As there are several definitions and nomenclatures used in the
literature, ``sumo-kgen`` simplifies this process, allowing for multiple high-symmetry k-point schemes.

.. contents:: Table of Contents
   :local:
   :backlinks: None

Usage
-----

*For simplicity, this tutorial will initially assume you are using VASP. Details for other codes are provided at the end.*

The full range of options supported by ``sumo-kgen`` are detailed in the `Command-Line Interface`_ section,
and be can be accessed using the command::

    sumo-kgen -h

To generate a set of k-points, simply run the following in a folder containing a VASP ``POSCAR`` file::

    sumo-kgen

For example, if we run this command in the ``sumo/tests/data/ZnO`` directory, the k-points will be written to
a file called ``KPOINTS_band``, with the terminal showing the following information::

   Structure information:
      Space group number: 152
      International symbol: P3_121
      Lattice type: hexagonal

   k-point path:
      \Gamma -> A -> L -> M -> \Gamma -> K -> H -> A

   k-points:
      \Gamma: 0.0 0.0 0.0
      A: 0.0 0.0 0.5
      L: 0.0 0.5 0.5
      M: 0.0 0.5 0.0
      K: -0.333 0.667 0.0
      H: -0.333 0.667 0.6

   k-point label indices:
      \Gamma: 1
      A: 21
      L: 74
      M: 94
      \Gamma: 147
      K: 208
      H: 232
      A: 293


Basic Options
~~~~~~~~~~~~~

As the path depends on the symmetry of the lattice, it is important that the symmetry is determined
correctly. The tolerance used for symmetry detection can be controlled using the ``--symprec`` option.
The default is 0.01 Angstrom.

By default, the paths used are those from Bradley and Cracknell [brad]_. To use the k-point paths provided
SeeK-path [seek]_ or pymatgen [curt]_. The options ``--seekpath`` or ``--pymatgen`` can be used.

``sumo-kgen`` automatically looks for a ``POSCAR`` file in the current directory. A different structure
file can be specified using the ``--poscar`` option.

The density of the k-points along the path can be controlled using the ``--density`` option. The default is
60.


Hybrid Band Structures
~~~~~~~~~~~~~~~~~~~~~~

By default, ``sumo-kgen`` generates ``KPOINTS`` files for use in non-self-consistent band structure calculations
(e.g. for use with generalised-gradient approximation functionals). To perform hybrid band structures, the
zero-weighted k-point scheme should be used. To generate ``KPOINTS`` files for use in hybrid band structures, an
``IBZKPT`` file must be located in the current directory (the generated k-points will be appended to those in
this file). Then simply run the following::

    sumo-kgen --hybrid

When generating hybrid k-points, the script will ask if you wish to split the k-points among multiple files,
due to the cost of hybrid band structures which often cannot finish under standard cluster walltimes.
Bear in mind that the total number of k-points per file you choose will not include the addition k-points
included in the ``IBZKPT`` file. To skip this prompt, the number of k-points per file can be specified
via the ``--split`` option.


Folder Generation
~~~~~~~~~~~~~~~~~

Often it is desirable to generate a new folder in which to run the band structure. ``sumo-kgen`` can automate
this procedure using the ``--folders`` (or ``-f``) option. The script will also copy in the required VASP
files from the current directory. Which files are copied depends on the mode of k-point generation.
For example, for non-self-consistent band structures, ``POSCAR``, ``INCAR``, ``POTCAR``, and ``CHGCAR``
will be copied.  For hybrid band structures, only the ``POSCAR``, ``INCAR``, and ``POTCAR`` files will be copied.

If you choose to split the hybrid k-points among a number of ``KPOINTS`` files, a separate folder will be
generated for each file. These will be named ``split-01``, ``split-02``, etc...

The ``sumo-bandplot`` command can automatically detect the presence of these folders and will
reconstruct the full band structure from the individual splits.


Custom k-Point Paths
~~~~~~~~~~~~~~~~~~~~

``sumo-kgen`` also supports generating k-points along custom k-point paths. This is controlled using the
``--kpoints`` option. The custom k-point path is specified as a string with commas separating the k-points.
For example, to generate the k-points along ``0. 0. 0. -> 0.5 0.5 0.5``, the usage is::

    sumo-kgen --kpoints "0 0 0, 0.5 0.5 0.5"

Breaks in the band structure can be indicated using the pipe character.
For example, the path ``0. 0. 0. -> 0.5 0.5 0.5 | 0. 0. 0. -> 0.5 0. 0.``, is specified as::

    sumo-kgen --kpoints "0 0 0, 0.5 0.5 0.5 | 0 0 0, 0.5 0 0"

Custom labels can also be provided using the ``--labels`` option. The syntax is the same as for the
``--kpoints`` option. For example, the labels for the above path are written as::

    sumo-kgen --kpoints "0 0 0, 0.5 0.5 0.5 | 0 0 0, 0.5 0 0" --labels "\Gamma, M | \Gamma, X"

Note: in all cases the arguments are surrounded in parentheses.

Other codes
-----------

CASTEP
~~~~~~

In CASTEP, band structure calculations also include SCF convergence so
two sets of k-points are set: a mesh for the SCF and a path for the
band structure. It is safe to provide *kgen* with a .cell file that
already contains e.g. a *KPOINTS_MP_GRID* tag::

    sumo-kgen --code castep -p seedname.cell

This will write a copy of the cell file to *band.cell*, including a
*BS_KPOINT_LIST* block with the high-symmetry path and with the
special-point labels included as comments. (These comments will help
*bandplot* prettify the x-axis).

Most *kgen* features will work as expected for CASTEP, but the
``--hybrid`` and ``--cartesian`` options are not relevant.
An extra feature is provided to aid phonon calculations with CASTEP: the
``--phonon`` option will write a *PHONON_FINE_KPOINT_LIST* block instead.

Questaal
~~~~~~~~

To perform LMTO band structure calculations the ``lmf`` program can be
given a file defining the band structure path. The crystal structure
is defined with an *init.ext* file (where *ext* is an identifier for
your system) or a *site.ext* file.
Questaal band structures will use the scale factor ALAT set in
*site.ext* which may have been modified from the initial setting,
so it is usually best to read from *site.ext*.
To read the crystal structure and create a band path::

    sumo-kgen --code questaal -p site.ext

will write a file named *syml.ext* ("symmetry lines"); by default this
will use lattice coordinates. To perform the band structure
calculation, specify this file with e.g.::

    lmf -vnit=1 --rs=1,0 --band~mq~fn=syml ext

where the ``~mq`` switch indicates that *syml.ext* is in fractional
coordinates. We recommend avoiding Cartesian coordinates for Questaal
band structures; it is tested and should work but between ALAT scaling
and Bohr units it can get a bit confusing.

Command-Line Interface
----------------------

.. argparse::
   :module: sumo.cli.kgen
   :func: _get_parser
   :prog: sumo-kgen
