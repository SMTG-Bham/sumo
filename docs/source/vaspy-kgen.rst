vaspy-kgen
==========


Background
----------

``vaspy-kgen`` is a program for generating VASP ``KPOINTS`` files containing k-points along high-symmetry
k-point paths, for electronic band structure calculations.
The path is unique for each of the 14 Bravais lattice types and, as such, will depend
on the symmetry of the unitcell. As there are several definitions and nomenclatures used in the
literature, ``vaspy-kgen`` simplifies this process, allowing for multiple high-symmetry k-point schemes.


Usage
-----

The full range of options supported by ``vaspy-kgen`` can be accessed using the command::

    vaspy-kgen -h

To generate a ``KPOINTS`` file, simply run the following in a folder containing a VASP ``POSCAR`` file::

    vaspy-kgen

For example, if we run this command in the ``vaspy/data/ZnO`` directory, the k-points will be written to
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

   k-point label indicies:
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
The default is 0.001 Angstrom.

By default, the paths used are those from Bradley and Cracknell [brad]_. To use the k-point paths provided
SeeK-path [seek]_ or pymatgen [curt]_. The options ``--seekpath`` or ``--pymatgen`` can be used.

``vaspy-kgen`` automatically looks for a ``POSCAR`` file in the current directory. A different structure
file can be specified using the ``--poscar`` option.

The denisty of the k-points along the path can be controlled using the ``--density`` option. The default is
60.


Hybrid Band Structures
~~~~~~~~~~~~~~~~~~~~~~

By default, ``vaspy-kgen`` generates ``KPOINTS`` files for use in non-selfconsistent band structure calculations (e.g. for use with generalised-gradient approximation functionals). To perform hybrid band structures, the
zero-weighted k-point scheme should be used. To generate ``KPOINTS`` files for use in hybrid band structures, an
``IBZKPT`` file must be located in the current directory (the generated k-points will be appended to those in
this file). Then simply run the following::

    vaspy-kgen --hybrid

When generating hybrid k-points, the script will ask if you wish to split the k-points among multiple files,
due to the cost of hybrid band structures which often cannot finish under standard cluster walltimes.
Bear in mind that the total number of k-points per file you choose will not include the addition k-points
included in the ``IBZKPT`` file. To skip this prompt, the number of k-points per file can be specified
via the ``--split`` option.


Folder Generation
~~~~~~~~~~~~~~~~~

Often it is desirable to generate a new folder in which to run the band structure. ``vaspy-kgen`` can automate
this procedure using the ``--folders`` (or ``-f``) option. The script will also copy in the required VASP
files from the current directory. Which files are copied depends on the mode of k-point generation.
For example, for non-selfconsistent band structures, ``POSCAR``, ``INCAR``, ``POTCAR``, and ``CHGCAR``
will be copied.  For hybrid band structures, only the ``POSCAR``, ``INCAR``, and ``POTCAR`` files will be copied.

If you choose to split the hybrid k-points among a number of ``KPOINTS`` files, a seperate folder will be
generated for each file. These will be named ``split-01``, ``split-02``, etc...

The ``vaspy-bandplot`` command can automatically detect the presence of these folders and will
reconstruct the full band structure from the individual splits.


Custom k-Point Paths
~~~~~~~~~~~~~~~~~~~~

``vaspy-kgen`` also supports generating k-points along custom k-point paths. This is controlled using the
``--kpoints`` option. The custom k-point path is specified as a string with commas separating the k-points.
For example, to generate the k-points along ``0. 0. 0. -> 0.5 0.5 0.5``, the usage is::

    vaspy-kgen --kpoints "0 0 0, 0.5 0.5 0.5"

Breaks in the band structure can be indicated using the pipe character.
For example, the path ``0. 0. 0. -> 0.5 0.5 0.5 | 0. 0. 0. -> 0.5 0. 0.``, is specified as::

    vaspy-kgen --kpoints "0 0 0, 0.5 0.5 0.5 | 0 0 0, 0.5 0 0"

Custom labels can also be provided using the ``--labels`` option. The syntax is the same as for the
``--kpoints`` option. For example, the labels for the above path are written as::

    vaspy-kgen --kpoints "0 0 0, 0.5 0.5 0.5 | 0 0 0, 0.5 0 0" --labels "\Gamma, M | \Gamma, X"

Note: in all cases the arguments are surrounded in parentheses.
