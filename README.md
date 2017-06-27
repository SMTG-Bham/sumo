y# vaspy
A collection of python scripts to make working with VASP a little easier.

There is also a small python library available which handles reading many
VASP input and OUTPUT scripts;
it is intended that this functionality will gradually by offloaded to the ASE project.

There are probably still some bugs. If you think you've find any,
please report them on the Issue Tracker. The worst case is that you've
found a need to improve the documentation!

Installation
===========

Check the SMTG wiki page for a more detailed account of setting up your SSH keys.

To install vaspy, first clone the repository to somewhere on the server.
This is source code, so a folder called `src` would be traditional.

```bash
cd ~/src
git clone git@github.com:smtg-ucl/vaspy.git
```

Make sure you have a proper Python setup including packaging tools. For systems used in SMTG, see the `module` commands in **Dependencies**, below.

Now we use `pip` to install to your `.local` folder. 
This is safer than a system-wide installation; in any case, on a remote cluster you only have permission to do a user install.


```bash
cd vaspy
pip install --user .
```

Your local `bin` folder and Python `lib` folder *should* already be on your path; check with `echo $PATH` and `echo $PYTHONPATH`. If not, then you should add them to your path with 

```bash
export PATH="$PATH:~/.local/bin"
export PYTHONPATH="$PYTHONPATH:~/.local/lib/Python2/site-packages"
```

(the exact location may vary depending on your Python setup.) Add these lines to your `.bashrc` so you don't have to type them every time you log in.

### Developer installation

Developers may prefer to install using `pip install --user -e .` which
creates an "editable" local installation. Instead of copying files,
this creates links to the source folder so that that tweaks to the
code in your source folder will be immediately reflected on the PATH.


Updating
--------

To update, we pull the main branch code from Git and re-copy the files
to the folders on your PATHs.

```bash
cd ~/src/vaspy
git pull
pip install --user .
```



Dependencies
============

This package requires:
 - Numpy
 - Scipy
 - Matplotlib
 - spglib
 - Atomic simulation environment
 - Pyhull (for calculating stability above hull)
 - Xmgrace (for visualising generated bandstructures)

On Archer these dependencies can be loaded by adding the following to your `.bashrc`:

```bash
module load xmgrace
module load anaconda
```

On Legion/Grace add:

```bash
module unload python
module load python2
module load grace
```

The non-standard Python libraries spglib and ASE will be installed automatically by `pip`.
