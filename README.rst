==============
ScriptsForVASP
==============
Making life easier using scripting languages (Bash and Python) to facilitate multiple VASP simulation jobs preparation, submission and analysis.

Introduction
============

There are multiple ways to carry out VASP calculations. Contrary to the Pythonic scheme the famous `ASE <https://wiki.fysik.dtu.dk/ase/>`_ uses
, this project does things in an old-fashioned way, providing a collection of scripts, written in Bash, and accompanying them with
scripts written in Python, in the same directory, with careful use of file names and input arguments, and letting users
add this single directory into their PATH environmental variable to have direct access. The scripts will be used as **executables**. This approach
is less comprehensive than the one ASE uses, but it is more file- and directory-wise convenient thanks to the fact that Bash is a shell
scripting language. Standalone Python scripts (with capitalized first letter) can be imported in an Python interpreter. The plotting Python
scripts (Plot***.py) can be imported and used to adjust the Matplotlib figures.

Scripts structure
=================

Python scripts
--------------

These can also be used alone, just type ``SCRIPT.py -h`` in the terminal you will get help. Additionally, they **can be readily imported and played with in a Python interpreter**.

* ``PlotTDOS.py`` - to plot total density of states, spin states are supported (above or below x axis)
* ``PlotLDOS.py`` - to plot local projected density of states of two ions (above or below x axis), spin states are supported (in different files)
* ``PlotBS.py`` - to plot the band structure along certain chosen points in the 1st BZ
* ``PlotCOHP.py`` - to plot the pCOHP from LOBSTER runs
* ``StrainApplier.py`` - to add strain to the cell lattice vectors in POSCAR for a series of runs, determining the value of a combination of elastic constants
* ``IonsRotator.py`` - to rotate ions in selected lines in POSCAR at a certain point, around a certain axis, to a certain angle

Bash convenient scripts
-----------------------

* ``CellInfo.sh [rwigs]`` - to show the important info from OUTCAR. With the optional rwigs, displays the RWIGS values for different species to 100 % fill the cell
* ``DelFiles.sh [DIRECTORY]`` - to delete the output files, default to current directory
* ``SwapFileNames.sh A B`` - to exchange the names of two files A and B
* ``R [DIRECTORY]`` - to submit the job in a specified directory or default to the current directory
* ``D [N]`` - to cancel the submitted N job, or default to 1 without argument

Bash job submission scripts
---------------------------

* ``Prepare.sh DIRECTORY`` - a fast way to create a directory with necessary input files
* ``SequenceTest.sh TEST_TYPE OPTIONS`` - to create a series of runs, performing a routine test
* ``SequenceDisp.sh TEST_TYPE`` - to do the post-analysis of a series of runs described above and determine physical quantities
* ``Elastic.sh test | disp-solve | solve`` - to do a full set of independent elastic constant determination runs
* ``Electronic.sh prepare | scrun | dosrun | bsrun | lobster kp | lobster test -n NBAND | lobster analysis`` - to carry out electronic-related runs
* ``Bader.sh test | analysis`` - to perform Bader charge transfer analysis
* ``BatchElements.sh`` - a wrapper script to execute the above routines and almost random Bash code for multiple species

Supportive scripts
------------------

**These are used by Bash scripts and better left unvisited.**

* ``_qsub_replacer.sh`` - to do some replacement to the job submission file so that it corresponds to the right working directory. Machine dependent
* ``_sequence_fit_plot.py`` - to perform different fitting techniques, including Birch-Murnaghan equation of state and polynomial
* ``_elastic_solver.py`` - to solve the independent elastic constants from the same number of combination runs with linear algebra
* ``_cell_info_solver.py`` - to obtain the RWIGS values for different species to 100 % fill the cell

General routines
----------------
    HAVEN'T FINALIZED YET BELOW
