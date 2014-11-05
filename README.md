ScriptsForVASP
==============

Making life easier using scripting languages (Bash and Python) to
facilitate multiple VASP simulation jobs preparation, submission and
analysis.

### NOTE

Due to the specific workflow related nature of this project, it did not 
and perhaps will not be of too much use to most VASP users. Since it has come 
to a level handy enough to suit my own research need, I will not actively 
maintain this project anymore, and the documentation will not be added either. 
However, I selected some of the generally useful features, like plotting, 
fitting, POSCAR manipulation and reorganized them, made a Python package out of it, 
named [pydass_vasp](http://terencezl.github.io/pydass_vasp/). It should be 
interesting to you. Just take a look.

Introduction
------------

There are multiple ways to carry out VASP calculations. Contrary to the
Pythonic scheme the famous [ASE](https://wiki.fysik.dtu.dk/ase/) uses,
this project does things in an old-fashioned way, providing a collection
of scripts, written in Bash, and accompanying them with scripts written
in Python, in the same directory, with careful use of file names and
input arguments, and letting users add this single directory into their
PATH environmental variable to have direct access. The scripts will be
used as **executables**. This approach is less comprehensive than the
one ASE uses, but it is more file- and directory-wise convenient thanks
to the fact that Bash is a shell scripting language. Standalone Python
scripts (with capitalized first letter) can be imported in an Python
interpreter.

Scripts structure
-----------------

### Python scripts

-   `StrainApplier.py` - to add strain to the cell lattice vectors in
    POSCAR for a series of runs, determining the value of a combination
    of elastic constants
-   `IonsRotator.py` - to rotate ions in selected lines in POSCAR at a
    certain point, around a certain axis, to a certain angle

### Bash convenient scripts

-   `CellInfo.sh [rwigs]` - to show the important info from OUTCAR. With
    the optional rwigs, displays the RWIGS values for different species
    to 100 % fill the cell
-   `DelFiles.sh [DIRECTORY]` - to delete the output files, default to
    current directory
-   `SwapFileNames.sh A B` - to exchange the names of two files A and B
-   `M [DIRECTORY]` - to submit the job in a specified directory or
    default to the current directory

### Bash job submission scripts

-   `Prepare.sh DIRECTORY` - a fast way to create a directory with
    necessary input files
-   `SequenceTest.sh TEST_TYPE [OPTIONS]` - to create a series of runs,
    performing a routine test
-   `SequenceDisp.sh TEST_TYPE` - to do the post-analysis of a series of
    runs described above and determine physical quantities
-   `Elastic.sh TEST_TYPE CRYST_SYS [OPTIONS]` - to do a full set of
    independent elastic constant determination runs
-   `Electronic.sh TEST_TYPE [OPTIONS]` - to carry out
    electronic-related runs
-   `BatchElements.sh [OPTIONS] 'COMMANDS'` - a wrapper script to
    execute the above routines and almost random Bash code for multiple
    species

### Supportive scripts

**These are used by Bash scripts and better left unvisited.**

-   `_qsub_replacer.sh` - to do some replacement to the job submission
    file so that it corresponds to the right working directory. Machine
    dependent
-   `_sequence_fit_plot.py` - to perform different fitting techniques,
    including Birch-Murnaghan equation of state and polynomial
-   `_elastic_solver.py` - to solve the independent elastic constants
    from the same number of combination runs with linear algebra
-   `_cell_info_solver.py` - to obtain the RWIGS values for different
    species to 100 % fill the cell
