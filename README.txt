Main Task Scripts:

ElasticPrepFire.sh       ElasticDispSolve.sh
       /      \              /           \
      /        \            /        ElasticSolve.sh
     /          \          / 
Prepare.sh    Fire.sh    Display.sh
                     |||
                  _files.py

----------------------------------------------------

Auxiliary Scripts:

Cellinfo.sh - Go into the deepest VASP working directory and execute
    |
 _flies.py

----------------------------------------------------

Small Quick Scripts:

PrepFast.sh
Delfiles.sh - Go into the deepest VASP working directory and execute
PrepRerun.sh
PrepElement.sh

####################################################################################

General Uses:

1. Create a directory for a specific compound, say PtN. This is the top working directory
2. Go into the directory and create folder INPUT, putting INCAR, KPOINTS, POSCAR, concatenated POTCAR and qsub.parallel into it
3. Customize the input files to prepare the run, including reserving the strings that need to be replaced as @R@, @N@, etc.
4. Go back to the top working directory. Use Prepare.sh xxx yyy to perform a set of runs, or PrepFast.sh xxx to simply create a folder xxx with proper substitutions
5. Use Fire.sh xxx to execute a set of qsub batch commands, or simply go into the folder prepared by PrepFast.sh and qsub qsub.parallel
6. Using Display.sh xxx to get the (step, energy) pairs written into xxx_output.txt and other instructive information. For lctest and elastic constant runs the fitting results and plot figures are also written
7. Rename the folder created during the process under the top working directory to whatever looks neat and tidy, i.e. labeling them with numbers

----------------------------------------------------

To perform a complete set of elastic constant runs:

1. Under the top working directory, execute ElasticPrepFire.sh xxx, i.e. cubic
2. ElasticDispSolve.sh xxx after every batched job is finished
3. Go check the elastic/output.txt for the results

----------------------------------------------------

To use a slightly different set of elastic constant equations (i.e. changing only one equation):

1. Perform an original complete run with ElasticPrepFire.sh, ElasticDispSolve.sh
2. Perform an additional run (with a specific combination of elastic consts) with Prepare.sh, Fire.sh and Display.sh
3. Change the dir_list in ElasticSolve.sh and one line of the matrix in _ElasticSolve_solver.py
4. Re-execute ElasticSolve.sh, the results will show on the screen. Save or redirect it if you wish
