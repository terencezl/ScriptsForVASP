#!/bin/bash
# Rerun.sh

cd $1 || exit 1
mkdir relax
mv * relax 2> /dev/null
cp relax/{INCAR,KPOINTS,POTCAR,CONTCAR,CHGCAR,qsub.parallel} .
rm -r relax
mv CONTCAR POSCAR
#qsub qsub.parallel
