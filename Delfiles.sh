#!/bin/bash

mkdir Temp
mv INCAR KPOINTS POTCAR POSCAR qsub.parallel CHGCAR Temp/
rm * 2>/dev/null
mv Temp/* .
rm -r Temp
