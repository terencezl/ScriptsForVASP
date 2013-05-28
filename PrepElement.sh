#!/bin/bash
# Use: In the TMN working directory
# Change element.dat first

element_list=$(cat INPUT/element.dat)
for n in $element_list
do
    mkdir -p $n/$1/INPUT
    cp -r INPUT/INCAR INPUT/KPOINTS INPUT/POSCAR INPUT/qsub.parallel $n/$1/INPUT
    cat INPUT/POTCAR_LDA/POTCAR_$n INPUT/POTCAR_LDA/POTCAR_N > $n/$1/INPUT/POTCAR
    cd $n/$1/INPUT
    sed -i "s/@N@/$n N ($1)/g" INCAR
    sed -i "s/@N@/$n N ($1)/g" POSCAR
    cd ../../..
done
