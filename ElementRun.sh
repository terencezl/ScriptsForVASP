#!/bin/bash
# Use: In the TMN working directory
# Change element.dat first

element_list=$(cat INPUT/element.dat)
echo "Tell me what you wish to do: prepare/lctest/display/elastic/solve"
read -r test_type
echo "What lattice type are you thinking of? zincblende/rocksalt/..."
read -r lattice_type

if [ $test_type == prepare ]; then
    for n in $element_list
    do
        mkdir -p $n/$lattice_type/INPUT
        cp -r INPUT/INCAR INPUT/KPOINTS INPUT/POSCAR_$lattice_type INPUT/qsub.parallel $n/$lattice_type/INPUT
        cat INPUT/POTCAR_LDA/POTCAR_$n INPUT/POTCAR_LDA/POTCAR_N > $n/$lattice_type/INPUT/POTCAR
        cd $n/$lattice_type/INPUT
        mv POSCAR_$lattice_type POSCAR
        sed -i "s/@N@/$n N ($lattice_type)/g" INCAR
        sed -i "s/@N@/$n N ($lattice_type)/g" POSCAR
        cd ../../..
    done

elif [ $test_type == lctest ]; then
    echo "Enter the starting scaling factor, ending scaling factor and the step, separated by a space."
    read -r start end step
    for n in $element_list
    do
        cd $n/$lattice_type
        Prepare.sh lctest $start $end $step
        Fire.sh lctest
        cd ../..
    done

elif [ $test_type == display ]; then
    for n in $element_list
    do
        cd $n/$lattice_type
        Display.sh lctest
        cd ../..
    done

elif [ $test_type == elastic ]; then
    echo "Tell me the crystallographic system the compound is in: cubic/tetragonal/orthorhombic/..."
    read -r cryst
    for n in $element_list
    do
        cd $n/$lattice_type
        scaling_factor=$(grep "The scaling factor is" lctest/lctest_output.txt | awk '{print $5}')
        sed -i "s/@R@/$scaling_factor/g" INPUT/POSCAR
        ElasticPrepFire.sh $cryst
        cd ../..
    done

elif [ $test_type == solve ]; then
    echo "Tell me the crystallographic system the compound is in: cubic/tetragonal/orthorhombic/..."
    read -r cryst
    for n in $element_list
    do
        cd $n/$lattice_type
        ElasticDispSolve.sh $cryst
        cd ../..
    done

fi
