#!/bin/bash
# Use: In the TMN working directory
# Change element.dat first

element_list=$(echo $(cat INPUT_ELEMENT/element.dat))
echo "Tell me what you wish to do: prepare / lctest / display / equi-relax / elastic / solve / scrun / dosrun / plot-ldos / plot-tdos / bader-prerun / bader"
read -r test_type

# create INPUT folder in the selected cryst structure of each element
if [ $test_type == prepare ]; then
    for n in $element_list
    do
        mkdir -p $n/INPUT
        cp -r INPUT_ELEMENT/{INCAR,KPOINTS,POSCAR,qsub.parallel} $n/INPUT
        cat ~/terencelz/library/PAW-GGA/POTCAR_{$n,N} > $n/INPUT/POTCAR
        cd $n/INPUT
        sed -i "s/@N@/$n N/g" INCAR
        sed -i "s/@N@/$n N/g" POSCAR
        cd ../..
    done

# basic lctest in the "top working directory"
elif [ $test_type == lctest ]; then
    echo "Enter the starting scaling factor, ending scaling factor and the step, separated by a space."
    read -r start end step
    for n in $element_list
    do
        cd $n || exit 1
        Prepare.sh lctest $start $end $step
        Fire.sh lctest
        cd ..
    done

# fit the lctest curve to Murnaghan eqn of state
elif [ $test_type == display ]; then
    for n in $element_list
    do
        cd $n || exit 1
        Display.sh lctest
        cd ..
    done

# create an equilibrium run
elif [ $test_type == equi-relax ]; then
    for n in $element_list
    do
        cd $n || exit 1
        Prepare.sh equi-relax
        F equi-relax
        cd ..
    done

# prep and fire the elastic constant runs
elif [ $test_type == elastic ]; then
    echo "Tell me the crystallographic system the compound is in: cubic / tetragonal / orthorhombic/..."
    read -r cryst_sys
    for n in $element_list
    do
        cd $n || exit 1
        Elastic.sh $cryst_sys prep-fire
        cd ..
    done

# display the fitting results and solve the elastic constants
elif [ $test_type == solve ]; then
    echo "Tell me the crystallographic system the compound is in: cubic / tetragonal / orthorhombic/..."
    read -r cryst_sys
    for n in $element_list
    do
        cd $n || exit 1
        Elastic.sh $cryst_sys disp-solve
        cd ..
    done

# scrun - do the self-consistent run to get the CHGCAR
# dosrun - do the dos run in the sc-dos-bs folder, putting the last run into a deeper folder called scrun
elif [ $test_type == scrun ] || [ $test_type == dosrun ] || [ $test_type == plot-ldos ] || [ $test_type == plot-tdos ]; then
    mkdir TDOS-at-Ef 2> /dev/null
    for n in $element_list
    do
        cd $n || exit 1
        Dos-bs.sh $test_type
        cat sc-dos-bs/TDOS@Ef-${n}N.txt >> ../../TDOS-at-Ef/TDOS@Ef.txt 2> /dev/null
        cd ..
    done

elif [ $test_type == bader-prerun ]; then
    for n in $element_list
    do
        cd $n || exit 1
        Prep-fast.sh bader
        scaling_factor=$(grep "Equilibrium scaling factor is" lctest/lctest_output.txt | head -1 | awk '{print $5}')
        cd bader
        sed -i "s/@R@/$scaling_factor/g" POSCAR
        echo 'LAECHG = .TRUE.' >> INCAR
        echo -e 'NGXF = 250\nNGYF = 250\nNGZF = 250' >> INCAR
        sed -i 's/LCHARG = .FALSE./#LCHARG = .FALSE./g' INCAR
        qsub qsub.parallel
        cd ../..
    done
    
elif [ $test_type == bader ]; then
    for n in $element_list
    do
        cd $n/bader || exit 1
        chgsum.pl AECCAR0 AECCAR2
        bader CHGCAR -ref CHGCAR_sum
        cd ../..
    done
    for n in $element_list
    do
        grep '' $n/bader/ACF*
    done

fi
