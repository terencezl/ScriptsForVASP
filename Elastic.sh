#!/usr/bin/env bash
# Use: In the top working directory
# S. K. R. Patil, S. V. Khare, B. R. Tuttle, J. K. Bording, and S. Kodambaka, Mechanical stability of possible structures of PtN investigated using first-principles calculations, PHYSICAL REVIEW B 73, 104118 2006, DOI: 10.1103/PhysRevB.73.104118
# Elastic.sh  start_test / disp-solve / solve  -y cubic / tetragonal /...

if [ $2 == cubic ]; then
    dir_list="c11-c12 c44"
elif [ $2 == cubic_A ]; then
    dir_list="A1 A2 A3 A4 A5 A6"
elif [ $2 == tetragonal ]; then echo
    dir_list="c11 c33 c44 5c11-4c12-2c13+c33 c11+c12-4c13+2c33 c11+c12-4c13+2c33+2c66"
elif [ $2 == orthorhombic ]; then echo
    dir_list="c11 c22 c33 c44 c55 c66 4c11-4c12-4c13+c22+2c23+c33 c11-4c12+2c13+4c22-4c23+c33 c11+2c12-4c13+c22-4c23+4c33"
elif [ $2 == hexagonal ]; then echo
elif [ $2 == trigonal ]; then echo
elif [ $2 == monoclinic ]; then echo
elif [ $2 == triclinic ]; then echo
fi

cd elastic 2> /dev/null
if [ $1 == start_test ]; then
    mkdir elastic && cd elastic
    cp -r ../INPUT .
    sed -i '/NSW/c NSW = 20' INPUT/INCAR
    for n in $dir_list
    do
        VASP_start_test.sh $n -y $2
    done

elif [ $1 == disp-solve ]; then
    for n in $dir_list
    do
        Display.sh $n
    done

    echo $PWD | tee elastic_output.txt
    cd ..
    Elastic.sh solve $2

elif [ $1 == solve ]; then
    echo >> elastic_output.txt
    Vpcell=$(grep 'volume of cell' ../equi-relax/OUTCAR | tail -1 | awk '{print $5;}')
    econst_raw=$(echo "print $(grep 'B0 =' ../lctest/lctest_output.txt | tail -1 | awk '{print $3}')/160.2*$Vpcell" | python)
    for n in $dir_list
    do
        if [ $2 == cubic ]; then
            econst_raw=$econst_raw" "$(grep "Fitting result" $n/$n"_output.txt" | tail -1 | awk '{print $5}')
        elif [ $2 == cubic_A ]; then
            declare ${n}1=$(grep "Fitting result" $n/$n"_output.txt" | tail -1 | awk '{print $8}')
            declare ${n}2=$(grep "Fitting result" $n/$n"_output.txt" | tail -1 | awk '{print $5}')
        fi
    done

#    econst_raw=$econst_raw' '$A11' '$A21' '$A61' '$A12' '$A22' '$A32' '$A42' '$A52' '$A62
    econst_raw=$(echo $econst_raw)
    econst_raw=[${econst_raw// /,}]
    _elastic_solver.py $2 $Vpcell $econst_raw | tee -a elastic_output.txt
fi
