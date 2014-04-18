#!/bin/bash
# Use: In the working directory of elastic consts
# S. K. R. Patil, S. V. Khare, B. R. Tuttle, J. K. Bording, and S. Kodambaka, Mechanical stability of possible structures of PtN investigated using first-principles calculations, PHYSICAL REVIEW B 73, 104118 2006, DOI: 10.1103/PhysRevB.73.104118
# Elastic-solve.sh cubic/tetragonal/... O/A

cd elastic
echo '' >> elastic_output.txt
if [ $2 == O ]; then
    if [ $1 == cubic ]; then
#        dir_list="c11+2c12 c11-c12 c44"
        dir_list="c11-c12 c44"
    elif [ $1 == cubic_A ]; then
        dir_list="A1 A2 A3 A4 A5 A6"
    elif [ $1 == tetragonal ]; then echo
        dir_list="c11 c33 c44 5c11-4c12-2c13+c33 c11+c12-4c13+2c33 c11+c12-4c13+2c33+2c66"
    elif [ $1 == orthorhombic ]; then echo
        dir_list="c11 c22 c33 c44 c55 c66 4c11-4c12-4c13+c22+2c23+c33 c11-4c12+2c13+4c22-4c23+c33 c11+2c12-4c13+c22-4c23+4c33"
    elif [ $1 == hexagonal ]; then echo
    elif [ $1 == trigonal ]; then echo
    elif [ $1 == monoclinic ]; then echo
    elif [ $1 == triclinic ]; then echo
    fi
    echo "Original strain set: $dir_list" | tee -a elastic_output.txt
elif [ $2 == A ]; then
    if [ $1 == cubic ]; then
        dir_list="c11+2c12 c11-c12 2c11+2c12+c44"
    elif [ $1 == tetragonal ]; then echo
        dir_list="c11 c33 c44 5c11-4c12-2c13+c33 c11+c12-4c13+2c33 c11+c12-4c13+2c33+2c66"
    elif [ $1 == orthorhombic ]; then echo
        dir_list="c11 c22 c33 c44 c55 c66 4c11-4c12-4c13+c22+2c23+c33 c11-4c12+2c13+4c22-4c23+c33 c11+2c12-4c13+c22-4c23+4c33"
    elif [ $1 == hexagonal ]; then echo
    elif [ $1 == trigonal ]; then echo
    elif [ $1 == monoclinic ]; then echo
    elif [ $1 == triclinic ]; then echo
    fi
    echo "A slight change of strain set: $dir_list" | tee -a elastic_output.txt
fi

Vpcell=$(grep 'volume of cell' $(cut -d' ' -f1 <<< $dir_list)/0.000/OUTCAR | tail -1 | awk '{print $5;}')
econst_raw=$(echo $(grep 'B0 =' ../lctest/lctest_output.txt | tail -1 | awk '{print $3}')/160.2*$Vpcell)
for n in $dir_list
do
    if [ $1 == cubic ]; then
        econst_raw=$econst_raw" "$(grep "Fitting result" $n/$n"_output.txt" | tail -1 | awk '{print $5}')
    elif [ $1 == cubic_A ]; then
        declare ${n}1=$(grep "Fitting result" $n/$n"_output.txt" | tail -1 | awk '{print $8}')
        declare ${n}2=$(grep "Fitting result" $n/$n"_output.txt" | tail -1 | awk '{print $5}')
    fi
done

econst_raw=$econst_raw' '$A11' '$A21' '$A61' '$A12' '$A22' '$A32' '$A42' '$A52' '$A62
#echo $econst_raw
econst_raw=$(echo $econst_raw)
econst_raw=[${econst_raw// /,}]
_Elastic-solve-solver.py $1 $Vpcell $econst_raw $2 | tee -a elastic_output.txt
