#!/bin/bash
# Use: In the top working directory,
# Display.sh TEST_TYPE (e.g. c11-c12)

if [[ "$1" == */ ]]; then test_type=${1%/}; else test_type=$1; fi
cd $test_type || exit 1
fname=$test_type"_output.txt"
ls_list=$(ls -F)
for n in $ls_list
do
    if [[ "$n" == */ ]]; then dir_list="$dir_list ${n%/}"; fi
done

if [[ $test_type == "entest" || $test_type == "kptest" ]]; then
    if [ $test_type == "entest" ]; then echo -e "\nENCUT(meV)\tE_LC1\tE_LC2\tDE=E_LC2-E_LC1\tdDE" >> $fname
    else echo -e "\nnKP\tE_LC1\tE_LC2\tDE=E_LC2-E_LC1\tdDE" >> $fname
    fi
    for n in $dir_list
    do
        E_LC1=$(grep sigma $n/OUTCAR | tail -1 | awk '{print $7;}')        # get the total energy when LC1
        E_LC2=$(grep sigma $n-1/OUTCAR | tail -1 | awk '{print $7;}')      # get the total energy when LC2, characterized by patterns like 160-1
        DE=$(echo "($E_LC2)-($E_LC1)" | bc)                         # echo the expr into bc (calculator) to get DE. bash doesn't support floats
        smallest=$(echo $dir_list | awk '{print $test_type}')
        if [ $n == $smallest ]; then DE_pre=$DE; fi                        # in the first cycle let dDE = 0
        dDE=$(echo "($DE)-($DE_pre)" | bc)                          # get dDE (the diff of energy difference between two adjacent rows)
        dDE=${dDE#-}                                                # return the absolute value of dDE
        echo -e "$n\t$E_LC1\t$E_LC2\t$DE\t$dDE" >> $fname          # output into a file
        DE_pre=$DE                                                  # set DE for this cycle as DE_pre for the next cycle to get the next dDE
    done

elif [[ $test_type == "lctest" || $test_type == "rttest" || $test_type == "mesh2d" ]]; then
    if [ $test_type == "lctest" ]; then echo -e "\nScaling Constant (Ams)\tE" >> $fname
    elif [ $test_type == "rttest" ]; then echo -e "\nRatio\tE" >> $fname
    else echo -e "\nScaling Constant (Ams)\tRatio\tE" >> $fname
    fi
    for n in $dir_list
    do
        echo -e "$n\t$(grep sigma $n/OUTCAR | tail -1 | awk '{print $7;}')" >> $fname
    done

elif [[ $test_type == "c11-c12" || $test_type == "c11+2c12" || $test_type == "c44" ]]; then
    echo -e "\nDelta (ratio)\tE" >> $fname
    for n in $dir_list
    do
        if [[ "$n" == *n ]]; then i=-${n%n}; else i=$n; fi
        echo -e "$i\t$(grep sigma $n/OUTCAR | tail -1 | awk '{print $7;}')" >> $fname
    done
fi

echo -e "\nMaximal time per run: \c" >> $fname
largest=$(echo $dir_list | awk '{print $NF}')
cat < $(find $largest -mindepth 1 -name "*$largest*") | grep real | awk '{print $2;}' >> $fname
