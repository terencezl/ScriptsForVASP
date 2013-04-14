#!/bin/bash
# Use: In the top working directory
# (the one that generates subfolders like entest),
# ./Display.sh TEST_TYPE min max step
# ./xxtest/xxtest_output is the output file in the form of
# ENCUT  E_LC1   E_LC2   DE=E_LC2-E_LC1   dDE, for entest
# nKP  E_LC1   E_LC2   DE=E_LC2-E_LC1   dDE, for kptest
# N(Ams)  E, for lctest

if [[ $1 == "entest" ||  $1 == "kptest" ]]; then
    fname=$1"_output.txt"
    cd $1                                                           # get into the test dir
    if [ $1 == "entest" ]; then echo "ENCUT(meV)   E_LC1   E_LC2   DE=E_LC2-E_LC1   dDE" >> $fname
    else echo -e "\\nnKP   E_LC1   E_LC2   DE=E_LC2-E_LC1   dDE" >> $fname
    fi
    for ((n=$2; n<=$3; n=n+$4))                                     # display each subfolder's OUTCAR
    do
        E_LC1=$(grep sigma $n/OUTCAR | tail -1 | cut -c 68-)        # get the total energy when LC1
        E_LC2=$(grep sigma $n-1/OUTCAR | tail -1 | cut -c 68-)      # get the total energy when LC2, characterized by patterns like 160-1
        DE=$(echo "($E_LC2)-($E_LC1)" | bc)                         # echo the expr into bc (calculator) to get DE. bash doesn't support floats
        if [ $n == $2 ]; then DE_pre=$DE; fi                        # in the first cycle let dDE = 0
        dDE=$(echo "($DE)-($DE_pre)" | bc)                          # get dDE (the diff of energy difference between two adjacent rows)
#       dDE=$(echo "a=$dDE;if(0>a)a*=-1;a" | bc)                    # return the absolute value of dDE
        dDE=${dDE#-}                                                # return the absolute value of dDE
        echo "$n   $E_LC1   $E_LC2   $DE   $dDE" >> $fname          # output into a file
        DE_pre=$DE                                                  # set DE for this cycle as DE_pre for the next cycle to get the next dDE
    done
    echo -e "\\nMaximal time per run: \c" >> $fname                 # record the maximal time elapsed
    cat < $(find $3 -mindepth 1 -name "*$3*") | grep real | cut -c 6- >> $fname
elif [ $1 == "lctest" ]; then
    fname=$1"_output.txt"
    cd $1                                                           # get into the test dir
    echo -e "\\nLC(Ams)   E" >> $fname
#    List=""                                                         # clear the initial variable
    for n in $(awk "BEGIN{for(i=$2;i<=$3;i+=$4)print i}")           # generate subfolders specified by float numbers
    do
        i=$(echo "scale=2;$n/1" | bc)                               # change decimal format from 5.1 to 5.10
        List=$List" "$i                                             # add List of float numbers up
    done
    for n in $List
    do echo "$n   $(grep sigma $n/OUTCAR | tail -1 | cut -c 68-)" >> $fname; done   # output into a file
    echo -e "\\nMaximal time per run: \c" >> $fname                 # record the maximal time elapsed
    i=$(echo "scale=2;$3/1" | bc)                                   # change the float format to match
elif [ $1 == "rttest" ]; then
    fname=$1"_output.txt"
    cd $1                                                           # get into the test dir
    echo -e "\\nRatio c/a   E" >> $fname
#    List=""                                                         # clear the initial variable
    for n in $(awk "BEGIN{for(i=$2;i<=$3;i+=$4)print i}")           # generate subfolders specified by float numbers
    do
        i=$(echo "scale=2;$n/1" | bc)                               # change decimal format from 5.1 to 5.10
        List=$List" "$i                                             # add List of float numbers up
    done
    for n in $List
    do echo "$n   $(grep sigma $n/OUTCAR | tail -1 | cut -c 68-)" >> $fname; done   # output into a file
    echo -e "\\nMaximal time per run: \c" >> $fname                 # record the maximal time elapsed
    i=$(echo "scale=2;$3/1" | bc)                                   # change the float format to match
    cat < $(find $i -mindepth 1 -name "*$i*") | grep real | cut -c 6- >> $fname   cat < $(find $i -mindepth 1 -name "*$i*") | grep real | cut -c 6- >> $fname
else
    echo "Specify what you are going to test! entest/kptest/lctest"
fi
