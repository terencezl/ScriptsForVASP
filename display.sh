#!/bin/bash
# Use: In the top working directory,
# ./display.sh TEST_TYPE min max step
# ./xxtest/xxtest_output is the output file in the form of
# ENCUT  E_LC1   E_LC2   DE=E_LC2-E_LC1   dDE, for entest
# nKP  E_LC1   E_LC2   DE=E_LC2-E_LC1   dDE, for kptest
# N(Ams)  E, for lctest

if [[ $1 == "entest" ||  $1 == "kptest" ]]; then
    fname=$1"_output.txt"
    cd $1
    if [ $1 == "entest" ]; then echo "ENCUT(meV)   E_LC1   E_LC2   DE=E_LC2-E_LC1   dDE" > $fname
    else echo "nKP   E_LC1   E_LC2   DE=E_LC2-E_LC1   dDE" > $fname
    fi
    for ((n=$2; n<=$3; n=n+$4))
    do
        E_LC1=$(grep sigma $n/OUTCAR | tail -1 | cut -c 68-)				# get the total energy when LC1
        E_LC2=$(grep sigma $n-1/OUTCAR | tail -1 | cut -c 68-)				# when LC2
        DE=$(echo "($E_LC2)-($E_LC1)" | bc)						# echoing the expr into bc to get DE
        if [ $n == $2 ]; then DE_pre=$DE; fi						# the first cycle let dDE = 0
        dDE=$(echo "($DE)-($DE_pre)" | bc)						# get dDE between two adjacent rows
#       if (($(echo $dDE/0.00001 | bc) < 0)); then dDE=$(echo "-($dDE)" | bc); fi	# return the int value to compare with 0, return the abs value of dDE
        dDE=$(echo "a=$dDE;if(0>a)a*=-1;a" | bc)                # return the absolute value of dDE
        echo "$n   $E_LC1   $E_LC2   $DE   $dDE" >> $fname				# output into a file
        DE_pre=$DE									# set DE_pre for the next cycle to get dDE
    done
elif [ $1 == "lctest" ]; then
    fname="lctest_output.txt"
    cd lctest
    echo "LC(Ams)   E" > $fname
    List=""
    for n in $(awk "BEGIN{for(i=$2;i<=$3;i+=$4)print i}")
    do i=$(echo "scale=2;$n/1" | bc)			# change decimal format from 5.1 to 5.10
    List=$List" "$i
    done
    for n in $List
    do echo "$n   $(grep sigma $n/OUTCAR | tail -1 | cut -c 68-)" >> $fname; done	# output into a file
else
    echo "Specify what you are going to test! entest/kptest/lctest"
fi
