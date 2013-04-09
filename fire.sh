#!/bin/bash
# Use: In the top working directory,
# ./fire.sh TEST_TYPE min max step

if [[ $1 == "entest" || $1 == "kptest" ]]; then
    cd $1
    for ((n=$2; n<=$3; n=n+$4))
    do
        qsub $n/qsub.parallel
        qsub $n-1/qsub.parallel
    done
elif [ $1 == "lctest" ]; then
    cd $1
    List=""
    for n in $(awk "BEGIN{for(i=$2;i<=$3;i+=$4)print i}")
    do i=$(echo "scale=2;$n/1" | bc)			# change decimal format from 5.1 to 5.10
    List=$List" "$i
    done
    for n in $List
    do
        qsub $n/qsub.parallel
    done
else
    echo "Specify what you are going to test! entest/kptest/lctest"
fi
