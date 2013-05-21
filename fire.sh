#!/bin/bash
# Use: In the top working directory,
# Fire.sh TEST_TYPE

cd $1 || exit 1
dir_list=$(ls -F)
for n in $dir_list 
do
    if [[ "$n" == */ ]]
    then
        (cd $n; qsub qsub.parallel)
    fi
done
