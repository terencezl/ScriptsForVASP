#!/bin/bash
# Use: In the top working directory,
# Fire.sh TEST_TYPE

directory_name=$1
cd "$directory_name" || exit 1
test_type=${directory_name%%_*}

for n in $(ls -F)
do
    if [[ "$n" == */ ]]; then dir_list=$dir_list" "${n%/}; fi
done

if [[ $test_type == lctest && $2 && $3 && $4 ]]; then
    unset dir_list
    for n in $(awk "BEGIN{for(i=$2;i<=$3;i+=$4)print i}")                   # generate subfolders specified by float numbers
    do
        i=$(echo "scale=3;$n/1" | bc)                                        # change decimal format from 5.1 to 5.10
        dir_list=$dir_list" "$i                                                         # add dir_list of float numbers up
    done
fi

for n in $dir_list 
do
    (cd $n; qsub qsub.parallel)
done
