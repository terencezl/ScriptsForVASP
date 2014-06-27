#!/usr/bin/env bash
# Use: In the top working directory
# Elastic.sh  start_test / disp-solve / solve  -y cubic / tetragonal /...

test_type="$1"
directory_name=elastic
fname="$directory_name"elastic_output.txt
cryst_sys="$2"
shfit 2

if [[ $cryst_sys == cubic ]]; then
    dir_list="c11-c12 c44"
elif [[ $cryst_sys == tetragonal ]]; then
    dir_list="c11 c33 c44 5c11-4c12-2c13+c33 c11+c12-4c13+2c33 c11+c12-4c13+2c33+2c66"
elif [[ $cryst_sys == orthorhombic ]]; then
    dir_list="c11 c22 c33 c44 c55 c66 4c11-4c12-4c13+c22+2c23+c33 c11-4c12+2c13+4c22-4c23+c33 c11+2c12-4c13+c22-4c23+4c33"
elif [[ $cryst_sys == hexagonal ]]; then echo
elif [[ $cryst_sys == trigonal ]]; then echo
elif [[ $cryst_sys == monoclinic ]]; then echo
elif [[ $cryst_sys == triclinic ]]; then echo
else echo "Expecting 2nd argument to be crystallographic system!"
fi

while getopts ":d:f" opt; do
    case $opt in
    d)
        directory_name=$OPTARG
        echo "-d specifies alternative directory name."
        ;;
    f)
        is_override=true
        echo "-f triggered overriding existing subdirectories."
        ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    :)
        echo "Option -$OPTARG requires an argument." >&2
        exit 1
        ;;
  esac
done

if [[ "$test_type" == start_test ]]; then
    if [[ -d "$directory_name" && $(ls -A "$directory_name") ]]; then
        echo -n "Directory contains files or sub-directories. "
        if [[ $is_override ]]; then
            echo "Overriding..."
        else
            echo "Escaping..."
            exit 1
        fi
    fi

    mkdir "$directory_name" && cd "$directory_name"
    cp -r ../INPUT .
    sed -i '/NSW/c NSW = 20' INPUT/INCAR
    for dir in $dir_list
    do
        Sequence_test.sh $dir $cryst_sys -f
    done


elif [[ "$test_type" == disp-solve ]]; then
    if [[ -d "$directory_name" && $(ls -A "$directory_name") ]]; then
        cd $directory_name
    else
        echo "The directory $directory_name does not exist!"
        exit 1
    fi
    for dir in $dir_list
    do
        Sequence_disp.sh $dir
    done

    echo $PWD | tee $fname
    cd ..
    Elastic.sh solve $cryst_sys


elif [[ "$test_type" == solve ]]; then
    if [[ -d "$directory_name" && $(ls -A "$directory_name") ]]; then
        cd $directory_name
    else
        echo "The directory $directory_name does not exist!"
        exit 1
    fi
    # get the info from elastic constant combination runs and/or lctest. 
    Vpcell=$(grep 'volume of cell' ../equi-relax/OUTCAR | tail -1 | awk '{print $5;}')
    B0=$(grep 'B0 =' ../lctest/lctest_output.txt | tail -1 | awk '{print $3}')
    econst_list=$B0

    for dir in $dir_list
    do
        if [[ $cryst_sys == cubic ]]; then
            econst=$(grep "Fitting result" $dir/$dir"_output.txt" | tail -1 | awk '{print $5}')
            econst=$(python -c "print($econst * 160.2 / $Vpcell)")
            econst_list=$econst_list" "$econst
        fi
    done
    # get rid of the leading or trailing spaces.
    econst_list=$(echo $econst_list)
    econst_list=[${econst_list// /,}]
    echo >> $fname
    _elastic_solver.py $cryst_sys $Vpcell $econst_list | tee -a elastic_output.txt

else
    echo "Specify what you are going to test!" >&2
    exit 1
fi
