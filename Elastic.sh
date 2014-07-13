#!/usr/bin/env bash
# Use: In the top working directory
# Elastic.sh  test / disp-solve / solve   cubic / tetragonal /...

function does_directory_exist {
    if [[ -d "$directory_name" && $(ls -A "$directory_name") ]]; then
        cd "$directory_name"
    elif [[ "$current_directory" == "$directory_name"* ]]; then
        echo "You are already in $current_directory/."
        directory_name="$current_directory"
    else
        echo "$directory_name/ does not exist!"
        exit 1
    fi
}

test_type="$1"
directory_name=elastic
current_directory=${PWD##*/}
fname="$directory_name"_output.txt
cryst_sys="$2"
equi_relax="equi-relax"
lctest="lctest"
shift 2

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
else
    echo "Expecting 2nd argument to be crystallographic system!"
    exit 1
fi

while getopts ":d:q:l:mf" opt; do
    case $opt in
    d)
        directory_name=$OPTARG
        ;;
    q)
        equi_relax=$OPTARG
        ;;
    l)
        lctest=$OPTARG
        ;;
    m)
        is_submit=true
        ;;
    f)
        is_override=true
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

if [[ "$test_type" == test ]]; then
    if [[ -d "$directory_name" && $(ls -A "$directory_name") ]]; then
        echo -n "$directory_name/ contains files or sub-directories. "
        if [[ $is_override ]]; then
            echo "Overriding..."
        else
            echo "Escaping..."
            exit 1
        fi
    fi

    if [[ $is_submit && $is_override ]]; then
        test_tag='-mf'
    elif [[ $is_submit && ! $is_override ]]; then
        test_tag='-m'
    elif [[ ! $is_submit && $is_override ]]; then
        test_tag='-f'
    else
        test_tag=''
    fi

    mkdir "$directory_name" 2> /dev/null
    cd "$directory_name"
    cp -r ../INPUT .
    sed -i '/NSW/c NSW = 20' INPUT/INCAR
    sed -i "/#PBS -l walltime/c #PBS -l walltime=01:00:00" INPUT/qsub.parallel
    for dir in $dir_list
    do
        SequenceTest.sh $dir $cryst_sys $test_tag
    done

    if [[ -d ../$equi_relax ]]; then
        echo "Will not create $directory_name/equi-relax/, but use the $equi_relax/ from outside $directory_name/."
    else
        echo "Did not find $equi_relax/ outside $directory_name. Will create $directory_name/equi-relax/."
        echo -n '  '
        Prepare.sh "equi-relax" $test_tag
    fi

elif [[ "$test_type" == disp-solve ]]; then
    does_directory_exist
    if [[ ! -d "equi-relax" ]]; then
        if [[ -d ../$equi_relax ]]; then
            cp -r ../$equi_relax "equi-relax"
        else
            echo "$directory_name/$equi_relax/ does not exist, neither does $equi_relax/ outside $directory_name/!"
            exit 1
        fi
    fi

    if [[ ! -f "lctest_output.txt" ]]; then
        if [[ -f ../$lctest/lctest_output.txt ]]; then
            cp ../$lctest/lctest_output.txt ./
        else
            echo "$directory_name/lctest_output.txt does not exist and cannot be found in $lctest/ either!"
            exit 1
        fi
    fi

    echo $PWD | tee $fname

    for dir in $dir_list
    do
        SequenceDisp.sh $dir
    done

    cd ..
    Elastic.sh solve $cryst_sys -d "$directory_name"


elif [[ "$test_type" == solve ]]; then
    does_directory_exist
    # get the info from elastic constant combination runs and/or lctest. 
    Vpcell=$(grep 'volume of cell' equi-relax/OUTCAR | tail -1 | awk '{print $5;}')
    B0=$(grep 'B0 =' lctest_output.txt | tail -1 | awk '{print $3}')
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
    _elastic_solver.py $cryst_sys $Vpcell $econst_list | tee -a $fname

else
    echo "Specify what you are going to test!" >&2
    exit 1
fi
