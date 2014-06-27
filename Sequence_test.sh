#!/usr/bin/env bash
# Sequence_test.sh entest -s start -e end -i interval -c scaling_const
# Sequence_test.sh kptest -s start -e end -i interval -c scaling_const
# Sequence_test.sh lctest -s start -e end -n num_points
# Sequence_test.sh rttest -s start -e end -n num_points
# Sequence_test.sh agltest -s start -e end -n num_points -r Ion_rotator_args
# Sequence_test.sh c11-c12 cubic

function create_copy_replace {
    mkdir "$1" 2> /dev/null
    if [[ $? != 0 ]]; then
        if [[ $is_override ]]; then
            echo "    $1 already exists. Overriding input files."
        else
            echo "    $1 already exists. Escaping input files."
            exit 1
        fi
    else
        echo "  Creating $1."
    fi
    cd "$1"
    cp ../../INPUT/INCAR .
    cp ../../INPUT/POSCAR .
    cp ../../INPUT/POTCAR .
    cp ../../INPUT/KPOINTS .
    cp ../../INPUT/qsub.parallel .
    _qsub_replacer.sh qsub.parallel
}

function change_dir_name_with_hyphen {
    for dir in $1
    do
        if [[ $dir == -* ]]; then dir=${dir#-}n; fi
        echo $dir
    done
}

function submission_trigger {
    if [[ $is_submit ]]; then
        echo -e '        \c'
        qsub qsub.parallel
    fi
}

function header_echo {
    echo -e "Creating test directory '$directory_name'..."
    if [[ -d "$directory_name" && "$(ls -A $directory_name)" ]]; then
        echo "  Directory contains files or sub-directories."
    fi
    mkdir "$directory_name" 2> /dev/null
    cd "$directory_name"
    echo "  Preparing "$test_type"..."
    fname=""$test_type""_output.txt
}

function argparse {
    while getopts ":s:e:n:i:c:r:mf" opt; do
        case $opt in
        s)
            start=$OPTARG
            ;;
        e)
            end=$OPTARG
            ;;
        n)
            num_points=$OPTARG
            ;;
        i)
            interval=$OPTARG
            ;;
        c)
            scaling_const=$OPTARG
            ;;
        r)
            ions_rotator_args=$OPTARG
            ;;
        m)
            is_submit=true
            echo "-m triggered job submission."
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
}

if [[ "$1" == */ ]]; then directory_name=${1%/}; else directory_name=$1; fi
test_type="${directory_name%%_*}"
shift 1


if [[ "$test_type" == "entest" || "$test_type" == "kptest" ]]; then
    argparse "$@"
    if [[ -z $start || -z $end || -z $interval || -z $scaling_const ]]; then
        echo "-s -e -i -c should be all set with valid values!" >&2
        exit 1
    fi
    header_echo
    lc1=$scaling_const
    lc2=$(echo "print($lc1+0.1)" | python)
    for ((dir=$start; dir<=$end; dir=$dir+$interval))
    do
        for i in "$dir" "$dir"-1
        do
            (
            create_copy_replace $i
            if [ "$test_type" == "entest" ]; then
                sed -i "s/.*ENCUT.*/ENCUT = $dir/g" INCAR
            else
                sed -i "4c $dir $dir $dir" KPOINTS
            fi
            if [ $i == $dir ]; then
                sed -i "2c $lc1" POSCAR
            else
                sed -i "2c $lc2" POSCAR
            fi
            submission_trigger
            )
        done
    done

elif [[ "$test_type" == "lctest" ]]; then
    argparse "$@"
    if [[ -z $start || -z $end || -z $num_points ]]; then
        echo "-s -e -n should be all set with valid values!" >&2
        exit 1
    fi
    header_echo
    dir_list=$(echo -e "import numpy as np\nfor i in np.linspace($start,$end,$num_points): print('{0:.3f}'.format(i))" | python)
    for dir in $dir_list
    do
        (
        create_copy_replace $dir
        sed -i "2c $dir" POSCAR
        submission_trigger
        )
    done

elif [[ "$test_type" == "rttest" ]]; then
    argparse "$@"
    if [[ -z $start || -z $end || -z $num_points ]]; then
        echo "-s -e -n should be all set with valid values!" >&2
        exit 1
    fi
    header_echo
    dir_list=$(echo -e "import numpy as np\nfor i in np.linspace($start,$end,$num_points): print('{0:.3f}'.format(i))" | python)
    dir_list=$(change_dir_name_with_hyphen "$dir_list")
    for dir in $dir_list
    do
        (
        create_copy_replace $dir
        if [[ "$dir" == *n ]]; then dir=-${dir%n}; fi
        sed -i "s/@R@/$dir/g" POSCAR
        submission_trigger
        )
    done

elif [[ "$test_type" == "agltest" ]]; then
    argparse "$@"
    if [[ -z $start || -z $end || -z $num_points || -z $ions_rotator_args ]]; then
        echo "-s -e -n -r should be all set with valid values!" >&2
        exit 1
    fi
    header_echo
    dir_list=$(echo -e "import numpy as np\nfor i in np.linspace($start,$end,$num_points): print('{0:.3f}'.format(i))" | python)
    dir_list=$(change_dir_name_with_hyphen "$dir_list")
    for dir in $dir_list
    do
        (
        create_copy_replace $dir
        if [[ "$dir" == *n ]]; then dir=-${dir%n}; fi
        Ions_rotator.py $ions_rotator_args -a $dir -p
        submission_trigger
        )
    done

elif [[ "$test_type" == *c[1-9][1-9]* ]]; then
    cryst_sys="$1"
    shift 1
    argparse "$@"
    if [[ -z "$cryst_sys" ]]; then
        echo "crystallographic system should be set!" >&2
        exit 1
    fi
    header_echo
    if [ "$test_type" == c44 ]; then
        dir_list="0.020 0.035 0.050n"
    elif [ "$test_type" == c11-c12 ]; then
        dir_list="0.020 0.030 0.040n"
    fi
    for dir in $dir_list
    do
        (
        create_copy_replace $dir
        if [[ "$dir" == *n ]]; then dir=-${dir%n}; fi
        Strain_applier.py "$test_type" "$cryst_sys" "$dir" -p
        submission_trigger
        )
    done

else
    echo "Specify what you are going to test!" >&2
    exit 1
fi
