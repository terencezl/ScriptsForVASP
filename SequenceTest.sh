#!/usr/bin/env bash
# SequenceTest.sh entest -b begin -e end -i interval -c scaling_const
# SequenceTest.sh kptest -b begin -e end -i interval -c scaling_const
# SequenceTest.sh lctest -b begin -e end -n num_points
# SequenceTest.sh rttest -b begin -e end -n num_points
# SequenceTest.sh agltest -b begin -e end -n num_points -r Ion_rotator_args
# SequenceTest.sh c11-c12 cubic

function create_copy_replace {
    if [[ -d "$1" ]]; then
        if [[ $is_override ]]; then
            echo "  $1/ already exists. Overriding input files."
        else
            echo "  $1/ already exists. Escaping input files."
            exit 1
        fi
    else
        echo "  Creating $1/."
    fi
    mkdir "$1" 2> /dev/null
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
        echo -n '    '
        qsub qsub.parallel
    fi
}

function info_echo {
    echo -e "Creating test directory $directory_name/..."
    if [[ -d "$directory_name" && "$(ls -A $directory_name)" ]]; then
        echo "  $directory_name/ contains files or sub-directories."
    fi
    mkdir "$directory_name" 2> /dev/null
    cd "$directory_name"
    echo "  Preparing "$test_type"..."
}

function argparse {
    while getopts ":b:e:n:i:c:r:mf" opt; do
        case $opt in
        b)
            begin=$OPTARG
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
}

if [[ "$1" == */ ]]; then directory_name=${1%/}; else directory_name=$1; fi
test_type="${directory_name%%_*}"
shift 1


if [[ "$test_type" == "entest" || "$test_type" == "kptest" ]]; then
    argparse "$@"
    if [[ -z $begin || -z $end || -z $interval ]]; then
        echo "-b -e -i should be all set with valid values!" >&2
        exit 1
    fi
    info_echo
    for ((dir=$begin; dir<=$end; dir=$dir+$interval))
    do
            (
            create_copy_replace $dir
            if [ "$test_type" == "entest" ]; then
                sed -i "s/.*ENCUT.*/ENCUT = $dir/g" INCAR
            else
                sed -i "4c $dir $dir $dir" KPOINTS
            fi
            [[ -n $scaling_const ]] && sed -i "2c $scaling_const" POSCAR
            submission_trigger
            )
    done

elif [[ "$test_type" == "lctest" ]]; then
    argparse "$@"
    if [[ -z $begin || -z $end || -z $num_points ]]; then
        echo "-b -e -n should be all set with valid values!" >&2
        exit 1
    fi
    info_echo
    dir_list=$(echo -e "import numpy as np\nfor i in np.linspace($begin,$end,$num_points): print('{0:.3f}'.format(i))" | python)
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
    if [[ -z $begin || -z $end || -z $num_points ]]; then
        echo "-b -e -n should be all set with valid values!" >&2
        exit 1
    fi
    info_echo
    dir_list=$(echo -e "import numpy as np\nfor i in np.linspace($begin,$end,$num_points): print('{0:.3f}'.format(i))" | python)
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
    if [[ -z $begin || -z $end || -z $num_points || -z $ions_rotator_args ]]; then
        echo "-b -e -n -r should be all set with valid values!" >&2
        exit 1
    fi
    info_echo
    dir_list=$(echo -e "import numpy as np\nfor i in np.linspace($begin,$end,$num_points): print('{0:.3f}'.format(i))" | python)
    dir_list=$(change_dir_name_with_hyphen "$dir_list")
    for dir in $dir_list
    do
        (
        create_copy_replace $dir
        if [[ "$dir" == *n ]]; then dir=-${dir%n}; fi
        IonsRotator.py $ions_rotator_args -a $dir -p
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
    info_echo

    # create dir_list of each sample point. To save computing resources, oen should manually choose the sample point.
    dir_list=$(echo -e "import numpy as np\nfor i in np.linspace(-0.040,0.040,5): print('{0:.3f}'.format(i))" | python)
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
        StrainApplier.py "$test_type" "$cryst_sys" $dir -p
        submission_trigger
        )
    done

else
    echo "Specify what you are going to test!" >&2
    exit 1
fi
