#!/usr/bin/env bash

function argparse {
    while getopts ":d:mfn:" opt; do
        case $opt in
        d)
            directory_name=$OPTARG
            ;;
        m)
            is_submit=true
            ;;
        f)
            is_override=true
            test_tag='-f'
            ;;
        n)
            nband=$OPTARG
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

function does_directory_exist {
    if [[ -d "$directory_name" ]]; then
        cd "$directory_name"
    elif [[ "$current_directory" == "$directory_name"* ]]; then
        echo "You are already in $current_directory/."
        directory_name="$current_directory"
    else
        echo "$directory_name/ does not exist!"
        exit 1
    fi
}

function subdirectory_check {
    if [[ -d "$subdir_name" && $(ls -A "$subdir_name") && -z $is_override ]]; then
        echo "$directory_name/$subdir_name/ contains files. Escaping..."
        exit 1
    fi
}

directory_name=electronic
current_directory=${PWD##*/}
if [[ "$1" == */ ]]; then subdir_name=${1%/}; else subdir_name=$1; fi
test_type="${subdir_name%%_*}"
test_type2=$2
shift 1

if [[ "$test_type" == prepare ]]; then
    argparse "$@"
    mkdir $directory_name 2> /dev/null
    cd $directory_name

    if [[ -d INPUT && $(ls -A INPUT) ]]; then
        echo -n "INPUT/ contains files. "
        if [[ $is_override ]]; then
            echo "Overriding..."
        else
            echo "Escaping..."
            exit 1
        fi
    fi

    cp -r ../INPUT .
#    sed -i "/PREC/c PREC = Accurate" INPUT/INCAR
#    sed -i "/NSW/c NSW = 0" INPUT/INCAR
    sed -i "/LCHARG/c LCHARG = .TRUE." INPUT/INCAR
    sed -i "/LMAXMIX/c LMAXMIX = 4" INPUT/INCAR

    sed -i "/NPAR/c NPAR = 8"  INPUT/INCAR
    sed -i "/#PBS -l walltime/c #PBS -l walltime=03:00:00" INPUT/qsub.parallel
    sed -i "/#PBS -l nodes/c #PBS -l nodes=1:ppn=8" INPUT/qsub.parallel

elif [[ "$test_type" == scrun ]]; then
    argparse "$@"
    does_directory_exist
    subdirectory_check
    Prepare.sh "$subdir_name" $test_tag
    cd scrun
    [[ $is_submit ]] && qsub qsub.parallel

elif [[ "$test_type" == dosrun ]]; then
    argparse "$@"
    does_directory_exist
    subdirectory_check
    Prepare.sh "$subdir_name" $test_tag
    cd dosrun

    if [[ -f CHGCAR ]]; then
        echo "Found CHGCAR under $directory_name/$subdir_name/. Will use it."
        sed -i "/ICHARG/c ICHARG = 11" INCAR
    elif [[ -f ../scrun/CHGCAR && -f ../scrun/CONTCAR ]]; then
        echo "Found CHGCAR and CONTCAR under $directory_name/scrun/. Will use them."
        cp ../scrun/CONTCAR POSCAR
        cp -l ../scrun/CHGCAR .
        sed -i "/ICHARG/c ICHARG = 11" INCAR
    fi

    sed -i '4c 17 17 17' KPOINTS
    sed -i "/NSW/c NSW = 0" INCAR
    sed -i "/ISMEAR/c ISMEAR = -5" INCAR
    sed -i "/NEDOS/c NEDOS = 1501" INCAR

    if [[ $2 == rwigs ]]; then
        rwigs=$(cd ../scrun; CellInfo.sh rwigs |awk '{print $4}')
        sed -i "/RWIGS/c RWIGS = ${rwigs//,/ }" INCAR
        sed -i "/NPAR/c NPAR = 1" INCAR
        sed -i "/LORBIT/c LORBIT = 0"  INCAR
    else
        sed -i "/NPAR/c NPAR = 8"  INCAR
        sed -i "/LORBIT/c LORBIT = 10"  INCAR
    fi

    sed -i "/#PBS -l walltime/c #PBS -l walltime=04:00:00" qsub.parallel
    sed -i "/#PBS -l nodes/c #PBS -l nodes=2:ppn=8" qsub.parallel
    [[ $is_submit ]] && qsub qsub.parallel

elif [[ "$test_type" == bsrun ]]; then
    argparse "$@"
    does_directory_exist
    subdirectory_check
    Prepare.sh "$subdir_name" $test_tag
    cd bsrun
    if [[ ! -f ../scrun/CHGCAR ]]; then
        echo "Didn't find $directory_name/scrun/CHGCAR. Aborted."
        exit 1
    fi
    cp -l ../scrun/CHGCAR .
    if [[ -f ../scrun/CONTCAR ]]; then
        echo "Found CONTCAR under $directory_name/scrun/. Will use it."
        cp ../scrun/CONTCAR POSCAR
    fi

    if [[ ! -f KPOINTS-bs ]]; then
        echo "You must manually change the KPOINTS file before submitting job!"
        exit 1
    fi
    mv KPOINTS-bs KPOINTS
    sed -i "/NSW/c NSW = 0" INCAR
    sed -i "/NEDOS/c NEDOS = 1501" INCAR
    sed -i "/ICHARG/c ICHARG = 11" INCAR
    sed -i "/LORBIT/c LORBIT = 10" INCAR

    sed -i "/NPAR/c NPAR = 8"  INCAR
    sed -i "/#PBS -l walltime/c #PBS -l walltime=04:00:00" qsub.parallel
    sed -i "/#PBS -l nodes/c #PBS -l nodes=2:ppn=8" qsub.parallel
    [[ $is_submit ]] && qsub qsub.parallel

elif [[ "$test_type" == lobster && "$test_type2" == kp ]]; then
    shift 1
    argparse "$@"
    does_directory_exist
    if [[ -d "$subdir_name"-kp && $(ls -A "$subdir_name"-kp) && -z $is_override ]]; then
        echo "$subdir_name-kp/ contains files. Escaping..."
        exit 1
    fi
    Prepare.sh "$subdir_name"-kp $test_tag -a qlobster.kp.serial
    cd "$subdir_name"-kp

    if [[ -f ../scrun/CONTCAR ]]; then
        echo "Found CONTCAR under $directory_name/scrun/. Will use it."
        cp ../scrun/CONTCAR POSCAR
    fi

#    sed -i '4c 11 11 11' KPOINTS
    sed -i "/NSW/c NSW = 0" INCAR
    sed -i "/ISYM/c ISYM = 0" INCAR
    sed -i "/LSORBIT/c LSORBIT = .TRUE." INCAR
    sed -i "/ISMEAR/c ISMEAR = -5" INCAR
    [[ $is_submit ]] && qsub qlobster.kp.serial

elif [[ "$test_type" == lobster && "$test_type2" == test ]]; then
    shift 1
    argparse "$@"
    does_directory_exist
    subdirectory_check
    if [[ -z "$nband" ]]; then
        echo "You must provide NBAND value for the lobster test by -n!"
        exit 1
    fi

    Prepare.sh "$subdir_name" $test_tag -a qlobster.parallel
    cd "$subdir_name"
    if [[ -f KPOINTS-full ]]; then
        echo "Found KPOINTS-full under $directory_name/$subdir_name/. Will use it."
        mv KPOINTS-full KPOINTS
    elif [[ -d lobster-kp ]]; then
        echo "Found lobster-kp/ under $directory_name/$subdir_name/. Will use it."
        cp lobster-kp/IBZKPT KPOINTS
    elif [[ -d ../lobster-kp ]]; then
        echo "Found lobster-kp/ under $directory_name/. Moving to $directory_name/$subdir_name/ for clarity..."
        mv ../lobster-kp .
        cp lobster-kp/IBZKPT KPOINTS
    else
        echo "Didn't find lobster-kp/ or KPOINTS-full. Did you have your own copied here?"
        exit 1
    fi

    if [[ -f CHGCAR ]]; then
        echo "Found CHGCAR under $directory_name/$subdir_name/. Will use it."
        sed -i "/ICHARG/c ICHARG = 11" INCAR
    elif [[ -f ../scrun/CHGCAR && -f ../scrun/CONTCAR ]]; then
        echo "Found CHGCAR and CONTCAR under $directory_name/scrun/. Will use them."
        cp ../scrun/CONTCAR POSCAR
        cp -l ../scrun/CHGCAR .
        sed -i "/ICHARG/c ICHARG = 11" INCAR
    fi

    sed -i "/NBANDS/c NBANDS = $nband" INCAR
    sed -i "/NSW/c NSW = 0" INCAR
    sed -i "/ISMEAR/c ISMEAR = -5" INCAR
    sed -i "/NEDOS/c NEDOS = 1501" INCAR
    sed -i "/LORBIT/c LORBIT = 10"  INCAR
    sed -i "/LWAVE/c LWAVE = .TRUE." INCAR

    sed -i "/NPAR/c NPAR = 8"  INCAR
    sed -i "/#PBS -l walltime/c #PBS -l walltime=04:00:00" qsub.parallel
    sed -i "/#PBS -l nodes/c #PBS -l nodes=2:ppn=8" qsub.parallel
    [[ $is_submit ]] && qsub qsub.parallel

elif [[ "$test_type" == lobster && "$test_type2" == analysis ]]; then
    shift 1
    argparse "$@"
    does_directory_exist
    if [[ -d "$subdir_name" && $(ls -A "$subdir_name") ]]; then
        cd "$subdir_name"
    else
        echo "$subdir_name/ does not exist!"
        exit 1
    fi
    [[ $is_submit ]] && qsub qlobster.parallel

elif [[ "$test_type" == bader && "$test_type2" == test ]]; then
    shift 1
    argparse "$@"
    does_directory_exist
    subdirectory_check
    Prepare.sh "$directory_name" $test_tag -a qbader.serial
    cd "$directory_name"
    sed -i '/LAECHG/c LAECHG = .TRUE.' INCAR
    sed -i '/NGXF/c NGXF = 250' INCAR
    sed -i '/NGYF/c NGYF = 250' INCAR
    sed -i '/NGZF/c NGZF = 250' INCAR
#    echo -e 'NGXF = 250\nNGYF = 250\nNGZF = 250' >> INCAR
    sed -i '/LCHARG/c LCHARG = .TRUE.' INCAR
    sed -i '/NSW/c NSW = 0' INCAR
    [[ $is_submit ]] && qsub qsub.parallel

elif [[ "$test_type" == bader && "$test_type2" == analysis ]]; then
    shift 1
    argparse "$@"
    does_directory_exist
    if [[ -d "$subdir_name" && $(ls -A "$subdir_name") ]]; then
        cd "$subdir_name"
    else
        echo "$subdir_name/ does not exist!"
        exit 1
    fi
    [[ $is_submit ]] && qsub qbader.serial

else
    echo "Specify what you are going to test!" >&2
    exit 1
fi

exit 0
