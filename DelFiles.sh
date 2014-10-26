#!/usr/bin/env bash

while getopts ":twc" opt; do
    case $opt in
    t)
        contcar=true
#        echo "-t also removes CONTCAR."
        ;;
    w)
        wavecar=true
#        echo "-w also removes WAVECAR."
        ;;
    c)
        chgcar=true
#        echo "-c also removes CHGCAR."
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

shift $(($OPTIND-1))
directory_name="$1"

if [[ -d "$directory_name" && $(ls -A "$directory_name") ]]; then
    cd "$directory_name"
elif [[ -z "$directory_name" ]]; then
    echo "Deleting files under the current directory!"
else
    echo "$directory_name/ does not exist!"
    exit 1
fi

rm IBZKPT XDATCAR PCDAT CHG EIGENVAL PROCAR vasprun.xml OUTCAR OSZICAR DOSCAR *.o* 2> /dev/null
[[ $contcar ]] && rm CONTCAR
[[ $wavecar ]] && rm WAVECAR
[[ $chgcar ]] && rm CHGCAR