#!/usr/bin/env bash
# Use: In the TMN working directory where there is a list of folders by the name of elements

POT_TYPE=PAW-GGA
POTENTIALS_DIR=$POTENTIALS_DIR

element_list_file=INPUT_ELEMENT/element.dat

while getopts ":e:c:" opt; do
    case $opt in
    e)
        element_list_file=$OPTARG
        ;;
    c)
        pot_combo=$OPTARG
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

element_list=$(cat $element_list_file)

if [[ "$test_type" == prepare ]]; then
    if [[ -z "$pot_combo" ]]; then
        echo "-c potential combination must be specified, comma separated, alternating element as X!"
    for element in $element_list
    do
        mkdir -p "$element"/INPUT
        cp -r INPUT_ELEMENT/* "$element"/INPUT
        rm "$element"/INPUT/*.dat
        pot_combo=${pot_combo//X/"element"}
        cat $POTENTIALS_DIR/$POT_TYPE/POTCAR_{$pot_combo} > "$element"/INPUT/POTCAR
        cd "$element"/INPUT
        sed -i "/SYSTEM/c SYSTEM = $pot_combo" INCAR
        sed -i "1c $pot_combo" POSCAR
        cd ../..
    done

else
    test_script=$1
    shift 1
    if ! type $test_script >/dev/null 2>&1; then
        echo "Command $test_script does not exist!"
        exit 1
    fi
    for element in $element_list
    do
        (
        if ! cd "$element"; then
            echo "Element $element directory does not exist!"
            exit 1
        fi
        $test_script "$@"
        )
    done
fi
