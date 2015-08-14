#!/bin/bash

echo "Converting input potential to GROMACS-like .xvg table"

INPUT=$2
NAME="`basename $2`"
TYPE="$1"

csg_call --options settings.xml --ia-type non-bonded --ia-name $TYPE convert_potential gromacs $INPUT new_$NAME.xvg
