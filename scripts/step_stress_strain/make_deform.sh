#!/bin/bash

# Deformation step $1
STRAIN=$1
TEMPLATE_FILE=grompp_deform.tpl

# Gets the box in z direction
STEPS=$(cat $TEMPLATE_FILE | grep nsteps | cut -f2 -d'=')
DT=$(cat $TEMPLATE_FILE | grep dt | cut -f2 -d'=')
echo Deformation steps $STEPS dt=$DT
Lz="`tail -n1 conf.gro | sed -e 's/.* \([0-9\.]*$\)/\1/g'`"
Vz=$(awk "BEGIN {printf \"%.16f\", ${STRAIN}*${Lz}/(${STEPS}*${DT})}")
cat grompp_deform.tpl | sed -e "s/V_DEFORM_Z/${Vz}/g;" > grompp_deform.mdp

echo "Lz: ${Lz}"
echo "Vz: ${Vz}"

