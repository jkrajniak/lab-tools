#!/bin/bash

# Gets the box in z direction
STEPS=5000000
Lz="`tail -n1 conf.gro | sed -e 's/.* \([0-9\.]*$\)/\1/g'`"
echo $Lz
Vz=$(awk "BEGIN {printf \"%.10f\", ${Lz}/${STEPS}*100}")
cat grompp_deform.tpl | sed -e "s/V_DEFORM_Z/${Vz}/g; s/V_STEPS/${STEPS}/g" > grompp_deform.mdp

echo "Lz: ${Lz}"
echo "Vz: ${Vz}"

