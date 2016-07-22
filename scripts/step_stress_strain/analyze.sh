#! /bin/sh
#
# analyze.sh
# Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

OUT="data.csv"
for p in 0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1; do
    echo $p
    cd pull_${p}
    echo "37 0 " | g_energy_mpi
    cd ..
    DATA="`g_analyze_mpi -f pull_${p}/energy.xvg | grep SS1 | sed -e 's/SS1//g'`"
    echo "${p} ${DATA}" >> $OUT
done
