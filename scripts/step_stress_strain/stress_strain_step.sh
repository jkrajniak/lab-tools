#! /bin/bash
#
# sc.sh
# Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

# Start from this step, it has to be prepared beforehand
last_step=pull_0.00
for s in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1; do
    echo "Step $s"
    mkdir pull_${s}
    cd pull_${s}
    cp ../${last_step}/confout.gro conf.gro
    cp ../topol.top .
    cp ../grompp.mdp .
    cp ../grompp_deform.mdp .
    # First run the deformation
    grompp_mpi -f grompp_deform.mdp
    mpiexec -np $n_proc $MDRUN
    # Now run NPT to collect data
    grompp_mpi -f grompp.mdp -c confout.gro
    mpiexec -np $n_proc $MDRUN
    last_step=pull_${s}
    cd ..
    echo "Done step $s"
done
