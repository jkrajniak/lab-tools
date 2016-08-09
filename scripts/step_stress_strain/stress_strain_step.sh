#! /bin/bash
#
# sc.sh
# Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

FIRST_STEP_FILES="conf.gro grompp.mdp topol.top topol.itp"
STEP_FILES="topol.top topol.itp grompp.mdp grompp_deform.mdp"

# Prepare deformation part, generate grompp_deform.mdp
bash make_deform.sh

# First run for pull_0.00
if [ -d "pull_0.00" ]; then
    echo "Dir pull_0.00 exists"
else
    mkdir pull_0.00
    for f in $FIRST_STEP_FILES; do
        cp -v $f pull_0.00/
    done
    cd pull_0.00/
    grompp_mpi -v
    [ "$?" != "0" ] && exit $?
    mpiexec -np $n_proc $MDRUN
    [ "$?" != "0" ] && exit $?
    cd ..
fi

# Now run rest of the stress-strain pulling
last_step=pull_0.00
for s in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1; do
    echo "Step $s"
    NEW_STEP_DIR="pull_${s}"
    if [ -d "$NEW_STEP_DIR" ]; then
        echo "Skip step $s"
        continue
    fi
    mkdir "$NEW_STEP_DIR"
    cp -v ${last_step}/confout.gro ${NEW_STEP_DIR}/conf.gro

    for sf in $STEP_FILES; do
        cp -v $sf ${NEW_STEP_DIR}/
    done
    cd "${NEW_STEP_DIR}"
    # First run the deformation
    grompp_mpi -f grompp_deform.mdp
    [ "$?" != "0" ] && exit $?
    mpiexec -np $n_proc $MDRUN
    [ "$?" != "0" ] && exit $?
    # Now run NPT to collect data
    grompp_mpi -f grompp.mdp -c confout.gro
    [ "$?" != "0" ] && exit $?
    mpiexec -np $n_proc $MDRUN
    [ "$?" != "0" ] && exit $?
    last_step=$NEW_STEP_DIR
    cd ..
    echo "Done step $s"
done
