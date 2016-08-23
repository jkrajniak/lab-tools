#! /bin/bash
#
# sc.sh
# Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

FIRST_STEP_FILES="conf.gro grompp.mdp topol.top topol.itp"
STEP_FILES="topol.top topol.itp grompp.mdp grompp_deform.mdp"
INITIAL_COORDINATE="conf.gro"

STRAIN_STEPS="0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.9 0.1"

function logg() {
   LOG_FILE="output_${PBS_JOBID}.log"
  echo ">>> $1" >> $LOG_FILE
}

# Prepare deformation part, generate grompp_deform.mdp
bash make_deform.sh

# First run for pull_0.00
if [ -d "pull_0.00" ] && [ -f "pull_0.00/done" ]; then
    logg "Dir pull_0.00 exists, skip it"
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

    touch done
    cd ..
fi

# Now run rest of the stress-strain pulling
last_step=pull_0.00
for s in $STRAIN_STEPS; do
    logg "Step $s"
    NEW_STEP_DIR="pull_${s}"
    if [ -d "$NEW_STEP_DIR" ]; then
        if [ -f "$NEW_STEP_DIR/done" ]; then
            logg "Skip step $s"
            last_step=$NEW_STEP_DIR
            continue
        else
            logg "Step $s not finished, clean up and run again"
            rm -rvf $NEW_STEP_DIR &>> $LOG_FILE
        fi
    fi
    mkdir "$NEW_STEP_DIR"
    cp -v ${last_step}/confout.gro ${NEW_STEP_DIR}/conf.gro

    for sf in $STEP_FILES; do
        cp -v $sf ${NEW_STEP_DIR}/
    done
    cd "${NEW_STEP_DIR}"
    # First run the deformation
    logg "=============== DEFORMATION ============"
    grompp_mpi -v -f grompp_deform.mdp &>> $LOG_FILE
    mpirun -np $n_proc $MDRUN &>> $LOG_FILE
    [ "$?" != "0" ] && exit $?
    logg "============== Collect data =============="
    # Now run NPT to collect data
    grompp_mpi -f grompp.mdp -c confout.gro -v &>> $LOG_FILE
    mpirun -np $n_proc $MDRUN &>> $LOG_FILE
    [ "$?" != "0" ] && exit $?

    touch "done"

    last_step=$NEW_STEP_DIR
    cd ..
    logg "================ Finished step $s ================="
done
