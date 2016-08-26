#! /bin/bash -e
#
# sc.sh
# Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

# Number of CPUs in the environment
NPROC=$(cat $PBS_NODEFILE | wc -l)

# COMMANDS
MDRUN="mpirun -n $NPROC mdrun_mpi"
GROMPP="grompp_mpi -v"

FIRST_STEP_FILES="conf.gro grompp.mdp topol.top topol.itp"
STEP_FILES="topol.top topol.itp grompp.mdp grompp_deform.mdp grompp_deform.tpl"

STRAIN_STEPS="0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1"

# Normally you do not need to modify the lines below

function logg() {
  LOG_FILE="output_${PBS_JOBID}.log"
  echo ">>> $1" >> $LOG_FILE
}

FIRST_STEP_FILES="conf.gro grompp.mdp topol.top topol.itp"
STEP_FILES="topol.top topol.itp grompp.mdp grompp_deform.mdp"

LOG_FILE="${PBS_O_WORKDIR}/output_${PBS_JOBID}.log"

function logg() {
  LOG_FILE="${PBS_O_WORKDIR}/output_${PBS_JOBID}.log"
  echo "==== $1" >> $LOG_FILE
}

# First run for pull_0.00
if [ -d "pull_0.00" ] && [ -f "pull_0.00/done" ]; then
    logg "Dir pull_0.00 exists"
else
    if [ -d "pull_0.00" ]; then
        logg "Dir pull_0.00 exists but it is not marked as done, remove it"
        rm -rvf pull_0.00 &>> $LOG_FILE
    fi
    logg "Step pull_0.00"
    mkdir pull_0.00
    for f in $FIRST_STEP_FILES; do
        cp -v $f pull_0.00/ >> $LOG_FILE
    done
    cd pull_0.00/
    $GROMPP -v  &>> $LOG_FILE
    $MDRUN &>> $LOG_FILE
    [ "$?" != "0" ] && exit $?

    touch "done"
    cd ..
fi

# Now run rest of the stress-strain pulling
last_step=pull_0.00
for s in $STRAIN_STEPS; do
    echo "Step $s"  &>> $LOG_FILE
    NEW_STEP_DIR="pull_${s}"
    if [ -d "$NEW_STEP_DIR" ]; then
        if [ -f "$NEW_STEP_DIR/done" ]; then
            echo "Skip step $s" &>> $LOG_FILE
            last_step=$NEW_STEP_DIR
            continue
        else
            logg "Step $s not finished, clean up and run again"
            rm -rvf $NEW_STEP_DIR &>> $LOG_FILE
        fi
    fi
    mkdir "$NEW_STEP_DIR"
    cp -v ${last_step}/confout.gro ${NEW_STEP_DIR}/conf.gro  &>> $LOG_FILE

    for sf in $STEP_FILES; do
        cp -v $sf ${NEW_STEP_DIR}/ &>> $LOG_FILE
    done

    bash make_deform.sh $s
    cd "${NEW_STEP_DIR}"
    # First run the deformation
    logg "=============== DEFORMATION $s ============"

    # Preapre deformation file from the template

    $GROMPP -v -f grompp_deform.mdp &>> $LOG_FILE
    $MDRUN &>> $LOG_FILE
    [ "$?" != "0" ] && exit $?
    logg "============== Collect data $s =============="
    # Now run NPT to collect data
    $GROMPP -f grompp.mdp -c confout.gro -v &>> $LOG_FILE
    $MDRUN &>> $LOG_FILE
    [ "$?" != "0" ] && exit $?

    touch "done"

    last_step=$NEW_STEP_DIR
    cd ..
    logg "================ Finished step $s ================="
done
