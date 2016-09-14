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
GROMPP="grompp_mpi"

source ./$1

function logg() {
  LOG_FILE="output_${PBS_JOBID}.log"
  echo ">>> $1" >> $LOG_FILE
}

WORKDIR=${PBS_O_WORKDIR}
if [ "X${PBS_O_WORKDIR}" = "X" ]; then
    WORKDIR="."
fi

LOG_FILE="${WORKDIR}/output_${PBS_JOBID}.log"

function logg() {
  LOG_FILE="${WORKDIR}/output_${PBS_JOBID}.log"
  echo "==== $1" >> $LOG_FILE
}

logg "Running tensil strain experiment, NPROC=$NPROC"

if [ "X$PREFIX" != "X" ]; then
    PREFIX="${PREFIX}_"
fi

# First run for pull_0.00
ZERO_DIR="${PREFIX}pull_0.00"

echo $ZERO_DIR

if [ -d "$ZERO_DIR" ] && [ -f "${ZERO_DIR}/done" ]; then
    logg "Dir pull_0.00 exists"
else
    if [ -d "${ZERO_DIR}" ]; then
        logg "Dir ${ZERO_DIR} exists but it is not marked as done, remove it"
        rm -rvf ${ZERO_DIR} &>> $LOG_FILE
    fi
    logg "Step ${ZERO_DIR}"
    mkdir ${ZERO_DIR}
    for f in $FIRST_STEP_FILES; do
        cp -v $f ${ZERO_DIR}/ &>> $LOG_FILE
    done
    logg "File copied"
    cd ${ZERO_DIR}
    $GROMPP  &>> $LOG_FILE
    $MDRUN &>> $LOG_FILE
    [ "$?" != "0" ] && exit $?

    touch "done"
    cd ..
fi

# Gets the current value of Lz of the box from pull_0.00/confout.gro
Lz="`tail -n1 pull_0.00/confout.gro | sed -e 's/.* \([0-9\.]*$\)/\1/g'`"

logg "Initial box z-size: $Lz"


# Now run rest of the stress-strain pulling
last_step=${ZERO_DIR}
last_strain=0.0

for s in $STRAIN_STEPS; do
    logg $s
    NEW_STEP_DIR="${PREFIX}pull_${s}"
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

    cd "${NEW_STEP_DIR}"

    # First run the deformation
    logg "=============== DEFORMATION $s ============"

    # Prepare deformation file from the template
    Lz="`tail -n1 conf.gro | sed -e 's/.* \([0-9\.]*$\)/\1/g'`"
    ds=$(awk "BEGIN { print ($s-$last_strain)}")
    LzFinal=$(awk "BEGIN { print ($Lz + ${ds}*${Lz})}")

    logg "Deformation by $ds from $Lz to ${LzFinal}"

    make_deform.sh $ds $Lz

    [ "$?" != "0" ] && exit $?

    $GROMPP -f grompp_deform.mdp &>> $LOG_FILE
    $MDRUN &>> $LOG_FILE
    [ "$?" != "0" ] && exit $?
    logg "============== Collect data $s =============="
    # Now run NPT to collect data
    $GROMPP -f grompp.mdp -c confout.gro &>> $LOG_FILE
    $MDRUN &>> $LOG_FILE
    [ "$?" != "0" ] && exit $?

    touch "done"

    last_step=$NEW_STEP_DIR
    last_strain=$s
    cd ..
    logg "================ Finished step $s ================="
done
