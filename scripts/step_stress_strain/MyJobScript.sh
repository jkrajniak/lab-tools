#!/bin/bash -l
#PBS -N MELF_Step_deform
#PBS -l walltime=72:00:00
#PBS -o Output.job
#PBS -j oe
#PBS -l nodes=1:ppn=20
#PBS -M jakub.krajniak@cs.kuleuven.be
#PBS -A  lp_polymer_goa_project


module purge
module load GROMACS
module remove impi
module load impi/5.0.1.035
MDRUN="mdrun_mpi"
GROMPP="grompp_mpi"

cd $PBS_O_WORKDIR

# Set up OpenMPI environment
n_proc=$(cat $PBS_NODEFILE | wc -l)
n_node=$(cat $PBS_NODEFILE | uniq | wc -l)
#mpdboot -f $PBS_NODEFILE -n $n_node -r ssh -v

# Prepare deformation part, generate grompp_deform.mdp
bash make_deform.sh

FIRST_STEP_FILES="conf.gro grompp.mdp topol.top topol.itp"
STEP_FILES="topol.top topol.itp grompp.mdp grompp_deform.mdp"

LOG_FILE="${PBS_O_WORKDIR}/output_${PBS_JOBID}.log"

function logg() {
  LOG_FILE="${PBS_O_WORKDIR}/output_${PBS_JOBID}.log"
  echo "==== $1" >> $LOG_FILE
}

# First run for pull_0.00
if [ -d "pull_0.00" ]; then
    echo "Dir pull_0.00 exists"  >> $LOG_FILE
else
    logg "Step pull_0.00"
    mkdir pull_0.00
    for f in $FIRST_STEP_FILES; do
        cp -v $f pull_0.00/ >> $LOG_FILE
    done
    cd pull_0.00/
    grompp_mpi -v  &>> $LOG_FILE
    mpirun -np $n_proc $MDRUN &>> $LOG_FILE
    [ "$?" != "0" ] && exit $?
    cd ..
fi

# Now run rest of the stress-strain pulling
last_step=pull_0.00
for s in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1; do
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

mpdallexit
