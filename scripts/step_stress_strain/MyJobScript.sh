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
MDRUN="mdrun_mpi"
GROMPP="grompp_mpi"

cd $PBS_O_WORKDIR

# Set up OpenMPI environment
n_proc=$(cat $PBS_NODEFILE | wc -l)
n_node=$(cat $PBS_NODEFILE | uniq | wc -l)
mpdboot -f $PBS_NODEFILE -n $n_node -r ssh -v

# Prepare deformation part, generate grompp_deform.mdp
bash make_deform.sh

# First run for pull_0.00
mkdir pull_0.00
cp conf.gro pull_0.00/
cp grompp.mdp pull_0.00/
cp topol.top pull_0.00/
cp topol.itp pull_0.00/
cd pull_0.00/
grompp_mpi -v
mpiexec -np $n_proc $MDRUN
cd ..

# Now run rest of the stress-strain pulling
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
    # Now run NVT to collect data
    grompp_mpi -f grompp.mdp -c confout.gro
    mpiexec -np $n_proc $MDRUN
    last_step=pull_${s}
    cd ..
    echo "Done step $s"
done

mpdallexit
