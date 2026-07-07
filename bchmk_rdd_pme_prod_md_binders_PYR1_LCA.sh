#!/bin/bash

#SBATCH --job-name=rdd56seq16_prod_md_PYR1_LCA
#SBATCH --output=output_benchmark_%j.out                  # Output file
#SBATCH --error=error_benchmark_%j.err                    # Error file
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

export TMPDIR=$SLURM_SCRATCH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load gcc
module load openmpi
module load anaconda
module load gromacs

conda activate IS_env

# Open MPI network settings (avoid openib on nodes without IB/verbs)
export SLURM_EXPORT_ENV=ALL
# Set some environment variables 
DIR=`pwd`
MDP=$DIR/MDP

# Get sequence value from command line
SEQ=$1
PME=16
RDD=1.2

# Production simulation
cd $DIR/benchmark/${SEQ}/${SLURM_CPUS_PER_TASK}ntomp/load_imbalance/PME
#mkdir prod_md_benchmark_${SLURM_NTASKS}np_${SLURM_CPUS_PER_TASK}omp_${SLURM_CPUS_PER_TASK}cpuspertask_0p95rdd
#cd prod_md_benchmark_${SLURM_NTASKS}np_${SLURM_CPUS_PER_TASK}omp_${SLURM_CPUS_PER_TASK}cpuspertask_0p95rdd
mkdir prod_md_benchmark_${SLURM_NTASKS}np_${SLURM_CPUS_PER_TASK}omp_${SLURM_CPUS_PER_TASK}cpuspertask_1p2rdd_${PME}pme
cd prod_md_benchmark_${SLURM_NTASKS}np_${SLURM_CPUS_PER_TASK}omp_${SLURM_CPUS_PER_TASK}cpuspertask_1p2rdd_${PME}pme
gmx_mpi grompp -f $MDP/prod_md_HMR_3dt_benchmark.mdp -c $DIR/binders/${SEQ}_binder/HMR/NPT/npt.gro -t $DIR/binders/${SEQ}_binder/HMR/NPT/npt.cpt -p $DIR/binders/${SEQ}_binder/HMR/${SEQ}_b_HMR.top -o prod_md.tpr
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md -ntomp $SLURM_CPUS_PER_TASK -rdd $RDD -npme $PME

# Production simulation
#cd $DIR/benchmark/${SLURM_NTASKS}rank
#mkdir prod_md_benchmark_${SLURM_NTASKS}np_${SLURM_CPUS_PER_TASK}omp
#cd prod_md_benchmark_${SLURM_NTASKS}np_${SLURM_CPUS_PER_TASK}omp
#gmx_mpi grompp -f $MDP/prod_md_HMR_benchmark.mdp -c $DIR/binders/${SEQ}_binder/HMR/NPT/npt.gro -t $DIR/binders/${SEQ}_binder/HMR/NPT/npt.cpt -p $DIR/binders/${SEQ}_binder/HMR/${SEQ}_b_HMR.top -o prod_md_100ns.tpr
#mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_100ns -ntomp $SLURM_CPUS_PER_TASK


