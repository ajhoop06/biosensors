#!/bin/bash

#SBATCH --output=output_benchmark_%j.out                  # Output file
#SBATCH --error=error_benchmark_%j.err                    # Error file
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

export TMPDIR=$SLURM_SCRATCH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SLURM_EXPORT_ENV=ALL

module purge
module load gcc
module load openmpi
module load anaconda
module load gromacs

conda activate IS_env

# Set some environment variables 
DIR=`pwd`
MDP=$DIR/MDP

# Get sequence value from command line
SEQ=$1
PME=16
RDD=1.2

D1=6
D2=4
D3=2

# Production simulation
cd $DIR/benchmark/${SEQ}/${SLURM_CPUS_PER_TASK}ntomp/dodecahedron/100ps/DD
#mkdir prod_md_benchmark_${SLURM_NTASKS}np_${SLURM_CPUS_PER_TASK}omp_${SLURM_CPUS_PER_TASK}cpuspertask_100ps_rep3
#cd prod_md_benchmark_${SLURM_NTASKS}np_${SLURM_CPUS_PER_TASK}omp_${SLURM_CPUS_PER_TASK}cpuspertask_100ps_rep3
#mkdir prod_md_benchmark_${SLURM_NTASKS}np_${SLURM_CPUS_PER_TASK}omp_${SLURM_CPUS_PER_TASK}cpuspertask_100ps_${PME}pme_rep3
#cd prod_md_benchmark_${SLURM_NTASKS}np_${SLURM_CPUS_PER_TASK}omp_${SLURM_CPUS_PER_TASK}cpuspertask_100ps_${PME}pme_rep3
mkdir prod_md_benchmark_${SLURM_NTASKS}np_${SLURM_CPUS_PER_TASK}omp_${SLURM_CPUS_PER_TASK}cpuspertask_${PME}pme_${D1}${D2}${D3}dd_rep5
cd prod_md_benchmark_${SLURM_NTASKS}np_${SLURM_CPUS_PER_TASK}omp_${SLURM_CPUS_PER_TASK}cpuspertask_${PME}pme_${D1}${D2}${D3}dd_rep5

# Binders
# dodecahedron
#gmx_mpi grompp -f $MDP/prod_md_HMR_3dt_benchmark.mdp -c $DIR/binders/${SEQ}_binder/HMR/dodecahedron/NPT/npt.gro -t $DIR/binders/${SEQ}_binder/HMR/dodecahedron/NPT/npt.cpt -p $DIR/binders/${SEQ}_binder/HMR/dodecahedron/${SEQ}_b_dodecahedron_HMR.top -o prod_md.tpr

# rectangular
#gmx_mpi grompp -f $MDP/prod_md_HMR_3dt_benchmark.mdp -c $DIR/binders/${SEQ}_binder/HMR/NPT/npt.gro -t $DIR/binders/${SEQ}_binder/HMR/NPT/npt.cpt -p $DIR/binders/${SEQ}_binder/HMR/${SEQ}_b_HMR.top -o prod_md.tpr

# Nonbinders
# dodecahedron
gmx_mpi grompp -f $MDP/prod_md_HMR_3dt_benchmark.mdp -c $DIR/nonbinders/${SEQ}_nb/HMR/dodecahedron/NPT/npt.gro -t $DIR/nonbinders/${SEQ}_nb/HMR/dodecahedron/NPT/npt.cpt -p $DIR/nonbinders/${SEQ}_nb/HMR/dodecahedron/${SEQ}_nb_dodecahedron_HMR.top -o prod_md.tpr

# rectangular
#gmx_mpi grompp -f $MDP/prod_md_HMR_3dt_benchmark.mdp -c $DIR/nonbinders/${SEQ}_nb/HMR/NPT/npt.gro -t $DIR/nonbinders/${SEQ}_nb/HMR/NPT/npt.cpt -p $DIR/nonbinders/${SEQ}_nb/HMR/${SEQ}_nb_HMR.top -o prod_md.tpr

mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md -ntomp $SLURM_CPUS_PER_TASK -npme $PME -dd $D1 $D2 $D3 #-rdd $RDD

# Production simulation
#cd $DIR/benchmark/${SLURM_NTASKS}rank
#mkdir prod_md_benchmark_${SLURM_NTASKS}np_${SLURM_CPUS_PER_TASK}omp
#cd prod_md_benchmark_${SLURM_NTASKS}np_${SLURM_CPUS_PER_TASK}omp
#gmx_mpi grompp -f $MDP/prod_md_HMR_benchmark.mdp -c $DIR/binders/${SEQ}_binder/HMR/NPT/npt.gro -t $DIR/binders/${SEQ}_binder/HMR/NPT/npt.cpt -p $DIR/binders/${SEQ}_binder/HMR/${SEQ}_b_HMR.top -o prod_md_100ns.tpr
#mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_100ns -ntomp $SLURM_CPUS_PER_TASK


