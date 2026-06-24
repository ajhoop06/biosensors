#!/bin/bash

#SBATCH --job-name=eq_PYR1_LCA
#SBATCH --output=output_%j.out                  # Output file
#SBATCH --error=error_%j.err                    # Error file
#SBATCH --account=ucb351_asc4
#SBATCH --partition=amilan
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

export TMPDIR=$SLURM_SCRATCH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load gcc
module load openmpi
module load anaconda
module load gromacs

# Set some environment variables
DIR=/projects/ivta1597/biosensors
MDP=$DIR/MDP

# Get sequence value from command line
ID=$1
SEQ_TYPE=$2 # binders | nonbinders | neg_fail_gate | neg_low_pkt
PREFIX=$3 # pair, bind

if [ "$SEQ_TYPE" == "binders" ]; then
    SUFFIX="binder"
elif [ "$SEQ_TYPE" == "nonbinders" ]; then
    SUFFIX="nb"
elif [ "$SEQ_TYPE" == "neg_fail_gate" ]; then
    SUFFIX="fail_gate"
elif [ "$SEQ_TYPE" == "neg_low_pkt" ]; then
    SUFFIX="low_pkt"
else
    echo "ERROR: Unknown SEQ_TYPE '$SEQ_TYPE'" >&2
    exit 1
fi

# dodecahedron unit cell
# NVT
cd $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}
mkdir NVT
cd NVT
gmx grompp -f $MDP/nvt.mdp -c $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/EM/em.gro -r $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/EM/em.gro -p $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/${PREFIX}_${ID}_dodecahedron_HMR.top -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT
cd $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}
mkdir NPT
cd NPT
gmx grompp -f $MDP/npt.mdp -c $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/NVT/nvt.gro -t $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/NVT/nvt.cpt -p $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/${PREFIX}_${ID}_dodecahedron_HMR.top -r $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/NVT/nvt.gro -o npt.tpr
gmx mdrun -deffnm npt

