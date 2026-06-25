#!/bin/bash      

#SBATCH --job-name=water_contact
#SBATCH --output=output_water_%j.out                  # Output file
#SBATCH --error=error_water_%j.err                    # Error file
#SBATCH --account=ucb351_asc3
#SBATCH --partition=amilan
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

module purge
module load anaconda

conda activate IS_env

ID="10"

TOP="/scratch/alpine/ivta1597/LCA_boltz_models/binders/seq${ID}_binder/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/medoid_system.pdb"
XTC="/scratch/alpine/ivta1597/LCA_boltz_models/binders/seq${ID}_binder/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_40_500ns.xtc"
OUTDIR="/scratch/alpine/ivta1597/LCA_boltz_models/binders/seq${ID}_binder/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/water_analysis_results"

RES="55 58 59 61 62 79 81 83 87 88 89 92 94 108 110 115 117 120 122 124 141 143 158 159 160 161 162 163 164 167"

PYTHON=${PYTHON:-python}

$PYTHON water_contact_analysis.py --top "$TOP" --traj "$XTC" --outdir "$OUTDIR" --pocket-residues $RES --stride 10

