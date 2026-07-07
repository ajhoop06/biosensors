#!/bin/bash
# contact_type_worker.sh
# ----------------------
# SLURM worker script for a single sequence.
# SEQ_ID is passed in via --export by submit_contact_analysis.sh

#SBATCH --job-name=contact_type
#SBATCH --output=logs/contact_type_%j.out
#SBATCH --error=logs/contact_type_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ─────────────────────────────────────────────
# USER VARIABLES
# ─────────────────────────────────────────────
PYTHON_SCRIPT="/scratch/alpine/ivta1597/LCA_boltz_models/LIG_contacts/contact_type_analysis.py"
CONDA_ENV="IT_env"
LOG_DIR="/scratch/alpine/ivta1597/LCA_boltz_models/LIG_contacts/logs"

# ─────────────────────────────────────────────
set -euo pipefail
mkdir -p "${LOG_DIR}"

if [[ -z "${SEQ_ID}" ]]; then
    echo "ERROR: SEQ_ID is not set. Submit via submit_contact_analysis.sh" >&2
    exit 1
fi

echo "──────────────────────────────────────────"
echo "Job ID     : ${SLURM_JOB_ID}"
echo "Seq ID     : ${SEQ_ID}"
echo "Node       : $(hostname)"
echo "Start time : $(date)"
echo "──────────────────────────────────────────"

module purge
module load anaconda
conda activate "${CONDA_ENV}"

python "${PYTHON_SCRIPT}" "${SEQ_ID}"

echo "Finished at: $(date)"
