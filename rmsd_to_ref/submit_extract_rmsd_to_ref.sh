#!/bin/bash

#SBATCH --job-name=rmsd_to_ref_extract
#SBATCH --output=output_%j.out
#SBATCH --error=error_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ============================================================
# Single job, not an array: extract_gate_latch_rmsd_feats.py already loops
# over every sequence in one seq list in a single Python process.
#
# Usage:
#   sbatch submit_extract_rmsd_to_ref.sh [seq_list] [extra extract_gate_latch_rmsd_feats.py args...]
#
# Examples:
#   sbatch submit_extract_rmsd_to_ref.sh                              # defaults to ../seq_ids_orig.txt
#   sbatch submit_extract_rmsd_to_ref.sh ../seq_ids_orig.txt --tag _500ns
#
# --time above is a rough estimate for 194 sequences; adjust after the
# first run if it under/over-shoots.
# ============================================================

set -euo pipefail

module purge
module load gcc
module load openmpi
module load anaconda
conda activate biosensors

cd "$(dirname "${BASH_SOURCE[0]}")"

SEQ_LIST="${1:-../seq_ids_orig.txt}"
if [[ $# -gt 0 ]]; then
    shift
fi

echo "Running extract_gate_latch_rmsd_feats.py against $SEQ_LIST"
python extract_gate_latch_rmsd_feats.py "$SEQ_LIST" "$@"
