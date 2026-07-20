#!/bin/bash

#SBATCH --job-name=hydration_calc
#SBATCH --output=output_hydration_%j.out
#SBATCH --error=error_hydration_%j.err
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
# Usage:
#   sbatch run_hydration_analysis.sh <seq_id> <seq_type> [start_ns] [end_ns] [hydration_cutoff] [stride]
#
# Arguments:
#   seq_id            - sequence identifier (e.g. pair_3059_binder)
#   seq_type          - directory group (binders | nonbinders | neg_low_pkt | neg_fail_gate)
#   start_ns          - start of analysis window in ns (default: 40)
#   end_ns            - end of analysis window in ns   (default: 500)
#   hydration_cutoff  - first-hydration-shell cutoff in Angstrom (default: 3.5)
#   stride            - frame stride (default: 10)
#
# Example:
#   sbatch run_hydration_analysis.sh pair_3059_binder binders
# ============================================================

set -euo pipefail

module purge
module load anaconda
conda activate biosensors

seq_id=$1
seq_type=$2
start_ns=${3:-40}
end_ns=${4:-500}
hydration_cutoff=${5:-3.5}
stride=${6:-10}

echo "============================================================"
echo "  Hydration-shell water count"
echo "  seq_id   : $seq_id"
echo "  seq_type : $seq_type"
echo "  window   : ${start_ns}-${end_ns} ns"
echo "  cutoff   : ${hydration_cutoff} A"
echo "  stride   : ${stride}"
echo "  start    : $(date)"
echo "============================================================"

python hydration_calc.py \
    --seq_id   "$seq_id"   \
    --seq_type "$seq_type" \
    --start-ns "$start_ns" \
    --end-ns   "$end_ns"   \
    --hydration-cutoff "$hydration_cutoff" \
    --stride   "$stride"

echo ""
echo "============================================================"
echo "  Done: $seq_id"
echo "  End: $(date)"
echo "============================================================"
