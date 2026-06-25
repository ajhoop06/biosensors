#!/bin/bash
# submit_contact_analysis.sh
# ---------------------------
# Submits one SLURM job per sequence in seq_ids.txt
#
# Usage:
#   bash submit_contact_analysis.sh

BASE="/scratch/alpine/ivta1597/LCA_boltz_models/LIG_contacts"
SEQ_IDS_FILE="${BASE}/seq_ids.txt"
WORKER="${BASE}/contact_type_worker.sh"

while read -r SEQ_ID _rest || [[ -n "${SEQ_ID}" ]]; do
    [[ -z "${SEQ_ID}" ]] && continue
    echo "Submitting job for: ${SEQ_ID}"
    sbatch --export=SEQ_ID="${SEQ_ID}" "${WORKER}"
done < "${SEQ_IDS_FILE}"
