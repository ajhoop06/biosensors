#!/bin/bash
# submit_contact_analysis.sh
# ---------------------------
# Submits one SLURM job per sequence in seq_list
#
# Usage:
#   bash submit_contact_analysis.sh                                   # seq_ids_orig.txt, full 40–500 ns
#   bash submit_contact_analysis.sh seq_ids_missing_water.txt         # custom list, full 40–500 ns
#   bash submit_contact_analysis.sh seq_ids_orig.txt 40 250           # 250 ns window
#   bash submit_contact_analysis.sh seq_ids_orig.txt 40 300           # 300 ns window
#
# Arguments:
#   $1  seq_list  - path to seq_ids.txt-style sequence list
#                   (default: /projects/ivta1597/biosensors/seq_ids_orig.txt)
#   $2  start_ns  (default: 40)
#   $3  end_ns    (default: 500)

BASE="/projects/ivta1597/biosensors/LIG_contacts"
WORKER="${BASE}/contact_type_worker.sh"

SEQ_IDS_FILE="${1:-/projects/ivta1597/biosensors/seq_ids_orig.txt}"
START_NS="${2:-40}"
END_NS="${3:-500}"

echo "============================================================"
echo "  Contact type analysis submission"
echo "  Seq list : $SEQ_IDS_FILE"
echo "  Window   : ${START_NS}–${END_NS} ns"
echo "  Output   : contact_type_results_${START_NS}_${END_NS}ns/"
echo "============================================================"

submitted=0

while read -r SEQ_ID _rest || [[ -n "${SEQ_ID}" ]]; do
    [[ -z "${SEQ_ID}" || "${SEQ_ID}" == \#* ]] && continue
    echo "Submitting: ${SEQ_ID}  window=${START_NS}–${END_NS}ns"
    sbatch --export=SEQ_ID="${SEQ_ID}",START_NS="${START_NS}",END_NS="${END_NS}" \
           "${WORKER}"
    ((submitted++))
done < "${SEQ_IDS_FILE}"

echo ""
echo "=== Done: submitted ${submitted} jobs ==="
