#!/bin/bash
# submit_hydration_analysis.sh
# ─────────────────────────────────────────────────────────────────────────────
# Reads seq_ids.txt and submits a separate run_hydration_analysis.sh SLURM job
# for each sequence. Sequences with a custom path (3rd column) are skipped
# since their trajectories live in a non-standard directory — run those
# manually with: sbatch run_hydration_analysis.sh <seq_id> <dir_type> [start_ns] [end_ns]
#
# Usage:
#   bash submit_hydration_analysis.sh                       # 40–500 ns (default)
#   bash submit_hydration_analysis.sh seq_ids.txt           # specify seq list
#   bash submit_hydration_analysis.sh seq_ids.txt 40 250    # 250 ns window
#   bash submit_hydration_analysis.sh seq_ids.txt 40 500 3.5 10   # explicit cutoff/stride
#
# seq_ids.txt format (tab-separated):
#   seq_id              seq_type (display)    optional_custom_path
#   pair_3069_binder    Binder
#   seq14_binder        Binder                /scratch/.../HMR/dodecahedron   ← skipped
# ─────────────────────────────────────────────────────────────────────────────

SEQ_LIST=${1:-/projects/ivta1597/biosensors/seq_ids_ngs_observed.txt}
START_NS=${2:-40}
END_NS=${3:-500}
HYDRATION_CUTOFF=${4:-3.5}
STRIDE=${5:-10}

if [ ! -f "$SEQ_LIST" ]; then
    echo "ERROR: seq list file not found: $SEQ_LIST"
    exit 1
fi

echo "============================================================"
echo "  Hydration-shell analysis submission"
echo "  Seq list : $SEQ_LIST"
echo "  Window   : ${START_NS}–${END_NS} ns"
echo "  Cutoff   : ${HYDRATION_CUTOFF} A"
echo "  Stride   : ${STRIDE}"
echo "  Out dir  : water_density_${START_NS}_${END_NS}ns/"
echo "============================================================"

# ── Map display seq_type → directory name used in the file system ─────────────
get_dir_type() {
    case "$1" in
        "Binder")         echo "binders"      ;;
        "False Positive") echo "nonbinders"   ;;
        "Low Confidence") echo "neg_low_pkt"  ;;
        "Fail Geometry")  echo "neg_fail_gate";;
        *)                echo "$1"           ;;   # fallback: use as-is
    esac
}

submitted=0
skipped=0

while IFS=$'\t' read -r seq_id seq_type custom_path || [[ -n "$seq_id" ]]; do

    # Skip empty lines and comments
    [[ -z "$seq_id" || "$seq_id" == \#* ]] && continue

    # Skip sequences with a custom path — their trajectories are in a
    # non-standard location that run_hydration_analysis.sh does not handle
    if [[ -n "$custom_path" ]]; then
        echo "SKIP (custom path): $seq_id"
        ((skipped++))
        continue
    fi

    dir_type=$(get_dir_type "$seq_type")

    echo "Submitting: $seq_id  [$seq_type → $dir_type]  window=${START_NS}–${END_NS}ns"
    sbatch run_hydration_analysis.sh "$seq_id" "$dir_type" "$START_NS" "$END_NS" "$HYDRATION_CUTOFF" "$STRIDE"
    ((submitted++))

done < "$SEQ_LIST"

echo ""
echo "=== Done ==="
echo "  Submitted : $submitted jobs"
echo "  Skipped   : $skipped sequences (run manually)"
echo ""
echo "  To run skipped sequences manually:"
echo "  sbatch run_hydration_analysis.sh <seq_id> <dir_type> $START_NS $END_NS $HYDRATION_CUTOFF $STRIDE"
