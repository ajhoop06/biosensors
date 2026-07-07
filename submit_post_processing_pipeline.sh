#!/bin/bash
# Reads seq_ids.txt and config.yaml, then submits one post_processing_pipeline_worker.sh
# SLURM job per sequence. Each job runs the full pipeline (PBC/fitting, medoid
# search, and all medoid-referenced analyses) for that sequence directory.
#
# Usage:
#   bash submit_post_processing_pipeline.sh                                    # seq_ids.txt, config.yaml, all phases, no force
#   bash submit_post_processing_pipeline.sh my_list.txt                        # custom sequence list
#   bash submit_post_processing_pipeline.sh my_list.txt my_cfg.yaml            # custom list and config
#   bash submit_post_processing_pipeline.sh my_list.txt my_cfg.yaml all true   # force-overwrite all existing outputs
#   bash submit_post_processing_pipeline.sh seq_ids.txt config.yaml 1          # phase 1 only
#   PHASE: 1, 2, 3, or all (default: all)
#   FORCE: true or false (default: false)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SEQ_LIST="${1:-${SCRIPT_DIR}/seq_ids.txt}"
CONFIG="${2:-${SCRIPT_DIR}/config.yaml}"
PHASE="${3:-all}"
FORCE="${4:-false}"
WORKER="${SCRIPT_DIR}/post_processing_pipeline_worker.sh"

for f in "$SEQ_LIST" "$CONFIG" "$WORKER"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f" >&2; exit 1; }
done

# ── HPC environment (needed for python3 + PyYAML below) ───────────────────────
# `conda activate` requires the shell function from `conda shell.bash hook`;
# this script runs as a plain (non-login) bash process via `bash script.sh`,
# so ~/.bashrc's `conda init` block is never sourced and the bare `conda`
# binary doesn't support activate/deactivate at all without it.
module purge
module load gcc
module load openmpi
module load anaconda
eval "$(conda shell.bash hook)"
conda activate biosensors

# ── Read base paths and seq_type → subdirectory map from config.yaml ──────────
eval "$(python3 << PYEOF
import yaml
with open("${CONFIG}") as f:
    d = yaml.safe_load(f)
print('BASE='   + repr(d['paths']['base']))
print('RUNREL=' + repr(d['paths']['runrel']))
labels  = d.get('seq_labels', {})
entries = ' '.join(f'[{repr(k)}]={repr(v)}' for k, v in labels.items())
print(f'declare -A LABEL_MAP=( {entries} )')
PYEOF
)"

submitted=0
skipped=0
missing=0

while IFS=$'\t' read -r seq_id seq_type custom_path || [[ -n "$seq_id" ]]; do
    [[ -z "$seq_id" || "$seq_id" == \#* ]] && continue

    # Sequences with a custom base path cannot be handled by the standard path
    # construction below; submit them manually with the explicit WORKDIR.
    if [[ -n "$custom_path" ]]; then
        echo "SKIP (custom path — submit manually): $seq_id"
        (( skipped++ )) || true
        continue
    fi

    subdir="${LABEL_MAP[$seq_type]:-}"
    if [[ -z "$subdir" ]]; then
        echo "WARNING: unknown seq_type '$seq_type' for $seq_id — skipping"
        (( skipped++ )) || true
        continue
    fi

    WORKDIR="${BASE}/${subdir}/${seq_id}/${RUNREL}"

    if [[ ! -d "$WORKDIR" ]]; then
        echo "WARNING: WORKDIR not found (job will fail at runtime): $WORKDIR"
        (( missing++ )) || true
    fi

    echo "Submitting: $seq_id  [$seq_type → $subdir]"
    echo "           WORKDIR: $WORKDIR"
    sbatch --job-name="pipeline_${seq_id}" "$WORKER" "$WORKDIR" "$CONFIG" 500000 "$PHASE" "$FORCE"
    (( submitted++ )) || true

done < "$SEQ_LIST"

echo ""
echo "=== Done ==="
echo "  Submitted : $submitted jobs"
echo "  Skipped   : $skipped sequences"
if (( missing > 0 )); then
    echo "  WARNING   : $missing submitted jobs have a missing WORKDIR — check error logs"
fi
