#!/usr/bin/env bash
# Check for expected trajectory files (EM: .trr, NVT/NPT/prod_md: .xtc)
# for every sequence listed in seq_ids_orig.txt.
#
# Usage: bash check_xtc_files.sh [seq_ids_file]
#   Default seq_ids file: seq_ids_orig.txt (relative to this script)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SEQ_IDS="${1:-$SCRIPT_DIR/seq_ids_orig.txt}"
BASE="/Users/ivanatang/OneDrive - UCB-O365/Shirts Lab/LCA_boltz_models"
PROD_DIR="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"

ok_count=0
issue_count=0
missing_dir_count=0

while IFS=$'\t ' read -r seq_id _label rest; do
    [[ -z "$seq_id" ]] && continue

    # Resolve subdirectory and full directory name from seq_id suffix
    if   [[ "$seq_id" == *"_binder"   ]]; then subdir="binders";      dir_name="$seq_id"
    elif [[ "$seq_id" == *"_nb"       ]]; then subdir="nonbinders";   dir_name="$seq_id"
    elif [[ "$seq_id" == *"_low_pkt"  ]]; then subdir="neg_low_pkt";  dir_name="$seq_id"
    elif [[ "$seq_id" == *"_fail_gate" ]]; then subdir="neg_fail_gate"; dir_name="$seq_id"
    elif [[ "$seq_id" == bind_*       ]]; then subdir="binders";      dir_name="${seq_id}_binder"
    elif [[ "$seq_id" == nonb_*       ]]; then subdir="nonbinders";   dir_name="${seq_id}_nb"
    else
        echo "UNKNOWN  $seq_id"
        continue
    fi

    seq_path="$BASE/$subdir/$dir_name"

    if [[ ! -d "$seq_path" ]]; then
        printf "MISSING_DIR  %-30s  %s\n" "$seq_id" "$seq_path"
        (( missing_dir_count++ ))
        continue
    fi

    issues=""

    # EM — GROMACS energy minimization writes .trr, not .xtc
    if [[ ! -d "$seq_path/EM" ]]; then
        issues+=" EM_dir_missing"
    elif [[ -z "$(find "$seq_path/EM" -maxdepth 1 -name "*.trr" 2>/dev/null)" ]]; then
        issues+=" EM_no_trr"
    fi

    # NVT
    if [[ ! -d "$seq_path/NVT" ]]; then
        issues+=" NVT_dir_missing"
    elif [[ -z "$(find "$seq_path/NVT" -maxdepth 1 -name "*.xtc" 2>/dev/null)" ]]; then
        issues+=" NVT_no_xtc"
    fi

    # NPT
    if [[ ! -d "$seq_path/NPT" ]]; then
        issues+=" NPT_dir_missing"
    elif [[ -z "$(find "$seq_path/NPT" -maxdepth 1 -name "*.xtc" 2>/dev/null)" ]]; then
        issues+=" NPT_no_xtc"
    fi

    # Production MD
    if [[ ! -d "$seq_path/$PROD_DIR" ]]; then
        issues+=" prod_dir_missing"
    else
        if [[ -z "$(find "$seq_path/$PROD_DIR" -maxdepth 1 -name "*.xtc" 2>/dev/null)" ]]; then
            issues+=" prod_no_xtc"
        fi

        # Check simulation reached 500 ns by reading the last Step/Time entry in the log
        prod_log="$seq_path/$PROD_DIR/prod_md_500ns.log"
        if [[ ! -f "$prod_log" ]]; then
            issues+=" prod_no_log"
        else
            last_time=$(grep -a -A1 "^ *Step *Time$" "$prod_log" \
                        | grep -v "Step\|--" | awk '{print $2}' | tail -1)
            if [[ -z "$last_time" ]]; then
                issues+=" prod_log_unreadable"
            elif ! awk -v t="$last_time" 'BEGIN{exit !(t >= 500000)}'; then
                printf_time=$(awk -v t="$last_time" 'BEGIN{printf "%.1f", t/1000}')
                issues+=" prod_incomplete(${printf_time}ns)"
            fi
        fi
    fi

    if [[ -z "$issues" ]]; then
        printf "OK       %s\n" "$seq_id"
        (( ok_count++ ))
    else
        printf "ISSUES   %-30s  %s\n" "$seq_id" "$issues"
        (( issue_count++ ))
    fi

done < "$SEQ_IDS"

echo ""
echo "Summary: $ok_count OK | $issue_count with issues | $missing_dir_count missing dirs"
