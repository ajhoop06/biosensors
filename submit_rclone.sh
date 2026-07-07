#!/bin/bash

BASE_SCRATCH=/scratch/alpine/ivta1597/LCA_boltz_models
BASE_ONEDRIVE="onedrive_ivana_cu:Shirts Lab/LCA_boltz_models"

declare -A SEQUENCES
SEQUENCES[pair_0272_nb]=nonbinders
#SEQUENCES[nonb_019_nb]=nonbinders

skipped=()

for seq_id in "${!SEQUENCES[@]}"; do
    type_dir=${SEQUENCES[$seq_id]}
    src="${BASE_SCRATCH}/${type_dir}/${seq_id}/"

    if [[ ! -d "$src" ]]; then
        skipped+=("$seq_id")
        continue
    fi

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=rclone_${seq_id}
#SBATCH --output=/projects/ivta1597/biosensors/rclone_logs/output_${seq_id}_%j.out
#SBATCH --error=/projects/ivta1597/biosensors/rclone_logs/error_${seq_id}_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=FAIL

module load slurm/alpine
module load rclone/1.58.0

rclone copy \
    "${src}" \
    "${BASE_ONEDRIVE}/${type_dir}/${seq_id}/" \
    --transfers 4 \
    --checkers 8 \
    --log-file /projects/ivta1597/rclone_logs/rclone_${seq_id}_\${SLURM_JOB_ID}.log \
    --log-level INFO
EOF

    echo "Submitted job for $seq_id"
done

if [[ ${#skipped[@]} -gt 0 ]]; then
    echo "Skipped (directory not found):"
    for s in "${skipped[@]}"; do
        echo "  $s"
    done
fi
