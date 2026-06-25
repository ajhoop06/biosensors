"""
aggregate_contact_features.py
------------------------------
After all SLURM jobs finish, run this on the login node to collect
per-sequence *_contact_summary.csv files into one table ready to
merge with feat_table.xlsx.

Usage:
    python aggregate_contact_features.py
"""

import os
import glob
import pandas as pd

# ─────────────────────────────────────────────
# PATHS — mirror contact_type_analysis.py
# ─────────────────────────────────────────────
base         = "/scratch/alpine/ivta1597/LCA_boltz_models/LIG_contacts"
results_dir  = os.path.join(base, "contact_type_results")
seq_ids_file = os.path.join(base, "seq_ids.txt")
out_path     = os.path.join(results_dir, "contact_features_all.csv")

# ─────────────────────────────────────────────
# LOAD
# ─────────────────────────────────────────────
summary_files = sorted(glob.glob(os.path.join(results_dir, "*_contact_summary.csv")))
print(f"Found {len(summary_files)} summary files")

if not summary_files:
    raise FileNotFoundError(f"No summary CSVs found in: {results_dir}")

dfs      = [pd.read_csv(f) for f in summary_files]
combined = pd.concat(dfs, ignore_index=True)

# ─────────────────────────────────────────────
# CHECK COVERAGE vs seq_ids.txt
# ─────────────────────────────────────────────
if os.path.exists(seq_ids_file):
    with open(seq_ids_file) as fh:
        all_ids = [l.split()[0] for l in fh if l.strip()]   # first column only
    missing = set(all_ids) - set(combined["seq_id"])
    if missing:
        print(f"WARNING: {len(missing)} sequences missing results: {missing}")
    else:
        print("All sequences accounted for.")

print(f"\nFeature table shape: {combined.shape}")
print(combined.to_string())

combined.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")
