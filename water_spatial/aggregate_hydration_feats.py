"""
aggregate_hydration_feats.py

Reads each sequence's {seq_id}_hydration_{TAG}.csv (produced by
hydration_calc.py, one row per frame: time_ns, hydration_count) and computes
per-sequence summary features, mirroring pkt_vol/extract_pkt_feats.py's
early/late/drift/slope column family for consistency with the rest of the
feature table.

Early/late windows are fractional (first/last 20% of frames, matching
extract_pkt_feats.py) rather than a fixed-ns window, since sequences can vary
in how many frames they actually cover after striding. Because hydration_calc.py
records real simulation time (unlike mdpocket's descriptors.txt), both a
fractional-position slope and an absolute-ns slope are emitted.

Output: water_density_feats.csv

Usage:
    python aggregate_hydration_feats.py [seq_ids_ngs_observed.txt] [--start-ns 40] [--end-ns 500]
                                         [--early-late-frac 0.2] [--out water_density_feats.csv]
                                         [--structure-source {ngs_observed,designed_assumed,all}]
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.stats import linregress

# ── Configurable paths ────────────────────────────────────────────────────────
BASE = "/scratch/alpine/ivta1597/LCA_boltz_models"
# ─────────────────────────────────────────────────────────────────────────────

TYPE_SUBDIR = {
    "Binder":         "binders",
    "False Positive": "nonbinders",
    "Low Confidence": "neg_low_pkt",
    "Fail Geometry":  "neg_fail_gate",
}

# md_candidate_guide.csv's md_group values use different suffixes than the
# seq_id naming convention used everywhere else in the repo.
MD_GROUP_SUFFIX = {
    "binder": "binder",
    "non_binder": "nb",
    "negative_low_pocket": "low_pkt",
    "negative_fail_gate": "fail_gate",
}


def load_source_ids(guide_path, source):
    """seq_ids from md_candidate_guide.csv where source == `source` (e.g.
    ngs_observed = sequencing-confirmed via Y2H/FACS sort-seq, vs.
    designed_assumed)."""
    guide = pd.read_csv(guide_path)
    matched = guide[guide["source"] == source].copy()
    matched["seq_id"] = matched.apply(
        lambda r: f"{r['pair_id']}_{MD_GROUP_SUFFIX.get(r['md_group'], r['md_group'])}",
        axis=1)
    return set(matched["seq_id"])


def extract_features(df, window_frac):
    """Compute summary statistics from a per-frame hydration_count series,
    including early-vs-late trajectory comparison (fractional window, since
    frame counts can vary across sequences after striding)."""
    t   = df["time_ns"].values
    cnt = df["hydration_count"].values
    n   = len(cnt)
    pct = int(round(window_frac * 100))

    k = max(1, int(round(window_frac * n)))
    early_cnt = cnt[:k]
    late_cnt  = cnt[n - k:]
    early_mean = float(early_cnt.mean())
    late_mean  = float(late_cnt.mean())

    pos = np.arange(n) / (n - 1) if n > 1 else np.zeros(n)
    frac_slope = float(linregress(pos, cnt).slope) if n > 1 else 0.0
    ns_slope   = float(linregress(t, cnt).slope) if n > 1 else 0.0

    return {
        "hydration_count_mean":        float(np.mean(cnt)),
        "hydration_count_std":         float(np.std(cnt)),
        "hydration_count_min":         float(np.min(cnt)),
        "hydration_count_max":         float(np.max(cnt)),
        f"hydration_count_early{pct}_mean": early_mean,
        f"hydration_count_late{pct}_mean":  late_mean,
        f"hydration_count_drift{pct}":      late_mean - early_mean,
        "hydration_count_slope":       frac_slope,
        "hydration_count_slope_per_ns": ns_slope,
        "n_frames":                    n,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("seq_list", nargs="?",
                        default="/projects/ivta1597/biosensors/seq_ids_ngs_observed.txt")
    parser.add_argument("--start-ns", type=float, default=40.0)
    parser.add_argument("--end-ns",   type=float, default=500.0)
    parser.add_argument("--early-late-frac", type=float, default=0.2,
                        help="Fraction of each sequence's own frames counted as "
                             "'early'/'late' for the drift comparison (default: 0.2).")
    parser.add_argument("--output", default="water_density_feats.csv",
                        help="Output CSV filename (default: water_density_feats.csv)")
    parser.add_argument("--base", default=BASE)
    parser.add_argument("--structure-source", default="all",
                        choices=["ngs_observed", "designed_assumed", "all"],
                        help="Filter sequences by md_candidate_guide.csv's source column. "
                             "'all' (default) applies no filtering.")
    parser.add_argument("--structure-guide",
                        default="/projects/ivta1597/biosensors/md_candidate_guide.csv",
                        help="Path to md_candidate_guide.csv (default: %(default)s)")
    args = parser.parse_args()

    if not os.path.exists(args.seq_list):
        print(f"ERROR: seq list not found: {args.seq_list}")
        sys.exit(1)

    source_ids = None
    if args.structure_source != "all":
        source_ids = load_source_ids(args.structure_guide, args.structure_source)
        print(f"Restricting to source == {args.structure_source}: "
              f"{len(source_ids)} in {args.structure_guide}")

    TAG = f"{int(args.start_ns)}_{int(args.end_ns)}ns"

    records = []
    missing = []
    skipped_wrong_source = []

    with open(args.seq_list) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            parts       = line.split("\t")
            seq_id      = parts[0].strip()
            seq_type    = parts[1].strip() if len(parts) > 1 else ""
            custom_path = parts[2].strip() if len(parts) > 2 else ""

            if source_ids is not None and seq_id not in source_ids:
                skipped_wrong_source.append(seq_id)
                continue

            if custom_path:
                run_dir = os.path.join(custom_path, f"water_density_{TAG}")
            else:
                dir_type = TYPE_SUBDIR.get(seq_type, seq_type)
                run_dir  = os.path.join(args.base, dir_type, seq_id, f"water_density_{TAG}")

            csv_file = os.path.join(run_dir, f"{seq_id}_hydration_{TAG}.csv")

            if not os.path.exists(csv_file):
                print(f"MISSING: {seq_id}  [{seq_type}]  ->  {csv_file}")
                missing.append(seq_id)
                continue

            try:
                df       = pd.read_csv(csv_file)
                features = extract_features(df, args.early_late_frac)
                features["seq_id"]   = seq_id
                features["seq_type"] = seq_type
                records.append(features)
                print(f"OK: {seq_id}  [{seq_type}]  "
                      f"mean={features['hydration_count_mean']:.2f}  "
                      f"n={features['n_frames']}")
            except Exception as e:
                print(f"ERROR: {seq_id}  -  {e}")
                missing.append(seq_id)

    if not records:
        print("\nNo sequences loaded - nothing to write.")
        sys.exit(1)

    feat_df = pd.DataFrame(records)
    pct = int(round(args.early_late_frac * 100))
    col_order = ["seq_id", "seq_type", "hydration_count_mean", "hydration_count_std",
                 "hydration_count_min", "hydration_count_max",
                 f"hydration_count_early{pct}_mean", f"hydration_count_late{pct}_mean",
                 f"hydration_count_drift{pct}",
                 "hydration_count_slope", "hydration_count_slope_per_ns", "n_frames"]
    feat_df = feat_df[col_order]
    feat_df.to_csv(args.output, index=False)

    print(f"\nFeatures written to: {args.output}")
    print(f"  Sequences processed : {len(records)}")
    print(f"  Sequences missing   : {len(missing)}")
    if missing:
        print(f"  Missing seq_ids     : {', '.join(missing)}")
    if source_ids is not None:
        print(f"  Skipped (source != {args.structure_source}) : {len(skipped_wrong_source)}")


if __name__ == "__main__":
    main()
