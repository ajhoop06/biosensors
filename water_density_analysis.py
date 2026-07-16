"""
water_density_analysis.py — Does water density in/around the LCA pocket
distinguish binders from nonbinders, and is it worth adding to the ML
feature table?

Two independently-computed water-density metrics are evaluated:
  hydration_count_mean        water_analysis/hydration_calc.py + aggregate_hydration_feats.py
                               (mdtraj: water oxygens within 3.5 A of ligand heavy atoms)
  pocket_water_density_mean   water_spatial/{water_spatial_prep,run}.sh + extract_water_spatial_feats.py
                               (gmx spatial: voxel density within 8 A of the ligand centroid)

This script does two things, in order:

  1. CORRELATION CHECK (always reported, regardless of outcome): are the two
     water metrics redundant with pocket_vol_mean (pkt_vol/pocket_volume_features.csv)
     or with the existing R-score water-mediated-contact signal
     (mean W across residues, water_analysis/dw_scores_*_ml.csv)? High
     correlation doesn't invalidate the science, but it does mean a metric
     adds little independent signal for automated feature selection.

  2. SIGNIFICANCE GATE (Mann-Whitney U, Binder vs each other group, with
     Benjamini-Hochberg FDR correction across all tested columns): per
     project decision, a water-density column is only a candidate for
     ML_classification.ipynb integration if it clears BH-adjusted q < 0.05
     against at least one other group. A feature that fails this bar is
     still a valid, reportable scientific result -- it's just not added to
     the model, since a non-discriminating feature only adds noise/
     overfitting risk to GroupAwareSelector rather than signal.

Usage:
    python water_density_analysis.py \
        --hydration-csv water_analysis/water_density_feats.csv \
        --spatial-csv   water_spatial/water_spatial_feats.csv \
        --pocket-vol-csv pkt_vol/pocket_volume_features.csv \
        --dw-scores-csv water_analysis/dw_scores_all_sequences_40_500ns_ml.csv \
        --seq-list      seq_ids_ngs_observed.txt \
        --out-dir       analysis/water_density
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr, pearsonr
from sklearn.metrics import roc_auc_score

GROUP_COLOR = {
    "Binder":         "#648FFF",
    "False Positive": "#DC267F",
    "Low Confidence": "#FE6100",
    "Fail Geometry":  "#FFB000",
}
GROUP_ORDER = ["Binder", "False Positive", "Low Confidence", "Fail Geometry"]
GROUP_A = "Binder"

SIG_THRESHOLD = 0.05

# Column families to gate: the two headline means plus their early/late/
# drift/slope companions (from aggregate_hydration_feats.py /
# extract_water_spatial_feats.py).
HYDRATION_COLS_PREFIX = "hydration_count_"
SPATIAL_COLS_PREFIX   = "pocket_water_density_"


# ── Stats helpers (mirrors core_vs_tail_regions.py) ─────────────────────────
def cohens_d(a, b):
    na, nb = len(a), len(b)
    pooled_var = ((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2)
    if pooled_var <= 0:
        return 0.0
    return (a.mean() - b.mean()) / np.sqrt(pooled_var)


def rank_auc(a, b):
    y = np.r_[np.ones(len(a)), np.zeros(len(b))]
    scores = np.r_[a, b]
    return roc_auc_score(y, scores)


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out


def compare_groups(df, col, group_a, group_b, min_n=5):
    a = df.loc[df["Group"] == group_a, col].dropna().values
    b = df.loc[df["Group"] == group_b, col].dropna().values
    if len(a) < min_n or len(b) < min_n:
        return None
    _, p = mannwhitneyu(a, b, alternative="two-sided")
    return dict(group_a=group_a, group_b=group_b, n_a=len(a), n_b=len(b),
                mean_a=a.mean(), mean_b=b.mean(),
                cohens_d=cohens_d(a, b), auc=rank_auc(a, b), p=p)


def box_jitter(ax, df, col, groups, rng=None):
    rng = rng or np.random.default_rng(42)
    xpos = 0
    for group in groups:
        vals = df.loc[df["Group"] == group, col].dropna().values
        if len(vals) == 0:
            xpos += 1
            continue
        color = GROUP_COLOR[group]
        ax.boxplot(vals, positions=[xpos], widths=0.6, patch_artist=True,
                   medianprops=dict(color="black", linewidth=1.5),
                   boxprops=dict(facecolor=color, alpha=0.5),
                   whiskerprops=dict(color=color), capprops=dict(color=color),
                   flierprops=dict(marker="", linestyle="none"))
        jitter = rng.uniform(-0.15, 0.15, len(vals))
        ax.scatter(xpos + jitter, vals, color=color, s=14, alpha=0.8, zorder=3)
        xpos += 1


# ── Loading ───────────────────────────────────────────────────────────────
def load_seq_type_map(seq_list_path):
    mapping = {}
    with open(seq_list_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            mapping[parts[0]] = parts[1] if len(parts) > 1 else "Unknown"
    return mapping


def load_optional(path, label):
    if path and os.path.exists(path):
        df = pd.read_csv(path)
        print(f"Loaded {label}: {path}  ({len(df)} rows)")
        return df
    print(f"  MISSING (skipping): {label}  ->  {path}")
    return None


def mean_w_score(dw_df):
    """Mean water-mediated-contact (W) score across all residue columns, one
    value per sequence -- the existing R-score signal to correlate against."""
    w_cols = [c for c in dw_df.columns if c.startswith("W_")]
    out = dw_df[["seq_id"]].copy()
    out["mean_W_score"] = dw_df[w_cols].mean(axis=1)
    return out


def gate_columns(prefix, df):
    return [c for c in df.columns if c.startswith(prefix)]


# ── Main ─────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hydration-csv", default="water_analysis/water_density_feats.csv")
    parser.add_argument("--spatial-csv", default="water_spatial/water_spatial_feats.csv")
    parser.add_argument("--pocket-vol-csv", default="pkt_vol/pocket_volume_features.csv")
    parser.add_argument("--dw-scores-csv",
                        default="water_analysis/dw_scores_all_sequences_40_500ns_ml.csv")
    parser.add_argument("--seq-list", default="seq_ids_ngs_observed.txt")
    parser.add_argument("--out-dir", default="analysis/water_density")
    parser.add_argument("--sig-threshold", type=float, default=SIG_THRESHOLD)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    seq_type_map = load_seq_type_map(args.seq_list)

    hydration_df = load_optional(args.hydration_csv, "hydration-shell counts")
    spatial_df   = load_optional(args.spatial_csv, "gmx-spatial pocket density")
    pocket_vol_df = load_optional(args.pocket_vol_csv, "pocket volume")
    dw_df        = load_optional(args.dw_scores_csv, "R-score D/W table")

    if hydration_df is None and spatial_df is None:
        raise FileNotFoundError(
            "Neither hydration-shell nor gmx-spatial water-density CSVs were "
            "found. Run the water_analysis/ or water_spatial/ pipelines first.")

    # ── Assemble one merged table, keyed on seq_id, with a Group column ─────
    base = pd.DataFrame({"seq_id": list(seq_type_map.keys())})
    base["Group"] = base["seq_id"].map(seq_type_map)
    merged = base

    for df, label in [(hydration_df, "hydration"), (spatial_df, "spatial")]:
        if df is not None:
            merged = merged.merge(df.drop(columns=["seq_type"], errors="ignore"),
                                   on="seq_id", how="left")

    if pocket_vol_df is not None:
        merged = merged.merge(
            pocket_vol_df[["seq_id", "pocket_vol_mean"]], on="seq_id", how="left")

    if dw_df is not None:
        w_df = mean_w_score(dw_df)
        merged = merged.merge(w_df, on="seq_id", how="left")

    merged_path = os.path.join(args.out_dir, "water_density_merged.csv")
    merged.to_csv(merged_path, index=False)
    print(f"\nSaved merged table -> {merged_path}  ({len(merged)} sequences)")

    # ── 1. Correlation check (always reported) ──────────────────────────────
    print("\n" + "=" * 70)
    print("CORRELATION CHECK — redundancy with existing pocket/R-score features")
    print("=" * 70)

    corr_pairs = []
    candidate_cols = {
        "hydration_count_mean":       "hydration_count_mean" in merged.columns,
        "pocket_water_density_mean":  "pocket_water_density_mean" in merged.columns,
    }
    reference_cols = {
        "pocket_vol_mean": "pocket_vol_mean" in merged.columns,
        "mean_W_score":    "mean_W_score" in merged.columns,
    }

    corr_targets = []
    if candidate_cols["hydration_count_mean"] and candidate_cols["pocket_water_density_mean"]:
        corr_targets.append(("hydration_count_mean", "pocket_water_density_mean"))
    for cand, present in candidate_cols.items():
        if not present:
            continue
        for ref, ref_present in reference_cols.items():
            if ref_present:
                corr_targets.append((cand, ref))

    for col_a, col_b in corr_targets:
        sub = merged[[col_a, col_b]].dropna()
        if len(sub) < 5:
            continue
        rho, p_s = spearmanr(sub[col_a], sub[col_b])
        r, p_p = pearsonr(sub[col_a], sub[col_b])
        corr_pairs.append(dict(feature_a=col_a, feature_b=col_b, n=len(sub),
                               spearman_rho=rho, spearman_p=p_s,
                               pearson_r=r, pearson_p=p_p))
        flag = "REDUNDANT (|r|>0.7)" if abs(r) > 0.7 else "complementary"
        print(f"  {col_a} vs {col_b}: spearman rho={rho:.3f}, pearson r={r:.3f}  [{flag}]")

    if corr_pairs:
        corr_df = pd.DataFrame(corr_pairs)
        corr_path = os.path.join(args.out_dir, "water_density_correlations.csv")
        corr_df.to_csv(corr_path, index=False)
        print(f"\nSaved {corr_path}")
    else:
        print("  Not enough overlapping data to compute correlations.")

    # ── 2. Significance gate (Binder vs each other group, BH-FDR) ───────────
    print("\n" + "=" * 70)
    print(f"SIGNIFICANCE GATE — {GROUP_A} vs each other group "
          f"(Mann-Whitney U, BH-FDR q < {args.sig_threshold})")
    print("=" * 70)

    gate_cols = gate_columns(HYDRATION_COLS_PREFIX, merged) + \
                gate_columns(SPATIAL_COLS_PREFIX, merged)
    gate_cols = [c for c in gate_cols if c not in
                 ("pocket_water_density_missing", "pocket_water_density_n_voxels")]

    other_groups = [g for g in GROUP_ORDER if g != GROUP_A]
    rows = []
    for col in gate_cols:
        for group_b in other_groups:
            res = compare_groups(merged, col, GROUP_A, group_b)
            if res is None:
                continue
            rows.append(dict(feature=col, **res))

    passed_features = set()
    if rows:
        stats_df = pd.DataFrame(rows)
        stats_df["q"] = bh_fdr(stats_df["p"].values)
        stats_path = os.path.join(args.out_dir, "water_density_significance_gate.csv")
        stats_df.sort_values("p").to_csv(stats_path, index=False)
        print(stats_df.sort_values("p").to_string(index=False))
        print(f"\nSaved {stats_path}")

        sig = stats_df[stats_df["q"] < args.sig_threshold]
        passed_features = set(sig["feature"].unique())
        if passed_features:
            print(f"\nPASSED gate ({len(passed_features)} column(s), q < {args.sig_threshold}):")
            for feat in sorted(passed_features):
                best = sig[sig["feature"] == feat].sort_values("q").iloc[0]
                print(f"  {feat}: {GROUP_A} vs {best.group_b}  "
                      f"q={best.q:.4g}  d={best.cohens_d:.2f}")
        else:
            print(f"\nNO columns cleared q < {args.sig_threshold}. "
                  "Per project decision, do NOT add these features to "
                  "ML_classification.ipynb -- report the trend descriptively instead.")
    else:
        print("  Not enough data to run the significance gate.")

    gate_summary_path = os.path.join(args.out_dir, "ml_integration_gate.txt")
    with open(gate_summary_path, "w") as f:
        if passed_features:
            f.write("Columns cleared for ML_classification.ipynb integration "
                    f"(q < {args.sig_threshold} vs at least one other group):\n")
            for feat in sorted(passed_features):
                f.write(f"  {feat}\n")
        else:
            f.write(f"No water-density columns cleared q < {args.sig_threshold} against "
                    f"{GROUP_A}. Do not add to ML_classification.ipynb.\n")
    print(f"Saved {gate_summary_path}")

    # ── Plots: headline metrics by group ─────────────────────────────────────
    headline_cols = [c for c in ["hydration_count_mean", "pocket_water_density_mean"]
                      if c in merged.columns]
    if headline_cols:
        fig, axes = plt.subplots(1, len(headline_cols),
                                  figsize=(6 * len(headline_cols), 4.5),
                                  dpi=300, constrained_layout=True, squeeze=False)
        for i, col in enumerate(headline_cols):
            ax = axes[0][i]
            box_jitter(ax, merged, col, GROUP_ORDER)
            ax.set_xticks(range(len(GROUP_ORDER)))
            ax.set_xticklabels(GROUP_ORDER, rotation=20, ha="right", fontsize=8)
            ax.set_ylabel(col)
            sig_marker = " *" if col in passed_features else ""
            ax.set_title(f"{col}{sig_marker}", fontsize=10)
            ax.grid(True, alpha=0.4)
        fig.suptitle(f"Water density by group ({GROUP_A} vs others)\n"
                     "* = cleared the BH-FDR significance gate", fontsize=12)
        plot_path = os.path.join(args.out_dir, "water_density_by_group.png")
        fig.savefig(plot_path)
        plt.close(fig)
        print(f"\nSaved {plot_path}")

    print("\nDone. All outputs written to:", args.out_dir)


if __name__ == "__main__":
    main()
