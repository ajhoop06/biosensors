"""
core_vs_tail_regions.py — Is there a distinguishable difference
between Binder and nonbinder sequences when looking at only the LCA
ligand's core (steroid ring system) or tail (C20-C24 carboxylate chain)?

Primary comparison: Binder vs False Positive. Both went through the same
computational pipeline; the only difference is the wet-lab outcome, which
makes this the most direct "binder vs nonbinder" contrast available. Low
Confidence and Fail Geometry sequences are shown in the plots for context
(their pocket/geometry didn't pass QC in the first place, which is a
different kind of exclusion, and there are only ~9-10 of each) but are not
part of the headline significance test.

Core and tail are analyzed COMPLETELY SEPARATELY — every plot and stats
table answers "is there a Binder vs False Positive difference in THIS
region alone?" with no cross-region encoding (no side-by-side panels, no
shared color axis). The two regions are never compared to each other here;
if you want that, see compare_ligand_regions.py instead. The only place
core and tail appear together is the final printed/saved verdict table,
which is plain text/CSV (region, n significant, top hit), not a chart.

For both the contact-type features and the per-residue R-scores, this
script, per region:
  1. Tests Binder vs False Positive (Mann-Whitney U) for every feature /
     residue, with a Cohen's d effect size and a rank-AUC (0.5 = no
     separation, 1.0 = perfect separation, <0.5 = False Positive has the
     higher value) alongside the p-value.
  2. Applies a Benjamini-Hochberg FDR correction (this matters — many
     residues get tested, so a handful of uncorrected "hits" is expected
     by chance alone; the q-value column tells you which survive that).
  3. Plots: one figure with every contact feature (Binder vs the other
     3 groups, headline test vs False Positive), one significance plot
     across all residues, and one figure with the top-hit residues —
     all scoped to that region only.

Inputs (produced by aggregate_r_scores.py / agg_contact_feats.py with
--ligand-region core|tail):
    r_scores_all_sequences_{TAG}_core.csv / _tail.csv
    contact_features_all_{TAG}_core.csv / _tail.csv

Usage:
    python core_vs_tail_regions.py \
        --r-scores-dir /projects/ivta1597/biosensors/water_analysis/agg_out \
        --contact-dir  /projects/ivta1597/biosensors/LIG_contacts \
        --seq-list     /projects/ivta1597/biosensors/seq_ids_orig.txt \
        --tag 40_500ns \
        --out-dir      /projects/ivta1597/biosensors/analysis/core_vs_tail
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

# ── Constants ─────────────────────────────────────────────────────────────────
GROUP_COLOR = {
    "Binder":         "#648FFF",
    "False Positive": "#DC267F",
    "Low Confidence": "#FE6100",
    "Fail Geometry":  "#FFB000",
}
GROUP_ORDER = ["Binder", "False Positive", "Low Confidence", "Fail Geometry"]

REGIONS = ["core", "tail"]

GROUP_A, GROUP_B = "Binder", "False Positive"   # the headline test

CONTACT_FEATURES = ["mean_frac_hydrophobic", "mean_n_total", "mean_n_hydrophobic",
                     "mean_n_polar", "mean_n_pos_charged", "mean_n_neg_charged"]


# ── Stats helpers ─────────────────────────────────────────────────────────────
def cohens_d(a, b):
    na, nb = len(a), len(b)
    pooled_var = ((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2)
    if pooled_var <= 0:
        return 0.0
    return (a.mean() - b.mean()) / np.sqrt(pooled_var)


def rank_auc(a, b):
    """AUC of using the feature value to separate group a (label 1) from
    group b (label 0). 0.5 = no separation; >0.5 = a tends higher."""
    y = np.r_[np.ones(len(a)), np.zeros(len(b))]
    scores = np.r_[a, b]
    return roc_auc_score(y, scores)


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR-adjusted q-values."""
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


def compare_groups(df, col, group_a=GROUP_A, group_b=GROUP_B, min_n=5):
    a = df.loc[df["seq_type"] == group_a, col].dropna().values
    b = df.loc[df["seq_type"] == group_b, col].dropna().values
    if len(a) < min_n or len(b) < min_n:
        return None
    _, p = mannwhitneyu(a, b, alternative="two-sided")
    return dict(n_binder=len(a), n_fp=len(b), mean_binder=a.mean(), mean_fp=b.mean(),
                cohens_d=cohens_d(a, b), auc=rank_auc(a, b), p=p)


# ── Loading ───────────────────────────────────────────────────────────────────
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


def load_r_scores(r_scores_dir, tag):
    out = {}
    for region in REGIONS:
        path = os.path.join(r_scores_dir, f"r_scores_all_sequences_{tag}_{region}.csv")
        if os.path.exists(path):
            out[region] = pd.read_csv(path)
            print(f"Loaded {region} R-scores: {path}  ({len(out[region])} sequences)")
        else:
            print(f"  MISSING: {path}")
    return out


def load_contact_features(contact_dir, tag, seq_type_map):
    out = {}
    for region in REGIONS:
        path = os.path.join(contact_dir, f"contact_features_all_{tag}_{region}.csv")
        if os.path.exists(path):
            df = pd.read_csv(path)
            df["seq_type"] = df["seq_id"].map(seq_type_map)
            out[region] = df
            print(f"Loaded {region} contact features: {path}  ({len(df)} sequences)")
        else:
            print(f"  MISSING: {path}")
    return out


# md_candidate_guide.csv's md_group values use different suffixes than the
# seq_id naming convention (pair_XXXX_binder / _nb / _low_pkt / _fail_gate)
# used everywhere else in the repo.
MD_GROUP_SUFFIX = {
    "binder": "binder",
    "non_binder": "nb",
    "negative_low_pocket": "low_pkt",
    "negative_fail_gate": "fail_gate",
}


def load_ngs_observed_ids(structure_source_path):
    """seq_ids from md_candidate_guide.csv where source == ngs_observed, i.e.
    sequencing-confirmed (Y2H/FACS sort-seq) rather than designed-and-assumed."""
    guide = pd.read_csv(structure_source_path)
    confirmed = guide[guide["source"] == "ngs_observed"].copy()
    confirmed["seq_id"] = confirmed.apply(
        lambda r: f"{r['pair_id']}_{MD_GROUP_SUFFIX.get(r['md_group'], r['md_group'])}",
        axis=1)
    return set(confirmed["seq_id"])


def filter_to_ngs_observed(data, ngs_ids, label):
    for region, df in data.items():
        before = len(df)
        data[region] = df[df["seq_id"].isin(ngs_ids)].reset_index(drop=True)
        print(f"  {label} {region}: {before} -> {len(data[region])} sequences "
              f"(sequencing-confirmed only)")
    return data


def resid_columns(df):
    return sorted((c for c in df.columns if c.startswith("R_")),
                  key=lambda c: int(c.split("_")[1]))


# ── Box + jitter helper (matches the repo's plot_Rg_sasa.py style) ─────────────
def box_jitter(ax, df, col, groups, xpos_start=0, rng=None):
    rng = rng or np.random.default_rng(42)
    xpos = xpos_start
    for group in groups:
        vals = df.loc[df["seq_type"] == group, col].dropna().values
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
    return xpos


# ── Contact-type features: Binder vs False Positive, ONE region at a time ─────
def contact_feature_region_analysis(df, region, out_dir):
    """Self-contained answer to: within `region` ALONE, does Binder differ
    from False Positive in contact-type composition? Nothing in this
    function's output references the other region — that's the point."""
    rows = []
    for feature in CONTACT_FEATURES:
        if feature not in df.columns:
            continue
        res = compare_groups(df, feature)
        if res is None:
            continue
        rows.append(dict(feature=feature, **res))
    if not rows:
        print(f"\n[Contact features - {region}] Not enough Binder/False Positive data — skipping.")
        return None
    stats_df = pd.DataFrame(rows)
    stats_df["q"] = bh_fdr(stats_df["p"].values)
    stats_path = os.path.join(out_dir, f"contact_features_{region}_binder_vs_fp_stats.csv")
    stats_df.sort_values("p").to_csv(stats_path, index=False)
    print(f"\nSaved {stats_path}")
    print(stats_df.sort_values("p").to_string(index=False))

    n_sig = int((stats_df["q"] < 0.05).sum())
    print(f"[Contact features - {region}] {n_sig}/{len(stats_df)} features "
          f"FDR-significant ({GROUP_A} vs {GROUP_B}, q<0.05)")

    # ── One figure, this region only: every feature, Binder vs False Positive ──
    features_present = stats_df["feature"].tolist()
    n = len(features_present)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows),
                              dpi=300, constrained_layout=True, squeeze=False)
    for i, feature in enumerate(features_present):
        ax = axes[i // ncols][i % ncols]
        box_jitter(ax, df, feature, GROUP_ORDER)
        ax.set_xticks(range(len(GROUP_ORDER)))
        ax.set_xticklabels(GROUP_ORDER, rotation=20, ha="right", fontsize=8)
        row = stats_df[stats_df.feature == feature].iloc[0]
        star = "  *SIGNIFICANT*" if row["q"] < 0.05 else ""
        ax.set_title(f"{feature.replace('mean_', '')}\n"
                     f"p={row.p:.2g}, q={row.q:.2g}, d={row.cohens_d:.2f}{star}", fontsize=9)
        ax.grid(True, alpha=0.4)
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle(f"{region.upper()} REGION ONLY — contact-type composition by group\n"
                 f"{GROUP_A} vs {GROUP_B}: {n_sig}/{len(stats_df)} features significant (FDR q<0.05)",
                 fontsize=13)
    path = os.path.join(out_dir, f"contact_features_{region}_binder_vs_nonbinder.png")
    fig.savefig(path)
    plt.close(fig)
    print(f"Saved {path}")

    top = stats_df.sort_values("p").iloc[0]
    return dict(region=region, n_features_tested=len(stats_df), n_features_significant=n_sig,
                top_feature=top["feature"], top_feature_p=top["p"], top_feature_d=top["cohens_d"])


# ── R-scores: Binder vs False Positive, per residue, ONE region at a time ─────
def r_score_region_analysis(df, region, out_dir, top_n=6):
    """Self-contained answer to: within `region` ALONE, which residues (if
    any) show a Binder vs False Positive R-score difference?"""
    rows = []
    for c in resid_columns(df):
        res = compare_groups(df, c)
        if res is None:
            continue
        rows.append(dict(resSeq=int(c.split("_")[1]), **res))
    if not rows:
        print(f"\n[R-score - {region}] Not enough Binder/False Positive data — skipping.")
        return None
    stats_df = pd.DataFrame(rows)
    stats_df["q"] = bh_fdr(stats_df["p"].values)
    stats_path = os.path.join(out_dir, f"r_score_{region}_binder_vs_fp_stats.csv")
    stats_df.sort_values("p").to_csv(stats_path, index=False)
    print(f"\nSaved {stats_path}")

    n_tested = len(stats_df)
    n_sig_raw = int((stats_df["p"] < 0.05).sum())
    n_sig_fdr = int((stats_df["q"] < 0.05).sum())
    print(f"[R-score - {region}] {n_tested} residues tested, {n_sig_raw} with p<0.05 "
          f"(uncorrected), {n_sig_fdr} with q<0.05 (FDR-corrected)")

    # ── Significance plot, this region only. Point color = which group is
    # higher (not region — there's only one region in this plot) ──
    fig, ax = plt.subplots(figsize=(max(10, stats_df["resSeq"].max() * 0.06), 5),
                            dpi=300, constrained_layout=True)
    colors = [GROUP_COLOR[GROUP_A] if d >= 0 else GROUP_COLOR[GROUP_B] for d in stats_df["cohens_d"]]
    ax.scatter(stats_df["resSeq"], -np.log10(stats_df["p"]), c=colors, s=25, alpha=0.85, zorder=2)
    sig = stats_df[stats_df["q"] < 0.05]
    if len(sig):
        sig_colors = [GROUP_COLOR[GROUP_A] if d >= 0 else GROUP_COLOR[GROUP_B] for d in sig["cohens_d"]]
        ax.scatter(sig["resSeq"], -np.log10(sig["p"]), facecolors="none", edgecolors=sig_colors,
                   s=110, linewidths=2, zorder=3)
        for _, row in sig.iterrows():
            ax.annotate(int(row["resSeq"]), (row["resSeq"], -np.log10(row["p"])),
                        textcoords="offset points", xytext=(0, 6), fontsize=7, ha="center")
    ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=1, label="p=0.05 (uncorrected)")
    ax.set_xlabel("Residue (resSeq)")
    ax.set_ylabel("-log10(p)")
    ax.set_title(f"{region.upper()} REGION ONLY — per-residue R-score, {GROUP_A} vs {GROUP_B}\n"
                 f"{n_sig_fdr}/{n_tested} residues significant (FDR q<0.05); "
                 "open circle = FDR-significant")
    ax.grid(True, alpha=0.4)
    legend_handles = [
        mpatches.Patch(facecolor=GROUP_COLOR[GROUP_A], label=f"{GROUP_A} higher"),
        mpatches.Patch(facecolor=GROUP_COLOR[GROUP_B], label=f"{GROUP_B} higher"),
    ]
    ax.legend(handles=legend_handles, fontsize=8)
    path = os.path.join(out_dir, f"r_score_{region}_significance.png")
    fig.savefig(path)
    plt.close(fig)
    print(f"Saved {path}")

    # ── Top-hit residues, this region only, one panel per residue ──
    top_hits = stats_df.sort_values("p").head(top_n)["resSeq"].tolist()
    print(f"[R-score - {region}] Top {len(top_hits)} residues by p-value: {top_hits}")
    if top_hits:
        n = len(top_hits)
        ncols = min(3, n)
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows),
                                  dpi=300, constrained_layout=True, squeeze=False)
        for i, resSeq in enumerate(top_hits):
            ax = axes[i // ncols][i % ncols]
            col = f"R_{resSeq}"
            box_jitter(ax, df, col, GROUP_ORDER)
            ax.set_xticks(range(len(GROUP_ORDER)))
            ax.set_xticklabels(GROUP_ORDER, rotation=20, ha="right", fontsize=8)
            ax.axhline(0, color="grey", linewidth=0.6, alpha=0.5)
            row = stats_df[stats_df.resSeq == resSeq].iloc[0]
            ax.set_title(f"Residue {resSeq}\np={row.p:.2g}, q={row.q:.2g}, d={row.cohens_d:.2f}", fontsize=9)
            ax.grid(True, alpha=0.4)
        for j in range(n, nrows * ncols):
            axes[j // ncols][j % ncols].axis("off")
        fig.suptitle(f"{region.upper()} REGION ONLY — top {n} residues, "
                     f"R-score by group ({GROUP_A} vs {GROUP_B})", fontsize=13)
        path = os.path.join(out_dir, f"r_score_{region}_top_residues.png")
        fig.savefig(path)
        plt.close(fig)
        print(f"Saved {path}")

    top = stats_df.sort_values("p").iloc[0]
    return dict(region=region, n_residues_tested=n_tested, n_residues_significant=n_sig_fdr,
                top_residue=int(top["resSeq"]), top_residue_p=top["p"], top_residue_d=top["cohens_d"])


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--r-scores-dir", default="/projects/ivta1597/biosensors/water_analysis/agg_out")
    parser.add_argument("--contact-dir", default="/projects/ivta1597/biosensors/LIG_contacts")
    parser.add_argument("--seq-list", default="/projects/ivta1597/biosensors/seq_ids_orig.txt")
    parser.add_argument("--tag", default="40_500ns")
    parser.add_argument("--structure-source", default=None,
                        help="Path to md_candidate_guide.csv. If given, restricts the analysis "
                             "to sequencing-confirmed sequences only (source == ngs_observed).")
    parser.add_argument("--out-dir", default="/projects/ivta1597/biosensors/analysis/core_vs_tail")
    parser.add_argument("--top-n", type=int, default=6,
                        help="Number of top-hit residues to plot individually (default: 6)")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    seq_type_map = load_seq_type_map(args.seq_list)

    print("=" * 70)
    print("Loading data (core/tail only)")
    print("=" * 70)
    r_data = load_r_scores(args.r_scores_dir, args.tag)
    contact_data = load_contact_features(args.contact_dir, args.tag, seq_type_map)

    if not r_data and not contact_data:
        raise FileNotFoundError(
            "No core/tail R-score or contact-feature CSVs found. Check "
            "--r-scores-dir / --contact-dir / --tag."
        )

    if args.structure_source:
        ngs_ids = load_ngs_observed_ids(args.structure_source)
        print(f"\nRestricting to sequencing-confirmed (source == ngs_observed) "
              f"sequences: {len(ngs_ids)} in {args.structure_source}")
        r_data = filter_to_ngs_observed(r_data, ngs_ids, "R-scores")
        contact_data = filter_to_ngs_observed(contact_data, ngs_ids, "Contact feats")

    contact_summaries, r_summaries = [], []
    for region in REGIONS:
        if region in contact_data:
            print("\n" + "=" * 70)
            print(f"CONTACT FEATURES — {region.upper()} region only: {GROUP_A} vs {GROUP_B}")
            print("=" * 70)
            s = contact_feature_region_analysis(contact_data[region], region, args.out_dir)
            if s:
                contact_summaries.append(s)

    for region in REGIONS:
        if region in r_data:
            print("\n" + "=" * 70)
            print(f"R-SCORE — {region.upper()} region only: {GROUP_A} vs {GROUP_B}, per residue")
            print("=" * 70)
            s = r_score_region_analysis(r_data[region], region, args.out_dir, top_n=args.top_n)
            if s:
                r_summaries.append(s)

    # ── Verdict: a plain table, not a chart — each region's result stands on
    # its own above; this just puts the headline numbers side by side ──
    print("\n" + "=" * 70)
    print(f"VERDICT: {GROUP_A} vs {GROUP_B} difference, by region")
    print("=" * 70)
    if contact_summaries:
        cdf = pd.DataFrame(contact_summaries)
        cdf.to_csv(os.path.join(args.out_dir, "verdict_contact_features.csv"), index=False)
        print("\nContact-type features:")
        print(cdf.to_string(index=False))
    if r_summaries:
        rdf = pd.DataFrame(r_summaries)
        rdf.to_csv(os.path.join(args.out_dir, "verdict_r_scores.csv"), index=False)
        print("\nPer-residue R-scores:")
        print(rdf.to_string(index=False))

    print("\nDone. All outputs written to:", args.out_dir)


if __name__ == "__main__":
    main()
