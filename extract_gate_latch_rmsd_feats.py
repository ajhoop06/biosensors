"""
extract_gate_latch_rmsd_feats.py

Computes gate and latch Ca RMSD-to-own-Boltz-reference (protein_<ID>.pdb),
independently for each region, for every sequence in seq_ids.txt.

This is a different axis from the existing gate-latch pairwise distance
(which conflates gate and latch motion into one relative coordinate) and
from per-region RMSF (which measures spread around a trajectory's own mean
structure, blind to where that mean sits). Here, each frame is superposed
on a "core" Ca selection (protein Ca atoms excluding gate, latch, Lb7a5,
and the C-terminal recoil helix) so gate/latch motion doesn't contaminate
the alignment, then gate and latch Ca RMSD are scored separately against
the per-sequence Boltz-predicted reference structure. This tests whether
each region's average position drifts toward or away from its modeled
starting pose over the trajectory.

Outputs:
    <run_dir>/gate_rmsd_to_ref.xvg    per-frame gate Ca RMSD to reference
    <run_dir>/latch_rmsd_to_ref.xvg   per-frame latch Ca RMSD to reference
    <REPO_DIR>/analysis/gate_latch_rmsd_to_ref_summary{TAG}.csv   per-sequence summary,
        including early/late window means (default 100 ns) and a
        full-trajectory linear regression slope (A/ns) for both gate and
        latch - the recommended ML features are the wide-window late mean
        and/or the regression slope, since those were the most robust
        binder/nonbinder discriminators found in this analysis.

Usage:
    python extract_gate_latch_rmsd_feats.py [seq_ids.txt] [--tag _500ns]
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import mdtraj as md
from scipy.stats import linregress

# ── Configurable paths / window tag ───────────────────────────────────────────
BASE     = "/Users/ivanatang/Library/CloudStorage/OneDrive-UCB-O365/Shirts Lab/LCA_boltz_models"
REPO_DIR = "/Users/ivanatang/Developer/biosensors"
RUNREL = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
TAG    = "_500ns"
NM_TO_ANG = 10.0
# ─────────────────────────────────────────────────────────────────────────────

TYPE_SUBDIR = {
    "Binder":         "binders",
    "False Positive": "nonbinders",
    "Low Confidence": "neg_low_pkt",
    "Fail Geometry":  "neg_fail_gate",
}

GATE   = (84, 90)
LATCH  = (114, 118)
LB7A5  = (148, 155)
RECOIL = (154, 166)

EXCLUDED_RESIDS = set()
for _lo, _hi in (GATE, LATCH, LB7A5, RECOIL):
    EXCLUDED_RESIDS.update(range(_lo, _hi + 1))


def seq_run_dir(seq_id, group_label):
    """Directory containing a sequence's per-frame trajectory outputs.
    Mirrors extract_rmsf_feats.py's rmsf_run_dir: newer pipeline runs write
    under runrel/500ns/, older ones write directly into runrel/."""
    run_dir = os.path.join(BASE, TYPE_SUBDIR[group_label], seq_id, RUNREL)
    nested_dir = os.path.join(run_dir, "500ns")
    if os.path.exists(os.path.join(nested_dir, "medoid_PL.pdb")):
        return nested_dir
    return run_dir


def find_reference_pdb(seq_id, group_label):
    """Per-sequence Boltz-predicted protein-only structure, protein_<ID>.pdb.
    Naming isn't fully uniform across naming conventions (e.g. protein_019.pdb,
    protein_pair3061.pdb, protein_seq10_binder.pdb), so glob rather than
    parse the ID out of seq_id."""
    seq_dir = os.path.join(BASE, TYPE_SUBDIR[group_label], seq_id)
    candidates = sorted(
        p for p in glob.glob(os.path.join(seq_dir, "protein_*.pdb"))
        if "fixed_H" not in os.path.basename(p)
    )
    return candidates[0] if candidates else None


def ca_index_map(topology):
    """resSeq -> atom index, for protein Ca atoms only."""
    return {
        res.resSeq: atom.index
        for res in topology.residues
        if res.is_protein
        for atom in res.atoms
        if atom.name == "CA"
    }


def region_indices(ca_map, lo, hi):
    return [ca_map[r] for r in range(lo, hi + 1) if r in ca_map]


def compute_region_rmsd(seq_id, group_label):
    ref_pdb = find_reference_pdb(seq_id, group_label)
    run_dir = seq_run_dir(seq_id, group_label)
    top_pdb = os.path.join(run_dir, "medoid_PL.pdb")
    xtc     = os.path.join(run_dir, "PL_only_40_500ns.xtc")

    if ref_pdb is None or not os.path.exists(top_pdb) or not os.path.exists(xtc):
        return None

    ref  = md.load(ref_pdb)
    traj = md.load(xtc, top=top_pdb)

    ref_ca  = ca_index_map(ref.topology)
    traj_ca = ca_index_map(traj.topology)

    core_resids = sorted((set(ref_ca) & set(traj_ca)) - EXCLUDED_RESIDS)
    if not core_resids:
        return None
    core_ref_idx  = [ref_ca[r]  for r in core_resids]
    core_traj_idx = [traj_ca[r] for r in core_resids]

    traj_sp = traj.superpose(ref, atom_indices=core_traj_idx, ref_atom_indices=core_ref_idx)

    def region_rmsd(lo, hi):
        ref_idx  = region_indices(ref_ca, lo, hi)
        traj_idx = region_indices(traj_ca, lo, hi)
        if not ref_idx or not traj_idx:
            return None
        diff = traj_sp.xyz[:, traj_idx, :] - ref.xyz[0, ref_idx, :]
        return np.sqrt(np.mean(np.sum(diff ** 2, axis=2), axis=1)) * NM_TO_ANG

    gate_rmsd  = region_rmsd(*GATE)
    latch_rmsd = region_rmsd(*LATCH)
    if gate_rmsd is None or latch_rmsd is None:
        return None

    time_ns = traj_sp.time / 1000.0  # mdtraj reports ps
    return time_ns, gate_rmsd, latch_rmsd


def write_xvg(path, time_ns, values, ylabel):
    with open(path, "w") as f:
        f.write('@    xaxis label "Time (ns)"\n')
        f.write(f'@    yaxis label "{ylabel}"\n')
        for t, v in zip(time_ns, values):
            f.write(f"{t:.4f}\t{v:.5f}\n")


def load_seq_ids(seq_list_path):
    """Yields (seq_id, group_label) tuples from a seq_ids.txt-style file."""
    with open(seq_list_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            seq_id      = parts[0].strip()
            group_label = parts[1].strip() if len(parts) > 1 else ""
            yield seq_id, group_label


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("seq_list", nargs="?", default="seq_ids.txt")
    parser.add_argument("--tag", default=TAG,
                        help=f"Suffix appended to output CSV filename (default: {TAG})")
    parser.add_argument("--wide-window-ns", type=float, default=100.0,
                        help="Window size (ns) for early/late mean comparison (default: 100.0), "
                             "used for the recommended ML feature.")
    args = parser.parse_args()
    tag = args.tag
    wide_window = args.wide_window_ns

    if not os.path.exists(args.seq_list):
        print(f"ERROR: seq list not found: {args.seq_list}")
        sys.exit(1)

    all_systems = list(load_seq_ids(args.seq_list))

    rows, missing = [], []
    for seq_id, group_label in all_systems:
        if group_label not in TYPE_SUBDIR:
            missing.append(seq_id)
            continue
        result = compute_region_rmsd(seq_id, group_label)
        if result is None:
            missing.append(seq_id)
            continue
        time_ns, gate_rmsd, latch_rmsd = result

        run_dir = seq_run_dir(seq_id, group_label)
        write_xvg(os.path.join(run_dir, "gate_rmsd_to_ref.xvg"), time_ns, gate_rmsd,
                  "Gate Ca RMSD to Boltz reference (A)")
        write_xvg(os.path.join(run_dir, "latch_rmsd_to_ref.xvg"), time_ns, latch_rmsd,
                  "Latch Ca RMSD to Boltz reference (A)")

        t_end = time_ns[-1]
        early_wide_mask = time_ns <= (time_ns[0] + wide_window)
        late_wide_mask  = time_ns >= (t_end - wide_window)

        row = {
            "Sequence": seq_id,
            "Group": group_label,
            "Gate RMSD mean (A)":        round(float(gate_rmsd.mean()), 4),
            "Gate RMSD SD (A)":          round(float(gate_rmsd.std(ddof=1)), 4),
            "Latch RMSD mean (A)":       round(float(latch_rmsd.mean()), 4),
            "Latch RMSD SD (A)":         round(float(latch_rmsd.std(ddof=1)), 4),
            "N frames": len(gate_rmsd),
        }

        for region_name, region_rmsd in (("Gate", gate_rmsd), ("Latch", latch_rmsd)):
            early_wide = float(region_rmsd[early_wide_mask].mean())
            late_wide  = float(region_rmsd[late_wide_mask].mean())
            slope      = float(linregress(time_ns, region_rmsd).slope)
            row[f"{region_name} RMSD early{int(wide_window)} mean (A)"] = round(early_wide, 4)
            row[f"{region_name} RMSD late{int(wide_window)} mean (A)"]  = round(late_wide, 4)
            row[f"{region_name} drift{int(wide_window)} (A)"]          = round(late_wide - early_wide, 4)
            row[f"{region_name} slope (A/ns)"]                          = round(slope, 6)

        rows.append(row)

    if missing:
        print(f"Skipped {len(missing)} sequences (missing reference/topology/trajectory): "
              f"{missing[:5]}{'...' if len(missing) > 5 else ''}")

    if not rows:
        print("\nNo sequences processed - nothing to write.")
        sys.exit(1)

    df = pd.DataFrame(rows).set_index("Sequence")
    out_dir = os.path.join(REPO_DIR, "analysis")
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, f"gate_latch_rmsd_to_ref_summary{tag}.csv")
    df.to_csv(csv_path)
    print(f"Saved ({len(df)} rows) -> {csv_path}")


if __name__ == "__main__":
    main()
