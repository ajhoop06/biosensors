#!/usr/bin/env python3
"""
Computes the ligand hydration-shell water count for a single sequence: the
per-frame number of water oxygens within a first-hydration-shell distance of
any LCA heavy atom.

Rationale
---------
LCA is strongly hydrophobic. If binders settle the ligand into the pocket in
a way that excludes water (the classic hydrophobic-effect driving force for
binding), their hydration-shell counts should trend lower than nonbinders',
whose pockets may stay more solvent-exposed. Complements R_score_calc.py's
per-residue water-mediated-CONTACT score (which only counts water bridging a
specific residue-ligand contact) with a bulk measure of how much water sits
near the ligand overall, whether or not it's bridging any particular contact.

Cutoff: 3.5 A (the standard first-hydration-shell O...O radius), tighter than
R_score_calc.py's 4 A heavy-atom CONTACT cutoff (sized for protein side
chains, not water first-shell distance). No normalization is applied: LCA is
the same molecule in every sequence, so raw per-frame counts are already
comparable across sequences.

Output
------
  {out_dir}/{seq_id}_hydration_{TAG}.csv   per-frame time (ns) and hydration count

Usage
-----
    conda activate biosensors
    python hydration_calc.py --seq_id pair_3059_binder --seq_type binders
"""

import os
import warnings
import numpy as np
import pandas as pd
import mdtraj as md
import argparse

warnings.filterwarnings("ignore", category=DeprecationWarning)

parser = argparse.ArgumentParser()
parser.add_argument('--seq_id',    required=True)
parser.add_argument('--seq_type',  required=True)
parser.add_argument('--start-ns',  type=float, default=40.0,
                    help='Start of analysis window in ns (default: 40)')
parser.add_argument('--end-ns',    type=float, default=500.0,
                    help='End of analysis window in ns (default: 500)')
parser.add_argument('--hydration-cutoff', type=float, default=3.5,
                    help='First-hydration-shell distance cutoff in Angstrom (default: 3.5)')
parser.add_argument('--stride',    type=int, default=10,
                    help='Frame stride (default: 10; set to 1 for publication quality)')
args = parser.parse_args()
seq_id   = args.seq_id
seq_type = args.seq_type

# ── Analysis window ────────────────────────────────────────────────
START_NS = args.start_ns
END_NS   = args.end_ns
START_PS = int(START_NS * 1000)
END_PS   = int(END_NS   * 1000)
TAG      = f"{int(START_NS)}_{int(END_NS)}ns"   # e.g. "40_500ns"

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG  ← edit these before running
# ─────────────────────────────────────────────────────────────────────────────
# Trajectory inputs are read from the PetaLibrary archive, not scratch --
# scratch auto-deletes after 90 days and older runs' xtc/gro may already be
# gone, and this script needs the raw solvated trajectory (water retained),
# which the trimmed/protein-only scratch intermediates used elsewhere don't
# have. Mirrors R_score_calc.py's input_base.
input_base  = "/pl/active/shirts_archive/IvanaTang/biosensors"
output_base = "/scratch/alpine/ivta1597/LCA_boltz_models"
prod   = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"

traj_path = os.path.join(input_base, seq_type, seq_id, prod, "prod_md_500ns.xtc")
top_path  = os.path.join(input_base, seq_type, seq_id, prod, "prod_md_500ns.gro")

out_dir = os.path.join(output_base, seq_type, seq_id, f"water_density_{TAG}")
os.makedirs(out_dir, exist_ok=True)

LIG_RESNAME    = "LIG"
WATER_RESNAMES = {"HOH", "WAT", "SOL"}

HYDRATION_CUT_NM = args.hydration_cutoff / 10.0   # Angstrom -> nm
STRIDE           = args.stride


# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD TRAJECTORY
# ─────────────────────────────────────────────────────────────────────────────
traj = md.load(traj_path, top=top_path, stride=STRIDE)

mask = (traj.time >= START_PS) & (traj.time <= END_PS)
traj = traj[mask]
print(f"  Time window: {START_NS:.0f}-{END_NS:.0f} ns  "
      f"({mask.sum()} of {mask.shape[0]} frames retained after striding)")

top = traj.topology
nf  = traj.n_frames
print(f"  {nf} frames  |  {traj.n_atoms} atoms")

# ─────────────────────────────────────────────────────────────────────────────
# 2. PARSE TOPOLOGY — ligand heavy atoms and water oxygens
# ─────────────────────────────────────────────────────────────────────────────
all_res = list(top.residues)
lig_res = [r for r in all_res if r.name == LIG_RESNAME]
wat_res = [r for r in all_res if r.name in WATER_RESNAMES]

if not lig_res:
    raise ValueError(
        f"No residue named '{LIG_RESNAME}' found in topology.\n"
        "Check LIG_RESNAME in the CONFIG block.")

lig_heavy_atoms = [a for r in lig_res for a in r.atoms if a.element.symbol != 'H']
lig_heavy = np.array([a.index for a in lig_heavy_atoms], dtype=int)

# Water O atoms — require exactly 2 H atoms in the residue to confirm it's a
# real water molecule, not a stray crystallographic oxygen. Mirrors
# R_score_calc.py's water-oxygen extraction.
wat_O = []
for r in wat_res:
    O = [a.index for a in r.atoms if a.element.symbol == 'O']
    H = [a.index for a in r.atoms if a.element.symbol == 'H']
    if O and len(H) == 2:
        wat_O.append(O[0])
wat_O = np.array(wat_O, dtype=int)

print(f"  Ligand heavy atoms : {len(lig_heavy)}")
print(f"  Water oxygens      : {len(wat_O)}")

if len(wat_O) == 0:
    raise ValueError(
        f"No water oxygens found (checked WATER_RESNAMES={WATER_RESNAMES}).\n"
        "Check the water residue naming in this topology.")

# ─────────────────────────────────────────────────────────────────────────────
# 3. PER-FRAME HYDRATION COUNT
# ─────────────────────────────────────────────────────────────────────────────
lig_xyz  = traj.xyz[:, lig_heavy, :]   # (F, L, 3)
watO_xyz = traj.xyz[:, wat_O,     :]   # (F, W, 3)

hydration_count = np.zeros(nf, dtype=int)

print("Computing hydration-shell water counts...")
for f in range(nf):
    if f % 200 == 0:
        print(f"  frame {f:>5d}/{nf}", end='\r', flush=True)

    lp = lig_xyz[f]    # (L, 3)
    wp = watO_xyz[f]   # (W, 3)

    d = np.linalg.norm(wp[:, None, :] - lp[None, :, :], axis=-1)  # (W, L)
    hydration_count[f] = int((d.min(axis=1) < HYDRATION_CUT_NM).sum())

print(f"\n  Hydration count: mean={hydration_count.mean():.2f}  "
      f"std={hydration_count.std():.2f}  "
      f"min={hydration_count.min()}  max={hydration_count.max()}")

# ─────────────────────────────────────────────────────────────────────────────
# 4. SAVE CSV
# ─────────────────────────────────────────────────────────────────────────────
time_ns = traj.time / 1000.0   # mdtraj reports ps
df = pd.DataFrame({"time_ns": time_ns, "hydration_count": hydration_count})

csv_path = os.path.join(out_dir, f"{seq_id}_hydration_{TAG}.csv")
df.to_csv(csv_path, index=False)
print(f"Hydration count table saved -> {csv_path}")
