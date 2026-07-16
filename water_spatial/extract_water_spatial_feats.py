"""
extract_water_spatial_feats.py

Parses each sequence's water_density_{seq_id}.cube (a Gaussian cube file
produced by `gmx spatial` — a 3D spatial density function of water-oxygen
occupancy, in the water_spatial_prep.sh/water_spatial_run.sh fitted-trajectory
reference frame) and computes a per-sequence scalar pocket water density: the
mean voxel density within a fixed radius of the ligand's mean position.

Ligand centroid is computed directly from this pipeline's own
fit_trim.xtc + fit_trim_ref.pdb (mean LIG heavy-atom position across the
fitted trajectory) rather than from the separate RMSD pipeline's
medoid_PL.pdb, to avoid depending on two independently-computed least-squares
fits staying bit-identical -- everything needed to place the region of
interest comes from this pipeline's own output.

Cube file format (Gaussian cube, as written by gmx spatial):
    Line 1-2   : comments
    Line 3     : NATOMS, origin_x, origin_y, origin_z   (bohr)
    Line 4-6   : NX/NY/NZ voxel counts + per-axis step vectors (bohr)
    Next NATOMS lines : atom records (ignored here)
    Remaining  : NX*NY*NZ density values, whitespace-separated across
        however many values-per-line gmx spatial happens to write (observed:
        one full Z-row per line, i.e. NZ values/line, not the 6-per-line
        convention some cube writers use) -- parse_cube() flattens all
        remaining whitespace-separated tokens regardless of line breaks, so
        this doesn't matter for correctness, only worth knowing if you're
        eyeballing the raw file.

1 bohr = 0.529177 Angstrom. Voxel step vectors are assumed axis-aligned
(gmx spatial writes an orthogonal grid); off-diagonal step components are
ignored. This whole parser (units, atom-record skip count, Z-fastest/Y/X
reshape order) was verified against a real GROMACS 2025.3 `gmx spatial`
run (bind_019_binder): computed mean/min/max of the parsed density array
matched gmx spatial's own reported "Raw data: average ..., min ..., max ..."
line exactly.

water_spatial_run.sh runs `gmx spatial -nodiv`, so voxel values here are raw
per-frame occupancy counts (not already normalized to bulk-water density --
`-div` turned out to be a boolean visualization-only flag, not a numeric
divisor, see water_spatial_run.sh). This script does the actual bulk-water
normalization below, in Angstrom^3 units for consistency with the rest of
the repo (everything else here uses Angstrom, not nm).

Output: water_spatial_feats.csv

Usage:
    python extract_water_spatial_feats.py [seq_ids_ngs_observed.txt]
                                           [--pocket-cutoff 8.0] [--out water_spatial_feats.csv]
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import mdtraj as md

BOHR_TO_ANG = 0.529177210903

# Bulk liquid water molecular number density, ~0.0334 waters/Angstrom^3
# (equivalent to the commonly-cited ~33.3 nm^-3 / ~1 g/cm^3 for water).
# Used only to express pocket_water_density_mean as a "fraction of bulk"
# for interpretability -- verify against this system's actual water model
# if that distinction matters (TIP3P vs SPC/E etc. differ slightly).
BULK_WATER_DENSITY_PER_ANG3 = 0.0334

TYPE_SUBDIR = {
    "Binder":         "binders",
    "False Positive": "nonbinders",
    "Low Confidence": "neg_low_pkt",
    "Fail Geometry":  "neg_fail_gate",
}

BASE   = "/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"


def parse_cube(cube_path):
    """Returns (origin_ang (3,), step_ang (3,), dims (3,) int, density (NX,NY,NZ))."""
    with open(cube_path) as f:
        lines = f.readlines()

    header = lines[2].split()
    natoms = int(header[0])
    origin_bohr = np.array([float(x) for x in header[1:4]])

    dims = np.zeros(3, dtype=int)
    step_bohr = np.zeros(3)
    for axis in range(3):
        parts = lines[3 + axis].split()
        dims[axis] = int(parts[0])
        # Only the diagonal step component is used -- see module docstring
        # caveat about the axis-aligned-grid assumption.
        step_bohr[axis] = float(parts[1 + axis])

    data_start = 6 + natoms
    values = []
    for line in lines[data_start:]:
        values.extend(float(x) for x in line.split())

    density = np.array(values).reshape(dims)  # Z fastest already matches this shape/order

    origin_ang = origin_bohr * BOHR_TO_ANG
    step_ang   = step_bohr * BOHR_TO_ANG
    return origin_ang, step_ang, dims, density


def ligand_centroid_ang(ref_pdb, fit_trim_xtc, lig_resname="LIG"):
    """Mean position (Angstrom) of the ligand's heavy atoms across the fitted,
    production-windowed trajectory -- the same reference frame the cube grid
    was generated in."""
    traj = md.load(fit_trim_xtc, top=ref_pdb)
    top  = traj.topology
    lig_atoms = [a.index for r in top.residues if r.name == lig_resname
                 for a in r.atoms if a.element.symbol != 'H']
    if not lig_atoms:
        raise ValueError(f"No residue named '{lig_resname}' found in {ref_pdb}")
    com_per_frame = traj.xyz[:, lig_atoms, :].mean(axis=1)   # (F, 3), nm
    return com_per_frame.mean(axis=0) * 10.0                  # -> Angstrom


def pocket_density(origin_ang, step_ang, density, centroid_ang, cutoff_ang):
    """Mean raw per-frame occupancy count within `cutoff_ang` of
    `centroid_ang`, plus the voxel volume (Angstrom^3) needed to convert
    that into an actual number density."""
    nx, ny, nz = density.shape
    ix, iy, iz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    voxel_xyz = origin_ang + np.stack([ix, iy, iz], axis=-1) * step_ang
    dist = np.linalg.norm(voxel_xyz - centroid_ang, axis=-1)
    mask = dist <= cutoff_ang
    n_voxels = int(mask.sum())
    voxel_volume_ang3 = float(step_ang[0] * step_ang[1] * step_ang[2])
    if n_voxels == 0:
        return np.nan, 0, voxel_volume_ang3
    return float(density[mask].mean()), n_voxels, voxel_volume_ang3


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("seq_list", nargs="?",
                        default="/projects/ivta1597/biosensors/seq_ids_ngs_observed.txt")
    parser.add_argument("--pocket-cutoff", type=float, default=8.0,
                        help="Radius in Angstrom around the ligand centroid to average "
                             "voxel density over (default: 8.0, matching "
                             "pkt_vol/select_binding_pocket.py's pocket-sphere cutoff).")
    parser.add_argument("--base", default=BASE)
    parser.add_argument("--output", default="water_spatial_feats.csv")
    args = parser.parse_args()

    if not os.path.exists(args.seq_list):
        print(f"ERROR: seq list not found: {args.seq_list}")
        sys.exit(1)

    records = []
    missing = []

    with open(args.seq_list) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts       = line.split("\t")
            seq_id      = parts[0].strip()
            seq_type    = parts[1].strip() if len(parts) > 1 else ""
            custom_path = parts[2].strip() if len(parts) > 2 else ""

            if custom_path:
                run_dir = os.path.join(custom_path, RUNREL, "water_spatial")
            else:
                dir_type = TYPE_SUBDIR.get(seq_type, seq_type)
                run_dir  = os.path.join(args.base, dir_type, seq_id, RUNREL, "water_spatial")

            cube_path = os.path.join(run_dir, f"water_density_{seq_id}.cube")
            ref_pdb   = os.path.join(run_dir, "fit_trim_ref.pdb")
            fit_xtc   = os.path.join(run_dir, "fit_trim.xtc")

            if not (os.path.exists(cube_path) and os.path.exists(ref_pdb)
                    and os.path.exists(fit_xtc)):
                print(f"MISSING: {seq_id}  [{seq_type}]  ->  {run_dir}")
                missing.append(seq_id)
                continue

            try:
                origin_ang, step_ang, dims, density = parse_cube(cube_path)
                centroid_ang = ligand_centroid_ang(ref_pdb, fit_xtc)
                mean_count, n_voxels, voxel_vol_ang3 = pocket_density(
                    origin_ang, step_ang, density, centroid_ang, args.pocket_cutoff)

                # mean_count is a raw per-frame occupancy count per voxel
                # (gmx spatial -nodiv). Convert to an actual number density
                # (waters/Angstrom^3) by dividing by voxel volume, then
                # express as a fraction of bulk for interpretability.
                density_per_ang3 = (mean_count / voxel_vol_ang3
                                     if voxel_vol_ang3 > 0 else np.nan)
                frac_bulk = (density_per_ang3 / BULK_WATER_DENSITY_PER_ANG3
                             if not np.isnan(density_per_ang3) else np.nan)

                records.append({
                    "seq_id":   seq_id,
                    "seq_type": seq_type,
                    "pocket_water_density_mean": mean_count,
                    "pocket_water_density_per_ang3": density_per_ang3,
                    "pocket_water_density_frac_bulk": frac_bulk,
                    "pocket_water_density_n_voxels":  n_voxels,
                    "grid_dims": f"{dims[0]}x{dims[1]}x{dims[2]}",
                })
                print(f"OK: {seq_id}  [{seq_type}]  "
                      f"mean_count={mean_count:.4f}  frac_bulk={frac_bulk:.4f}  "
                      f"n_voxels={n_voxels}")
            except Exception as e:
                print(f"ERROR: {seq_id}  -  {e}")
                missing.append(seq_id)

    if not records:
        print("\nNo sequences loaded - nothing to write.")
        sys.exit(1)

    feat_df = pd.DataFrame(records)
    feat_df["pocket_water_density_missing"] = False
    feat_df.to_csv(args.output, index=False)

    print(f"\nFeatures written to: {args.output}")
    print(f"  Sequences processed : {len(records)}")
    print(f"  Sequences missing   : {len(missing)}")
    if missing:
        print(f"  Missing seq_ids     : {', '.join(missing)}")


if __name__ == "__main__":
    main()
