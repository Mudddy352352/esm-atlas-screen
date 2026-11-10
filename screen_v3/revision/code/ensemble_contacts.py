#!/usr/bin/env python3
"""The contacts of Figures 9 and 10, over every pose instead of the best one.

Sections 2.11 and 3.9 read the contacts of a single pose per receptor, the
best-scoring one, and both figure captions say so. A referee is right that one
pose from a rigid docking run is a thin thing to read chemistry off, and that
the honest version is a frequency over the ensemble: in what fraction of the
27 poses does a given residue touch the ligand at all.

This runs the same contact code over every pose of every seed for the four
receptors in those figures, and reports, per residue, how often each kind of
contact occurs. Nothing about the method changes; only the number of poses it
is asked about.

Usage: ensemble_contacts.py
"""
import os
import sys

import numpy as np
import pandas as pd

FIG = _ROOT + "/manuscript/figures"
sys.path.insert(0, FIG)
import ligand_interactions as li      # noqa: E402
import make_pose_figures as mpf       # noqa: E402

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
_ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

OUT = _ROOT + "/screen_v3/revision/results/"
SCRATCH = "/tmp/ensemble_contacts"
SEEDS = [7, 999331, 20261021]
MODES = range(1, 10)


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(SCRATCH, exist_ok=True)
    mpf.fetch_6d4o()
    chem, rings = li.ligand_chemistry()
    sugar = mpf.sugar_serials()
    all_poses = pd.read_csv(
        _ROOT + "/screen_v2/results/docking/pose_scan_all_poses.tsv",
        sep="\t")

    rows, per_pose = [], []
    for p in mpf.PANELS:
        key = p["key"]
        rec = os.path.join(SCRATCH, "rec_%s.pdb" % key)
        mpf.write_receptor_pdb(p["pdb"], p["chain"], rec)
        cat = set(p["cat"])
        seen = {}
        n_pose = 0
        for seed in SEEDS:
            src = (_ROOT + "/screen_v2/results/docking/poses/"
                   "%s__SN38G__s%s.pdbqt" % (key, seed))
            if not os.path.exists(src):
                continue
            for mode in MODES:
                dest = os.path.join(SCRATCH, "pose.pdb")
                try:
                    xyz = mpf.write_pose_pdb(key, seed, mode, dest)
                except SystemExit:
                    continue
                r = li.interactions(rec, xyz, chem, rings, sugar)
                n_pose += 1
                kinds = set()
                for c in r["contacts"]:
                    k = (c["res_name"], c["res_num"], c["kind"])
                    e = seen.setdefault(k, dict(poses=set(), best=c["dist"],
                                                n_contacts=0, parts=set()))
                    # one pose can make several atom-level contacts of the same
                    # kind to the same residue; the frequency counts poses, not
                    # contacts.
                    e["poses"].add((seed, mode))
                    e["n_contacts"] += 1
                    e["best"] = min(e["best"], c["dist"])
                    e["parts"].add(c["lig_part"])
                    kinds.add((c["res_num"], c["kind"]))
                cat_touch = any(rn in cat for rn, _ in kinds)
                row = all_poses[(all_poses.receptor == key)
                                & (all_poses.seed == seed)
                                & (all_poses["mode"] == mode)]
                per_pose.append(dict(
                    receptor=key, seed=seed, mode=mode,
                    affinity=float(row.affinity.iloc[0]) if len(row) else np.nan,
                    sugar_to_catalytic=(float(row.sugar_to_catalytic.iloc[0])
                                        if len(row) else np.nan),
                    catalytic_contact=cat_touch,
                    n_contacts=len(r["contacts"]),
                    n_pocket_residues=len(r["residues"])))
        for (name, num, kind), e in seen.items():
            rows.append(dict(receptor=key, residue=li.res_label(name, num),
                             res_num=num, kind=li.KIND_LABEL.get(kind, kind),
                             catalytic=num in cat,
                             n_poses_with=len(e["poses"]), n_poses=n_pose,
                             frac_poses=round(len(e["poses"]) / n_pose, 2),
                             contacts_total=e["n_contacts"],
                             best_dist=round(e["best"], 2),
                             ligand_parts=";".join(sorted(e["parts"]))))

    tab = pd.DataFrame(rows).sort_values(
        ["receptor", "frac_poses", "res_num"], ascending=[True, False, True])
    tab.to_csv(OUT + "ensemble_contacts.tsv", sep="\t", index=False)
    pp = pd.DataFrame(per_pose)
    pp.to_csv(OUT + "ensemble_contacts_per_pose.tsv", sep="\t", index=False)

    for key, g in tab.groupby("receptor", sort=False):
        n = int(g.n_poses.iloc[0])
        print("\n%s, %d poses" % (key, n))
        print("  residues contacted in at least half the poses:")
        for _, r in g[g.frac_poses >= 0.5].iterrows():
            print("    %-8s %-22s %4.0f%% of poses, closest %.2f A  [%s]%s"
                  % (r.residue, r.kind, 100 * r.frac_poses, r.best_dist,
                     r.ligand_parts, "  <- catalytic" if r.catalytic else ""))
        cat = g[g.catalytic]
        if len(cat):
            print("  catalytic residues, any contact kind:")
            for _, r in cat.iterrows():
                print("    %-8s %-22s %4.0f%% of poses, closest %.2f A"
                      % (r.residue, r.kind, 100 * r.frac_poses, r.best_dist))
        else:
            print("  neither catalytic residue contacts the ligand in any pose")
        sub = pp[pp.receptor == key]
        print("  poses with a contact to either catalytic residue: %d of %d (%.0f%%)"
              % (int(sub.catalytic_contact.sum()), len(sub),
                 100 * sub.catalytic_contact.mean()))
        print("  contacts per pose: median %d, range %d to %d"
              % (sub.n_contacts.median(), sub.n_contacts.min(),
                 sub.n_contacts.max()))


if __name__ == "__main__":
    main()
