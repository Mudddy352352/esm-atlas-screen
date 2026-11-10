#!/usr/bin/env python3
"""Is an open pocket, rather than a wrong site, why the gut proteins fail?

Section 3.9 reports that the two proteins sampled from the human large
intestine never place the glucuronide's sugar against their catalytic
carboxylates. A referee offers an alternative reading: several GH2
beta-glucuronidases close a loop over the substrate, the rigid predicted models
cannot do that, and a model whose pocket is left open would hold the ligand
badly whatever its catalytic residues look like. If that is what is happening,
it is not evidence about the protein.

The claim is testable. From the centre of each docking box, cast rays in many
directions and ask what fraction of them run into the protein: that is how
enclosed the site is, and an open, loop-less pocket should score low. The two
positive control crystals have their loops; the predicted models do not
necessarily. If the gut proteins are as enclosed as the controls, the open
pocket is not the explanation. If they are markedly more open, the referee is
right and the paper has to say so.

Usage: pocket_enclosure.py
"""
import os

import numpy as np
import pandas as pd

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
_ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

V1 = _ROOT + "/screen_v1/"
V2 = _ROOT + "/screen_v2/results/"
DOCK = V2 + "docking/"
FIG = _ROOT + "/manuscript/figures/"
XTAL = _ROOT + "/screen_v3/revision/structures/"
OUT = _ROOT + "/screen_v3/revision/results/"

N_RAYS = 1000
RAY_MAX = 15.0        # how far out a direction is followed, angstrom
RAY_TOL = 2.6         # how close an atom must be to the ray to block it
STEP = 0.5


def heavy_atoms(path, chain=None):
    xyz, plddt = [], []
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        if line[76:78].strip() == "H" or line[12:16].strip().startswith("H"):
            continue
        if chain and line[21] != chain:
            continue
        xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        try:
            plddt.append(float(line[60:66]))
        except ValueError:
            plddt.append(np.nan)
    return np.array(xyz), np.array(plddt)


def fibonacci(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi), np.cos(phi)])


def cavity_point(xyz, centre, reach=6.0, step=0.5):
    """The most open point near `centre`.

    The docking boxes were centred differently for models and for crystals: on
    the catalytic residues for a model, which puts the centre on the side chains
    themselves, and on the bound ligand for a crystal, which puts it in the
    cavity. Casting rays from those two kinds of point measures two different
    things. So the probe is moved to the point of greatest clearance within
    `reach` of the given centre, which is the same definition for every
    receptor, and the clearance is reported alongside."""
    g = np.arange(-reach, reach + step, step)
    pts = np.array([[x, y, z] for x in g for y in g for z in g])
    pts = centre + pts[np.linalg.norm(pts, axis=1) <= reach]
    near = xyz[np.linalg.norm(xyz - centre, axis=1) <= reach + 12.0]
    if not len(near):
        return centre, np.nan
    best, best_d = centre, -1.0
    for i in range(0, len(pts), 4000):
        chunk = pts[i:i + 4000]
        d = np.linalg.norm(chunk[:, None, :] - near[None, :, :], axis=-1).min(axis=1)
        k = int(np.argmax(d))
        if d[k] > best_d:
            best_d, best = float(d[k]), chunk[k]
    return best, best_d


def enclosure(xyz, centre):
    """Fraction of directions from `centre` that run into the protein."""
    dirs = fibonacci(N_RAYS)
    near = xyz[np.linalg.norm(xyz - centre, axis=1) <= RAY_MAX + RAY_TOL]
    if not len(near):
        return 0.0, np.nan
    rel = near - centre
    steps = np.arange(STEP, RAY_MAX + STEP, STEP)
    blocked = np.zeros(len(dirs), bool)
    first = np.full(len(dirs), np.nan)
    for k, s in enumerate(steps):
        pts = dirs * s
        d = np.linalg.norm(pts[:, None, :] - rel[None, :, :], axis=-1).min(axis=1)
        hit = (d <= RAY_TOL) & ~blocked
        first[hit] = s
        blocked |= hit
        if blocked.all():
            break
    return float(blocked.mean()), float(np.nanmedian(first))


def receptor_path(name):
    cand = [os.path.join(V1, "atlas_structures", name + ".pdb"),
            DOCK + "atlas_structures_extra/" + name + ".pdb",
            os.path.join(V1, "structures", name + ".pdb"),
            FIG + "pose_work/rec_" + name + ".pdb",
            XTAL + name.split("_")[0] + ".pdb"]
    for p in cand:
        if os.path.exists(p):
            return p
    return None


def main():
    os.makedirs(OUT, exist_ok=True)
    boxes = pd.read_csv(DOCK + "boxes_final.tsv", sep="\t")
    conf = pd.read_csv(DOCK + "pose_scan_with_confidence.tsv", sep="\t")
    final = pd.read_csv(DOCK + "final_table.tsv", sep="\t")
    frac = dict(zip(final.receptor, final.frac_poses_catalytic))

    rows = []
    for b in boxes.itertuples():
        if b.kind == "neg_offsite":
            continue
        path = receptor_path(b.receptor)
        if path is None:
            rows.append(dict(receptor=b.receptor, kind=b.kind,
                             note="structure not on disk"))
            continue
        chain = None
        if b.receptor.endswith(("_A",)) and "MGYP" not in b.receptor:
            chain = "A"
        xyz, plddt = heavy_atoms(path, chain)
        box = np.array([b.cx, b.cy, b.cz])
        centre, clearance = cavity_point(xyz, box)
        d_near = float(np.linalg.norm(xyz - box, axis=1).min())
        enc, depth = enclosure(xyz, centre)
        shell = np.linalg.norm(xyz - centre, axis=1) <= 12.0
        rows.append(dict(receptor=b.receptor, kind=b.kind,
                         enclosure=round(enc, 3),
                         clearance=round(clearance, 2),
                         probe_moved=round(float(np.linalg.norm(centre - box)), 1),
                         median_wall_dist=round(depth, 1) if depth == depth else np.nan,
                         atoms_within_12A=int(shell.sum()),
                         nearest_atom=round(d_near, 2),
                         frac_poses_catalytic=frac.get(b.receptor),
                         structure=os.path.basename(path)))
    tab = pd.DataFrame(rows)
    pl = conf[conf.target.notna()][["target", "plddt_site", "plddt_mean"]]
    tab = tab.merge(pl.rename(columns={"target": "receptor"}), on="receptor",
                    how="left")
    tab = tab.sort_values(["kind", "enclosure"], ascending=[True, False])
    tab.to_csv(OUT + "pocket_enclosure.tsv", sep="\t", index=False)
    print("How enclosed each docking site is, %d rays, %.0f A, %.1f A tolerance"
          % (N_RAYS, RAY_MAX, RAY_TOL))
    print(tab.to_string(index=False))

    pos = tab[tab.kind.isin(["pos_xtal", "pos_pred"])].enclosure.dropna()
    gut = tab[tab.receptor.isin(["MGYP001533818318", "MGYP000245278657"])]
    cand = tab[tab.kind == "candidate"].enclosure.dropna()
    print("\npositive controls: enclosure %.3f to %.3f (median %.3f)"
          % (pos.min(), pos.max(), pos.median()))
    print("all candidates:    enclosure %.3f to %.3f (median %.3f)"
          % (cand.min(), cand.max(), cand.median()))
    for _, r in gut.iterrows():
        print("  %s (human large intestine): enclosure %.3f, %d atoms within 12 A, "
              "pose fraction %s" % (r.receptor, r.enclosure, r.atoms_within_12A,
                                    r.frac_poses_catalytic))
    ok = tab[tab.enclosure.notna() & tab.frac_poses_catalytic.notna()]
    if len(ok) > 3:
        from scipy.stats import spearmanr
        rho, p = spearmanr(ok.enclosure, ok.frac_poses_catalytic)
        print("\nenclosure against the fraction of catalytically placed poses: "
              "Spearman rho = %.2f, p = %.3f, n = %d" % (rho, p, len(ok)))


if __name__ == "__main__":
    main()
