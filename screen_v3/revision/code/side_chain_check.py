#!/usr/bin/env python3
"""The active-site test asked where the backbone is. This asks about the side chains.

Sections 2.4 and 3.3 score a site residue on three things: that the alignment
puts the same amino acid there, that its CA lands close to the query's CA after
superposition, and that the site residues sit the same distances apart. All
three are backbone or register questions. A referee is right that the text reads
as though catalytic geometry had been checked, and it had not: two glutamates
can have their CAs 1 A apart and point their carboxylates in opposite
directions.

This fills the gap for the thirteen shortlisted proteins, where full-atom
models are on disk. For each curated site residue it superposes the candidate
on the query using the Foldseek alignment, then measures

  * the distance between the functional atoms of the two side chains
    (carboxylate midpoint for Asp and Glu, OH for Tyr, NZ for Lys, ring
    centroid for Trp, and so on), and
  * the chi-1 torsion of each, and their difference.

Usage: side_chain_check.py
"""
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, _ROOT + "/screen_v2/code")
from active_site_check import THREE_TO_ONE, best_fit, walk  # noqa: E402

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
_ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

V1 = _ROOT + "/screen_v1/"
V2 = _ROOT + "/screen_v2/results/"
EXTRA = V2 + "docking/atlas_structures_extra/"
OUT = _ROOT + "/screen_v3/revision/results/"

# the atoms that do the chemistry, by residue
FUNCTIONAL = {
    "D": ["OD1", "OD2"], "E": ["OE1", "OE2"], "N": ["OD1", "ND2"],
    "Q": ["OE1", "NE2"], "Y": ["OH"], "K": ["NZ"], "R": ["NE", "NH1", "NH2"],
    "H": ["ND1", "NE2"], "S": ["OG"], "T": ["OG1"], "C": ["SG"],
    "W": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "F": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "M": ["SD"], "G": ["CA"], "A": ["CB"], "P": ["CG"],
    "L": ["CD1", "CD2"], "I": ["CD1"], "V": ["CG1", "CG2"],
}
# chi-1 is N-CA-CB-X, X being the first side-chain atom after CB
CHI1_LAST = {
    "D": "CG", "E": "CG", "N": "CG", "Q": "CG", "Y": "CG", "F": "CG",
    "W": "CG", "H": "CG", "K": "CG", "R": "CG", "L": "CG", "M": "CG",
    "P": "CG", "S": "OG", "T": "OG1", "C": "SG", "I": "CG1", "V": "CG1",
}


def read_pdb(path):
    """residue index (1-based, in file order) -> (one-letter, {atom: xyz})."""
    res, order = {}, []
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        key = (line[21], line[22:27])
        if key not in res:
            res[key] = [THREE_TO_ONE.get(line[17:20].strip(), "X"), {}]
            order.append(key)
        res[key][1][line[12:16].strip()] = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return [res[k] for k in order]


def centroid(atoms, names):
    pts = [atoms[n] for n in names if n in atoms]
    return np.mean(pts, axis=0) if pts else None


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def chi1(one, atoms):
    last = CHI1_LAST.get(one)
    if not last:
        return np.nan
    need = ["N", "CA", "CB", last]
    if any(n not in atoms for n in need):
        return np.nan
    return dihedral(*[atoms[n] for n in need])


def angdiff(a, b):
    if a != a or b != b:
        return np.nan
    d = abs(a - b) % 360.0
    return min(d, 360.0 - d)


def candidate_path(target):
    for p in (os.path.join(V1, "atlas_structures", target + ".pdb"),
              EXTRA + target + ".pdb"):
        if os.path.exists(p):
            return p
    return None


def load_alignments(want):
    """(query window, target) -> the Foldseek alignment record for that pair."""
    queries = {q for q, _ in want}
    aligns = {}
    for jf in sorted(glob.glob(os.path.join(V1, "foldseek_raw", "*.json"))):
        qname = os.path.basename(jf).replace(".pdb.json", "")
        if qname not in queries:
            continue
        d = json.load(open(jf))
        for res in d["results"]:
            if res["db"] != "mgnify_esm30":
                continue
            al = res["alignments"]
            al = al[0] if al and isinstance(al[0], list) else al
            for a in al:
                tgt = str(a["target"]).split()[0].replace(".pdb.gz", "")
                if (qname, tgt) in want:
                    aligns[(qname, tgt)] = a
    return aligns


def superpose(a, qres, tres):
    """Rotation and translation taking the candidate onto the query model."""
    qxyz = np.array([x[1]["CA"] for x in qres if "CA" in x[1]])
    txyz = np.array([x[1]["CA"] for x in tres if "CA" in x[1]])
    pairs = walk(a)
    qi = [i for i in pairs if i <= len(qxyz) and pairs[i][0] <= len(txyz)]
    if len(qi) < 3:
        return None, None, pairs
    R, t, _ = best_fit(np.array([qxyz[i - 1] for i in qi]),
                       np.array([txyz[pairs[i][0] - 1] for i in qi]),
                       a["qLen"])
    return R, t, pairs


def main():
    os.makedirs(OUT, exist_ok=True)
    t1 = pd.read_csv(V2 + "tier1_candidates.tsv", sep="\t")
    sites = pd.read_csv(V2 + "query_sites.tsv", sep="\t")
    aligns = load_alignments(set(zip(t1["query"], t1.target)))

    rows = []
    for r in t1.itertuples():
        a = aligns.get((r.query, r.target))
        path = candidate_path(r.target)
        if a is None or path is None:
            rows.append(dict(target=r.target, query=r.query, site="",
                             note="alignment or model not on disk"))
            continue
        acc = r.query.split("_")[0]
        m = re.search(r"_w(\d+)-(\d+)$", r.query)
        offset = int(m.group(1)) if m else 1
        qres = read_pdb(os.path.join(V1, "structures", r.query + ".pdb"))
        tres = read_pdb(path)
        qxyz = np.array([x[1]["CA"] for x in qres if "CA" in x[1]])
        txyz = np.array([x[1]["CA"] for x in tres if "CA" in x[1]])

        pairs = walk(a)
        qi = [i for i in pairs if i <= len(qxyz) and pairs[i][0] <= len(txyz)]
        if len(qi) < 3:
            continue
        R, t, _ = best_fit(np.array([qxyz[i - 1] for i in qi]),
                           np.array([txyz[pairs[i][0] - 1] for i in qi]),
                           a["qLen"])
        if R is None:
            continue

        s = sites[sites.acc == acc].copy()
        s["wpos"] = s.pos - offset + 1
        s = s[(s.wpos >= 1) & (s.wpos <= len(qres))].sort_values("wpos")
        for sr in s.itertuples():
            hit = pairs.get(int(sr.wpos))
            tag = "%s%d" % (sr.res, int(sr.pos))
            row = dict(target=r.target, query=r.query, site=tag,
                       kind=sr.kind, catalytic=sr.kind in ("Active site",
                                                           "Modified residue"))
            if hit is None or hit[0] > len(tres):
                row.update(note="not aligned")
                rows.append(row)
                continue
            tpos, _ = hit
            q_one, q_at = qres[int(sr.wpos) - 1]
            t_one, t_at = tres[tpos - 1]
            row["candidate_res"] = t_one
            row["same_res"] = (t_one == sr.res)
            qca, tca = q_at.get("CA"), t_at.get("CA")
            if qca is not None and tca is not None:
                row["ca_dist"] = round(float(np.linalg.norm(
                    R @ tca + t - qca)), 2)
            qf = centroid(q_at, FUNCTIONAL.get(q_one, []))
            tf = centroid(t_at, FUNCTIONAL.get(t_one, []))
            if qf is not None and tf is not None:
                row["func_dist"] = round(float(np.linalg.norm(
                    R @ tf + t - qf)), 2)
            row["chi1_query"] = round(chi1(q_one, q_at), 1)
            row["chi1_candidate"] = round(chi1(t_one, t_at), 1)
            row["chi1_diff"] = round(angdiff(row["chi1_query"],
                                             row["chi1_candidate"]), 1)
            rows.append(row)

    tab = pd.DataFrame(rows)
    cols = ["target", "query", "site", "kind", "catalytic", "candidate_res",
            "same_res", "ca_dist", "func_dist", "chi1_query",
            "chi1_candidate", "chi1_diff", "note"]
    for c in cols:
        if c not in tab:
            tab[c] = np.nan
    tab = tab[cols]
    tab.to_csv(OUT + "side_chain_check.tsv", sep="\t", index=False)

    print("Side-chain geometry at the curated site residues of the thirteen")
    print(tab.to_string(index=False))

    ok = tab[tab.func_dist.notna()]
    cat = ok[ok.catalytic.fillna(False).astype(bool)]
    print("\nfunctional-atom distance, all site residues: median %.2f A, "
          "range %.2f to %.2f (n = %d)"
          % (ok.func_dist.median(), ok.func_dist.min(), ok.func_dist.max(), len(ok)))
    if len(cat):
        print("catalytic residues only:                    median %.2f A, "
              "range %.2f to %.2f (n = %d)"
              % (cat.func_dist.median(), cat.func_dist.min(),
                 cat.func_dist.max(), len(cat)))
    print("\nhow far the side chain moves while the CA stays close:")
    both = ok[ok.ca_dist.notna()]
    print("  CA within 2 A but functional atoms more than 3 A apart: %d of %d"
          % (int(((both.ca_dist <= 2) & (both.func_dist > 3)).sum()), len(both)))
    print("  chi-1 more than 60 degrees apart: %d of %d"
          % (int((tab.chi1_diff > 60).sum()), int(tab.chi1_diff.notna().sum())))
    per = ok.groupby("target").func_dist.agg(["median", "max", "count"])
    print("\nby protein (functional-atom distance, A):")
    print(per.round(2).to_string())
    per.round(2).to_csv(OUT + "side_chain_by_protein.tsv", sep="\t")


if __name__ == "__main__":
    main()
