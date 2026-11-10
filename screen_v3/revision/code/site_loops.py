#!/usr/bin/env python3
"""Does any candidate carry an extra loop over its site, or the query one?

The second half of the loop question. Pocket enclosure says the candidates'
sites are not more open than the controls'; this says whether there is a loop
there at all. For each shortlisted protein it reads the Foldseek alignment and
lists every insertion of five residues or more, in either direction, whose
flanking aligned positions sit within 12 A of the site: a candidate-only
insertion is an extra segment the candidate could close over a substrate, and a
query-only insertion is a segment of the characterised enzyme the candidate does
not have.

Usage: site_loops.py
"""
import os

import numpy as np
import pandas as pd

from side_chain_check import V1, V2, load_alignments, read_pdb, superpose
from pocket_enclosure import DOCK

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
_ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

OUT = _ROOT + "/screen_v3/revision/results/"
MIN_LEN = 5
NEAR = 12.0


def runs(flags):
    """[(start, end)] of consecutive True in a 1-based dict-like flag list."""
    out, start = [], None
    for i, f in enumerate(flags):
        if f and start is None:
            start = i
        elif not f and start is not None:
            out.append((start, i - 1))
            start = None
    if start is not None:
        out.append((start, len(flags) - 1))
    return out


def main():
    t1 = pd.read_csv(V2 + "tier1_candidates.tsv", sep="\t")
    boxes = pd.read_csv(DOCK + "boxes_final.tsv", sep="\t").set_index("receptor")
    aligns = load_alignments(set(zip(t1["query"], t1.target)))

    rows = []
    for r in t1.itertuples():
        a = aligns.get((r.query, r.target))
        if a is None or r.target not in boxes.index:
            continue
        qres = read_pdb(os.path.join(V1, "structures", r.query + ".pdb"))
        from side_chain_check import candidate_path
        tres = read_pdb(candidate_path(r.target))
        R, t, pairs = superpose(a, qres, tres)
        if R is None:
            continue
        b = boxes.loc[r.target]
        site = np.array([b.cx, b.cy, b.cz])
        offset = int(r.query.split("_w")[1].split("-")[0]) if "_w" in r.query else 1

        # a query position is "covered" if the alignment pairs it with a residue
        qcov = [pairs.get(i) is not None for i in range(1, len(qres) + 1)]
        tpaired = {v[0] for v in pairs.values()}
        tcov = [i in tpaired for i in range(1, len(tres) + 1)]

        for label, cov, res, is_query in [("query only", qcov, qres, True),
                                          ("candidate only", tcov, tres, False)]:
            for s, e in runs([not c for c in cov]):
                n = e - s + 1
                if n < MIN_LEN:
                    continue
                # where is it, relative to the site, in the candidate's frame
                pts = [res[i][1]["CA"] for i in range(s, e + 1) if "CA" in res[i][1]]
                if not pts:
                    continue
                pts = np.array(pts)
                if is_query:
                    pts = np.array([R.T @ (p - t) for p in pts])  # query -> candidate
                d = float(np.linalg.norm(pts - site, axis=1).min())
                if d > NEAR:
                    continue
                rows.append(dict(target=r.target, query=r.query, side=label,
                                 length=n,
                                 span=("%d-%d" % (s + offset, e + offset)
                                       if is_query else "%d-%d" % (s + 1, e + 1)),
                                 min_dist_to_site=round(d, 1)))

    tab = pd.DataFrame(rows)
    tab.to_csv(OUT + "site_loops.tsv", sep="\t", index=False)
    print("Insertions of %d residues or more within %.0f A of the site"
          % (MIN_LEN, NEAR))
    if not len(tab):
        print("  none")
        return
    print(tab.sort_values(["side", "target"]).to_string(index=False))
    print("\nby side:")
    print(tab.groupby("side").agg(proteins=("target", "nunique"),
                                  segments=("length", "size"),
                                  longest=("length", "max")).to_string())
    gut = tab[tab.target.isin(["MGYP001533818318", "MGYP000245278657"])]
    print("\nthe two gut-sampled proteins:")
    print(gut.to_string(index=False) if len(gut) else "  no segment near the site")
    cand_only = set(tab[tab.side == "candidate only"].target)
    print("\ncandidates carrying an extra segment near their own site: %d of %d"
          % (len(cand_only), t1.target.nunique()))


if __name__ == "__main__":
    main()
