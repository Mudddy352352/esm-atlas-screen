#!/usr/bin/env python3
"""Does a candidate keep the residues that read the glucuronide's carboxylate?

Section 3.5 already says that the catalytic glutamate pair places a protein in
glycoside hydrolase family 2 and does not tell a glucuronidase from a
galactosidase. What distinguishes the two in GH2 is a small set of residues
that stabilise the carboxyl group at C6 of the sugar: the Asn-x-Lys motif in
loop 8, assisted by an arginine upstream in the same loop and a tyrosine in
loop 6 (Lebreton et al. 2025, and the structures cited there). In the
Escherichia coli enzyme used as the query, P05804, those are Asn566, Lys568,
Arg562 and Tyr472.

UniProt annotates only two of the four, Tyr472 and Lys568, so the screen's
active-site check never asked the specificity question separately: a candidate
could pass on the generic GH2 machinery alone. This scores all four, by residue
identity and by where the side chain actually sits, for the ten shortlisted
proteins found through the glucuronidase queries.

Usage: specificity_motif.py
"""
import os

import numpy as np
import pandas as pd

from side_chain_check import (FUNCTIONAL, OUT, V1, V2, candidate_path, centroid,
                              load_alignments, read_pdb, superpose)

# E. coli beta-glucuronidase P05804 numbering
MOTIF = [(562, "R", "loop 8 arginine, assists the carboxylate"),
         (472, "Y", "loop 6 tyrosine"),
         (566, "N", "NxK asparagine, loop 8"),
         (568, "K", "NxK lysine, loop 8")]
GUS_QUERIES = ("P05804",)


def main():
    os.makedirs(OUT, exist_ok=True)
    t1 = pd.read_csv(V2 + "tier1_candidates.tsv", sep="\t")
    t1 = t1[t1["query"].str.startswith(GUS_QUERIES)].copy()
    aligns = load_alignments(set(zip(t1["query"], t1.target)))

    rows = []
    for r in t1.itertuples():
        a = aligns.get((r.query, r.target))
        path = candidate_path(r.target)
        if a is None or path is None:
            continue
        offset = int(r.query.split("_w")[1].split("-")[0]) if "_w" in r.query else 1
        qres = read_pdb(os.path.join(V1, "structures", r.query + ".pdb"))
        tres = read_pdb(path)
        R, t, pairs = superpose(a, qres, tres)
        if R is None:
            continue
        for pos, expect, what in MOTIF:
            wpos = pos - offset + 1
            row = dict(target=r.target, motif_residue="%s%d" % (expect, pos),
                       role=what)
            if not (1 <= wpos <= len(qres)):
                row["note"] = "outside the query window"
                rows.append(row)
                continue
            q_one, q_at = qres[wpos - 1]
            if q_one != expect:
                row["note"] = "query model carries %s here, not %s" % (q_one, expect)
                rows.append(row)
                continue
            hit = pairs.get(wpos)
            if hit is None or hit[0] > len(tres):
                # how long is the unaligned stretch the residue falls in? A
                # residue missing because the query carries an insertion the
                # candidate does not have is a different fact from a residue
                # substituted, and loop 8 is where that happens here.
                lo = hi = wpos
                while lo - 1 >= 1 and pairs.get(lo - 1) is None:
                    lo -= 1
                while hi + 1 <= len(qres) and pairs.get(hi + 1) is None:
                    hi += 1
                row.update(candidate_res="-", same_res=False,
                           note="falls in a %d-residue stretch of the query "
                                "(%d-%d) with no counterpart in the candidate"
                                % (hi - lo + 1, lo + offset - 1, hi + offset - 1))
                rows.append(row)
                continue
            t_one, t_at = tres[hit[0] - 1]
            row["candidate_res"] = t_one
            row["same_res"] = (t_one == expect)
            qf = centroid(q_at, FUNCTIONAL.get(q_one, []))
            tf = centroid(t_at, FUNCTIONAL.get(t_one, []))
            if qf is not None and tf is not None:
                row["func_dist"] = round(float(np.linalg.norm(R @ tf + t - qf)), 2)
            rows.append(row)

    tab = pd.DataFrame(rows)
    for c in ["candidate_res", "same_res", "func_dist", "note"]:
        if c not in tab:
            tab[c] = np.nan
    tab = tab[["target", "motif_residue", "role", "candidate_res", "same_res",
               "func_dist", "note"]]
    tab.to_csv(OUT + "specificity_motif.tsv", sep="\t", index=False)
    print("The glucuronide-specificity residues in the ten GUS-fold shortlisted")
    print(tab.to_string(index=False))

    w = tab.pivot_table(index="target", columns="motif_residue",
                        values="same_res", aggfunc="first")
    w = w.reindex(columns=["R562", "Y472", "N566", "K568"])
    w["n_kept"] = w.sum(axis=1)
    print("\nwhich of the four each protein keeps")
    print(w.to_string())
    w.to_csv(OUT + "specificity_motif_wide.tsv", sep="\t")
    print("\nkept by all four: %d of %d proteins"
          % int((w.n_kept == 4).sum()) if False else
          "\nkept all four: %d of %d proteins" % (int((w.n_kept == 4).sum()), len(w)))
    print("kept three of four: %d; two: %d; one or none: %d"
          % (int((w.n_kept == 3).sum()), int((w.n_kept == 2).sum()),
             int((w.n_kept <= 1).sum())))
    for res in ["R562", "Y472", "N566", "K568"]:
        col = w[res].dropna()
        print("  %s kept by %d of %d" % (res, int(col.sum()), len(col)))
    gap = tab[tab.note.fillna("").str.contains("no counterpart")]
    if len(gap):
        print("\nmotif residues lost to a gap rather than a substitution:")
        for _, r in gap.iterrows():
            print("  %-18s %-5s %s" % (r.target, r.motif_residue, r.note))


if __name__ == "__main__":
    main()
