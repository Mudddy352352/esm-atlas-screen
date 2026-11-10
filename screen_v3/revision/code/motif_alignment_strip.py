#!/usr/bin/env python3
"""The structural alignment around the glucuronide-reading residues.

specificity_motif.py says which of the four residues each shortlisted
glucuronidase-fold protein keeps. It does not show why Arg562 is never one of
them. This writes the aligned stretch of each candidate against the E. coli
query over two windows, loop 6 around Tyr472 and loop 8 around Arg562, Asn566
and Lys568, with the alignment's own gaps preserved, so a figure can show that
Arg562 falls inside a deletion rather than being substituted.

The strip is taken from the same Foldseek alignment specificity_motif.py used,
and the four motif columns are checked against its table before anything is
written, so the two cannot disagree.

Usage: motif_alignment_strip.py
"""
import os

import pandas as pd

from side_chain_check import OUT, V1, V2, load_alignments

# E. coli beta-glucuronidase P05804 numbering, as in specificity_motif.py
MOTIF = {562: "R", 472: "Y", 566: "N", 568: "K"}
WINDOWS = [("loop 6", 466, 478), ("loop 8", 555, 572)]
GUS_QUERIES = ("P05804",)


def strip(a, lo, hi, offset):
    """Query and target characters over query positions lo..hi, gaps kept."""
    qi, ti = a["qStartPos"] - 1, a["dbStartPos"] - 1
    q, t = [], []
    for qc, tc in zip(a["qAln"], a["dbAln"]):
        # the query position this column sits at, in query-model numbering
        if qc != "-":
            pos = qi + 1 + offset - 1
            if lo <= pos <= hi:
                q.append(qc)
                t.append(tc)
            qi += 1
        if tc != "-":
            ti += 1
    return "".join(q), "".join(t)


def main():
    os.makedirs(OUT, exist_ok=True)
    t1 = pd.read_csv(V2 + "tier1_candidates.tsv", sep="\t")
    t1 = t1[t1["query"].str.startswith(GUS_QUERIES)].copy()
    aligns = load_alignments(set(zip(t1["query"], t1.target)))
    measured = pd.read_csv(OUT + "specificity_motif.tsv", sep="\t")

    rows, checked = [], 0
    for r in t1.itertuples():
        a = aligns.get((r.query, r.target))
        if a is None:
            continue
        offset = int(r.query.split("_w")[1].split("-")[0]) if "_w" in r.query else 1
        for name, lo, hi in WINDOWS:
            qs, ts = strip(a, lo, hi, offset)
            rows.append(dict(target=r.target, window=name, q_start=lo, q_end=hi,
                             query_seq=qs, candidate_seq=ts))
            # cross-check every motif residue that falls in this window
            for pos, expect in MOTIF.items():
                if not lo <= pos <= hi:
                    continue
                qsingle, tsingle = strip(a, pos, pos, offset)
                assert qsingle == expect, (
                    "query model does not carry %s%d, it carries %r"
                    % (expect, pos, qsingle))
                want = measured[(measured.target == r.target)
                                & (measured.motif_residue == "%s%d" % (expect, pos))]
                assert len(want) == 1, (r.target, pos)
                kept = bool(want.same_res.iloc[0])
                if (tsingle == expect) != kept:
                    raise SystemExit(
                        "strip disagrees with specificity_motif.tsv at %s %s%d: "
                        "alignment gives %r, table says kept=%s"
                        % (r.target, expect, pos, tsingle, kept))
                checked += 1

    df = pd.DataFrame(rows).sort_values(["window", "target"])
    out = OUT + "motif_alignment_strip.tsv"
    df.to_csv(out, sep="\t", index=False)
    print("%s  (%d rows, %d motif columns agree with specificity_motif.tsv)"
          % (out, len(df), checked))

    # the deletion each candidate carries where Arg562 sits
    l8 = df[df.window == "loop 8"]
    for r in l8.itertuples():
        gaps = sum(1 for c in r.candidate_seq if c == "-")
        print("  %s  loop 8 gaps %2d   %s" % (r.target, gaps, r.candidate_seq))


if __name__ == "__main__":
    main()
