#!/usr/bin/env python3
"""Pool the M-CSA scan into the generalised retention statement.

The screen's side result was: of same-fold hits a single-sequence search also
finds, 83.9% keep their catalytic residues, against 42.3% of the hits only
structure search reaches. That rested on nine query enzymes. This asks the same
question of every M-CSA family scanned, so the number becomes a statement about
folds in general rather than about nine proteins.

Usage: analyse_retention.py <work_dir>
"""
import sys

import numpy as np
import pandas as pd

BANDS = [(0, 30), (30, 50), (50, 100.1)]


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def rate(df):
    n = len(df)
    k = int(df["catalytic_kept"].sum())
    lo, hi = wilson(k, n)
    return n, k, 100.0 * k / n if n else np.nan, 100.0 * lo, 100.0 * hi


def main(work):
    res = work.rstrip("/") + "/results/"
    df = pd.read_csv(res + "mcsa_hits.tsv", sep="\t")
    df["catalytic_kept"] = df["catalytic_kept"].astype(str).str.lower().eq("true")
    df["seq_detectable"] = df["seq_detectable"].astype(str).str.lower().eq("true")
    df = df.drop_duplicates(["mcsa_id", "target"])

    # One row per family, so a huge family cannot carry the pooled number.
    per = []
    for (mid, name), g in df.groupby(["mcsa_id", "enzyme"]):
        d = g[g.seq_detectable]
        s = g[~g.seq_detectable]
        per.append(dict(
            mcsa_id=mid, enzyme=name, n_catalytic=int(g.n_catalytic.iloc[0]),
            n_hits=len(g), n_seq=len(d), n_structure_only=len(s),
            pct_kept_seq=round(100.0 * d.catalytic_kept.mean(), 1) if len(d) else np.nan,
            pct_kept_structure_only=round(100.0 * s.catalytic_kept.mean(), 1) if len(s) else np.nan,
            median_seq_id_structure_only=round(s.seq_id_pct.median(), 1) if len(s) else np.nan))
    per = pd.DataFrame(per).sort_values("mcsa_id")
    per.to_csv(res + "retention_by_family.tsv", sep="\t", index=False)

    lines = []
    def say(s=""):
        lines.append(s)
        print(s)

    say("Families scored: %d   hits scored: %d" % (per.mcsa_id.nunique(), len(df)))
    say()

    say("Pooled over all hits")
    say("%-24s %8s %8s %8s %18s" % ("group", "hits", "kept", "pct", "95% CI"))
    for label, sub in [("sequence-detectable", df[df.seq_detectable]),
                       ("structure-only", df[~df.seq_detectable]),
                       ("all same-fold hits", df)]:
        n, k, p, lo, hi = rate(sub)
        say("%-24s %8d %8d %7.1f%% %8.1f - %4.1f%%" % (label, n, k, p, lo, hi))
    say()

    both = per.dropna(subset=["pct_kept_seq", "pct_kept_structure_only"])
    diff = both.pct_kept_seq - both.pct_kept_structure_only
    say("Per family, both groups present: %d families" % len(both))
    say("  median retention, sequence-detectable : %.1f%%" % both.pct_kept_seq.median())
    say("  median retention, structure-only      : %.1f%%" % both.pct_kept_structure_only.median())
    say("  median within-family difference       : %+.1f points" % diff.median())
    say("  families where sequence-detectable is higher: %d of %d"
        % (int((diff > 0).sum()), len(both)))
    if len(both) >= 6:
        from scipy.stats import wilcoxon
        stat, p = wilcoxon(both.pct_kept_seq, both.pct_kept_structure_only)
        say("  Wilcoxon signed-rank on the paired family rates: W=%.1f, p=%.2g" % (stat, p))
    say()

    say("By sequence identity to the query, pooled")
    say("%-14s %8s %8s %8s %18s" % ("identity band", "hits", "kept", "pct", "95% CI"))
    rows = []
    for lo_b, hi_b in BANDS:
        sub = df[(df.seq_id_pct >= lo_b) & (df.seq_id_pct < hi_b)]
        n, k, p, lo, hi = rate(sub)
        say("%-14s %8d %8d %7.1f%% %8.1f - %4.1f%%"
            % ("%d-%d%%" % (lo_b, min(hi_b, 100)), n, k, p, lo, hi))
        rows.append(dict(band="%d-%d" % (lo_b, min(hi_b, 100)), n_hits=n,
                         n_kept=k, pct_kept=round(p, 1),
                         ci_low=round(lo, 1), ci_high=round(hi, 1)))
    pd.DataFrame(rows).to_csv(res + "retention_by_identity_band.tsv", sep="\t", index=False)
    say()

    say()
    say("The confound: structure-only hits are remote by construction, so the gap")
    say("above could be nothing but a sequence-identity effect. Within band:")
    say("%-14s %22s %22s" % ("identity band", "sequence-detectable", "structure-only"))
    rows = []
    for lo_b, hi_b in BANDS:
        sub = df[(df.seq_id_pct >= lo_b) & (df.seq_id_pct < hi_b)]
        a, b = sub[sub.seq_detectable], sub[~sub.seq_detectable]
        na, ka, pa, loa, hia = rate(a)
        nb, kb, pb, lob, hib = rate(b)
        say("%-14s %9d hits %6.1f%% %9d hits %6.1f%%"
            % ("%d-%d%%" % (lo_b, min(hi_b, 100)), na, pa, nb, pb))
        rows.append(dict(band="%d-%d" % (lo_b, min(hi_b, 100)),
                         n_seq=na, pct_kept_seq=round(pa, 1),
                         n_structure_only=nb, pct_kept_structure_only=round(pb, 1)))
    pd.DataFrame(rows).to_csv(res + "retention_band_by_detectability.tsv",
                              sep="\t", index=False)

    say("Retention against site size (more curated residues is a harder test)")
    say("%-14s %10s %8s %8s" % ("n catalytic", "families", "hits", "pct kept"))
    for size, g in df.groupby("n_catalytic"):
        if len(g) < 50:
            continue
        say("%-14d %10d %8d %7.1f%%"
            % (size, g.mcsa_id.nunique(), len(g), 100.0 * g.catalytic_kept.mean()))

    open(res + "RETENTION_SUMMARY.txt", "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
