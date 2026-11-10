#!/usr/bin/env python3
"""Does the shortlist survive a different choice of threshold?

The screen's three geometric cut-offs (TM-score, site RMSD, the rotation-
invariant MAD) were set once, by eye, on the query set. A referee asked the
obvious question: how much of the result is the data and how much is the
cut-off? This sweeps each threshold one at a time and together, and reports
both what it does to the shortlist and what it does to the precision measured
on the labelled benchmark, so the two move in view of each other.

Usage: threshold_sensitivity.py
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, _ROOT + "/screen_v3/benchmark/code")
from baseline_comparison import build_name_set, rule_scores, wilson  # noqa: E402

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
_ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

V2 = _ROOT + "/screen_v2/results/"
BENCH = _ROOT + "/screen_v3/benchmark/results/"
OUT = _ROOT + "/screen_v3/revision/results/"

TM = [0.5, 0.55, 0.6, 0.65, 0.7]
RMSD = [1.5, 1.75, 2.0, 2.5, 3.0]
MAD = [1.0, 1.25, 1.5, 2.0, 2.5]
PUBLISHED = (0.6, 2.0, 1.5)


def shortlist(c, tm, rmsd, mad):
    """The tier-1 rule as rank_candidates.py applies it, with the cut-offs free."""
    kept = (c.n_catalytic > 0) & (c.n_catalytic_identical == c.n_catalytic)
    geo = (c.site_rmsd <= rmsd) & (c.site_dist_mad.isna() | (c.site_dist_mad <= mad))
    return c[kept & geo & (c.tm_score_qnorm >= tm) & c.phmmer_evalue.isna()]


def bench_precision(n, tm, rmsd, mad, band):
    sub = n if band == "all" else (n[n.seq_id_pct < 30] if band == "0-30"
                                  else n[n.seq_id_pct >= 30])
    kept = np.where(sub.n_catalytic > 0,
                    sub.n_catalytic_identical == sub.n_catalytic, False)
    geo = (sub.site_rmsd <= rmsd) & (sub.site_dist_mad.isna()
                                     | (sub.site_dist_mad <= mad))
    call = pd.Series(kept.astype(bool), index=sub.index) & geo & (
        sub.tm_score_qnorm >= tm)
    s = sub.assign(screen_call=call)
    tp, fp, prec, rec = rule_scores(s, "screen_call")
    lo, hi = wilson(tp, tp + fp)
    return dict(precision=prec, recall=rec, n_called=tp + fp, n_true=tp,
                ci_low=lo, ci_high=hi)


def main():
    os.makedirs(OUT, exist_ok=True)
    c = pd.read_csv(V2 + "active_site_candidates.tsv", sep="\t")
    n = build_name_set(BENCH)

    rows = []
    grids = [("TM-score", TM, 0), ("site RMSD", RMSD, 1), ("site MAD", MAD, 2)]
    for name, values, slot in grids:
        for v in values:
            thr = list(PUBLISHED)
            thr[slot] = v
            sl = shortlist(c, *thr)
            r = dict(varied=name, value=v, tm=thr[0], rmsd=thr[1], mad=thr[2],
                     shortlist=len(sl), published=(tuple(thr) == PUBLISHED))
            for band in ["all", "0-30"]:
                b = bench_precision(n, *thr, band)
                r["prec_" + band] = (round(b["precision"], 3)
                                     if b["precision"] == b["precision"] else np.nan)
                r["recall_" + band] = (round(b["recall"], 3)
                                       if b["recall"] == b["recall"] else np.nan)
                r["calls_" + band] = b["n_called"]
            rows.append(r)

    # all three moved together, strict and loose
    for label, thr in [("all three strict", (0.7, 1.5, 1.0)),
                       ("all three loose", (0.5, 3.0, 2.5))]:
        sl = shortlist(c, *thr)
        r = dict(varied=label, value=np.nan, tm=thr[0], rmsd=thr[1], mad=thr[2],
                 shortlist=len(sl), published=False)
        for band in ["all", "0-30"]:
            b = bench_precision(n, *thr, band)
            r["prec_" + band] = (round(b["precision"], 3)
                                 if b["precision"] == b["precision"] else np.nan)
            r["recall_" + band] = (round(b["recall"], 3)
                                   if b["recall"] == b["recall"] else np.nan)
            r["calls_" + band] = b["n_called"]
        rows.append(r)

    tab = pd.DataFrame(rows)
    tab.to_csv(OUT + "threshold_sensitivity.tsv", sep="\t", index=False)
    print("Shortlist size and benchmark precision against each cut-off")
    print("%-18s %6s %5s %5s %5s %10s %9s %8s %10s %8s"
          % ("varied", "value", "TM", "RMSD", "MAD", "shortlist",
             "prec all", "rec all", "prec 0-30", "rec 0-30"))
    for _, r in tab.iterrows():
        mark = "  <- published" if r.published else ""
        print("%-18s %6s %5s %5s %5s %10d %9s %8s %10s %8s%s"
              % (r.varied, r.value, r.tm, r.rmsd, r.mad, r.shortlist,
                 r.prec_all, r.recall_all, r["prec_0-30"], r["recall_0-30"], mark))

    # which members of the shortlist are stable across the whole grid
    import itertools
    seen = {}
    combos = list(itertools.product(TM, RMSD, MAD))
    for thr in combos:
        for t in shortlist(c, *thr).target:
            seen[t] = seen.get(t, 0) + 1
    stab = pd.DataFrame(sorted(seen.items(), key=lambda kv: -kv[1]),
                        columns=["target", "n_threshold_sets"])
    stab["frac"] = (stab.n_threshold_sets / len(combos)).round(3)
    stab.to_csv(OUT + "shortlist_stability.tsv", sep="\t", index=False)
    print("\nHow often each protein is shortlisted across all %d threshold sets"
          % len(combos))
    print(stab.to_string(index=False))
    pub = set(shortlist(c, *PUBLISHED).target)
    always = set(stab[stab.frac == 1.0].target)
    print("\n  shortlisted at the published thresholds: %d" % len(pub))
    print("  shortlisted under every threshold set:    %d" % len(always))
    print("  published members that are not always in: %s"
          % ", ".join(sorted(pub - always)) or "none")


if __name__ == "__main__":
    main()
