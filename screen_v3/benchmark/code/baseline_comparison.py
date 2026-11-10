#!/usr/bin/env python3
"""What would a simpler method have done?

The retrieval evaluation reports AUROC per feature, which says how well each
feature ranks, and the screen rule's precision and recall, which says how the
pipeline actually behaves. A reviewer will ask the question those two leave
open: at the recall the screen achieves, is the screen's precision better than
simply thresholding sequence identity, or TM-score, or taking the Pfam family?

For every labelled set and identity band this fixes recall at the screen rule's
own value and reports each baseline's precision there, so the comparison is
like for like. The Pfam baseline has no threshold to tune, so it is reported at
its own operating point.

Every precision is reported with the number of calls it is computed over and a
95% Wilson interval, because several of these operating points rest on fewer
than twenty calls. The `--sensitivity` flag repeats the name-labelled
comparison under four different rules for which UniProt names carry no label,
so the reader can see how much the low-identity numbers depend on that choice.

Usage: baseline_comparison.py <benchmark_dir> [--sensitivity]
"""
import os
import sys

import numpy as np
import pandas as pd

GUS = "P05804"
GUS_NAME = "glucuronid"
UNANNOTATED = (r"deleted|uncharacteri|domain-containing|hypothetical"
               r"|^glycosi|^glycosyl|duf\d|^glycoside hydrolase")


def wilson(k, n, z=1.96):
    """95% Wilson score interval for k successes in n trials."""
    if not n:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def precision_at_recall(score, label, target_recall):
    """Highest precision any threshold on `score` reaches at >= target_recall.

    Returns the precision, the recall actually reached, the threshold, the
    number of calls the threshold makes and how many of them are true."""
    s = np.asarray(score, float)
    y = np.asarray(label, bool)
    ok = ~np.isnan(s)
    s, y = s[ok], y[ok]
    if not y.any() or not len(s) or target_recall != target_recall:
        return np.nan, np.nan, np.nan, 0, 0
    order = np.argsort(-s)
    y = y[order]
    tp = np.cumsum(y)
    fp = np.cumsum(~y)
    rec = tp / y.sum()
    prec = tp / (tp + fp)
    reach = rec >= target_recall - 1e-9
    if not reach.any():
        return np.nan, float(rec.max()), np.nan, 0, 0
    i = int(np.argmax(reach))
    return (float(prec[i]), float(rec[i]), float(s[order][i]),
            int(tp[i] + fp[i]), int(tp[i]))


def rule_scores(sub, col):
    tp = int((sub[col] & sub.label).sum())
    fp = int((sub[col] & ~sub.label).sum())
    fn = int((~sub[col] & sub.label).sum())
    prec = tp / (tp + fp) if tp + fp else np.nan
    recl = tp / (tp + fn) if tp + fn else np.nan
    return tp, fp, prec, recl


def ci_cols(k, n):
    lo, hi = wilson(k, n)
    return dict(n_called=int(n), n_true=int(k),
                ci_low=round(lo, 3) if lo == lo else np.nan,
                ci_high=round(hi, 3) if hi == hi else np.nan)


def compare(sub, set_name, band):
    tp, fp, prec, rec = rule_scores(sub, "screen_call")
    rows = [dict(labelled_set=set_name, band=band, n_pos=int(sub.label.sum()),
                 n_neg=int((~sub.label).sum()), method="screen rule",
                 precision=round(prec, 3) if prec == prec else np.nan,
                 recall=round(rec, 3) if rec == rec else np.nan,
                 threshold="", **ci_cols(tp, tp + fp),
                 note="catalytic residues kept, geometry ok, TM >= 0.6")]
    for name, col in [("sequence identity", "seq_id_pct"),
                      ("TM-score", "tm_score_qnorm"),
                      ("site retention", "frac_identical")]:
        p, r, thr, called, true = precision_at_recall(sub[col], sub.label, rec)
        rows.append(dict(labelled_set=set_name, band=band,
                         n_pos=int(sub.label.sum()), n_neg=int((~sub.label).sum()),
                         method=name + ", at matched recall",
                         precision=round(p, 3) if p == p else np.nan,
                         recall=round(r, 3) if r == r else np.nan,
                         threshold=round(thr, 3) if thr == thr else "",
                         **ci_cols(true, called),
                         note="threshold chosen to match the screen rule's recall"))
    tp, fp, prec, r = rule_scores(sub, "pfam_same_family")
    rows.append(dict(labelled_set=set_name, band=band, n_pos=int(sub.label.sum()),
                     n_neg=int((~sub.label).sum()), method="Pfam same family",
                     precision=round(prec, 3) if prec == prec else np.nan,
                     recall=round(r, 3) if r == r else np.nan,
                     threshold="", **ci_cols(tp, tp + fp), note="no threshold to tune"))
    return rows


def shares_pfam(row):
    if not isinstance(row.pfam, str) or not isinstance(row.query_pfam, str):
        return False
    return bool({x for x in row.pfam.split(";") if x}
                & {x for x in row.query_pfam.split(";") if x})


EXCLUSION_VARIANTS = {
    "as published": UNANNOTATED,
    "deleted records only": r"deleted",
    "no family-name exclusion": r"deleted|uncharacteri|domain-containing|hypothetical",
    "also drop putative and probable": UNANNOTATED + r"|putative|probable",
}


def build_name_set(res, exclusion=UNANNOTATED):
    """The name-labelled set, under a given rule for which names carry no label.

    UniProt entries whose record has been deleted, and entries whose name states
    no activity (uncharacterised, "glycoside hydrolase family 2", DUF and
    domain-containing names), carry no label either way: counting them as
    negatives is what made an earlier, unreproducible version of this analysis
    read far worse than the data support. They are excluded, not scored as
    negatives. `exclusion` is the regex that decides which, so the choice can be
    varied and its effect measured.
    """
    s = pd.read_csv(os.path.join(res, "labelled_hits_scored.tsv"), sep="\t")
    lab = pd.read_csv(os.path.join(res, "uniprot_labels.tsv"), sep="\t")
    ql = pd.read_csv(os.path.join(res, "query_labels.tsv"), sep="\t").rename(
        columns={"Entry": "acc", "EC number": "query_ec", "Pfam": "query_pfam"})
    s = s[(s.db == "afdb50") & (s.acc == GUS)].copy()
    s["uniprot"] = s.target.str.extract(r"AF-([A-Z0-9]+)-F")
    n = s.merge(lab[["uniprot", "protein_name", "pfam"]], on="uniprot", how="left")
    n = n.merge(ql[["acc", "query_pfam"]], on="acc", how="left")
    n = n.sort_values("tm_score_qnorm", ascending=False).drop_duplicates(["uniprot"])
    n = n[n.protein_name.notna()].copy()
    name = n.protein_name.str.lower()
    n = n[~name.str.contains(exclusion, regex=True)].copy()
    n["label"] = n.protein_name.str.lower().str.contains(GUS_NAME)
    n["catalytic_kept"] = np.where(
        n.n_catalytic > 0, n.n_catalytic_identical == n.n_catalytic, False)
    n["geometry_ok"] = (n.site_rmsd <= 2.0) & (
        n.site_dist_mad.isna() | (n.site_dist_mad <= 1.5))
    n["screen_call"] = (n.catalytic_kept.astype(bool) & n.geometry_ok
                        & (n.tm_score_qnorm >= 0.6))

    n["pfam_same_family"] = n.apply(shares_pfam, axis=1)
    return n


def sensitivity(res):
    """How much the name-labelled numbers depend on the exclusion rule."""
    rows = []
    for label, rx in EXCLUSION_VARIANTS.items():
        n = build_name_set(res, rx)
        for band, sub in [("all", n), ("0-30", n[n.seq_id_pct < 30]),
                          ("30+", n[n.seq_id_pct >= 30])]:
            if sub.label.sum() == 0:
                continue
            for r in compare(sub, "name-labelled", band):
                r["exclusion"] = label
                rows.append(r)
    tab = pd.DataFrame(rows)
    cols = ["exclusion", "band", "method", "n_pos", "n_neg", "precision",
            "recall", "ci_low", "ci_high", "n_called", "n_true"]
    tab[cols].to_csv(os.path.join(res, "baseline_sensitivity.tsv"),
                     sep="\t", index=False)
    for band in ["0-30", "all", "30+"]:
        g = tab[tab.band == band]
        if not len(g):
            continue
        print("\nidentity %s, name-labelled, under four exclusion rules" % band)
        print("%-32s %-36s %9s %7s %18s"
              % ("exclusion rule", "method", "precision", "calls", "95% CI"))
        for _, r in g.iterrows():
            ci = ("%.3f to %.3f" % (r.ci_low, r.ci_high)
                  if r.ci_low == r.ci_low else "")
            print("%-32s %-36s %9s %7s %18s"
                  % (r.exclusion, r.method, r.precision, r.n_called, ci))
    return tab


def main(bench):
    res = os.path.join(bench, "results")
    d = pd.read_csv(os.path.join(res, "retrieval_labelled_hits.tsv"), sep="\t")
    d = d[d.acc == GUS].copy()

    out = []
    for band, sub in [("all", d), ("0-30", d[d.seq_id_pct < 30]),
                      ("30+", d[d.seq_id_pct >= 30])]:
        if sub.label.sum() == 0:
            continue
        out += compare(sub, "EC-labelled", band)

    # The name-labelled set is wider than the EC-curated one, and the only one
    # with more than three positives below 30% identity. build_name_set says
    # which names are excluded and why.
    n = build_name_set(res)
    n.to_csv(os.path.join(res, "name_labelled_gus.tsv"), sep="\t", index=False)

    for band, sub in [("all", n), ("0-30", n[n.seq_id_pct < 30]),
                      ("30+", n[n.seq_id_pct >= 30])]:
        if sub.label.sum() == 0:
            continue
        out += compare(sub, "name-labelled", band)

    tab = pd.DataFrame(out)
    tab.to_csv(os.path.join(res, "baseline_comparison.tsv"), sep="\t", index=False)
    for (ls, band), g in tab.groupby(["labelled_set", "band"], sort=False):
        print("\n%s, identity %s  (%d positives, %d negatives)"
              % (ls, band, g.n_pos.iloc[0], g.n_neg.iloc[0]))
        print("%-34s %10s %8s %10s %18s %7s"
              % ("method", "precision", "recall", "threshold", "95% CI", "calls"))
        for _, r in g.iterrows():
            ci = ("%.3f to %.3f" % (r.ci_low, r.ci_high)
                  if r.ci_low == r.ci_low else "")
            print("%-34s %10s %8s %10s %18s %7s"
                  % (r.method, r.precision, r.recall, r.threshold, ci, r.n_called))


if __name__ == "__main__":
    args = sys.argv[1:]
    bench = args[0]
    main(bench)
    if "--sensitivity" in args:
        sensitivity(os.path.join(bench, "results"))
