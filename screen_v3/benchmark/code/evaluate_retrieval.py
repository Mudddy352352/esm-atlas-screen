#!/usr/bin/env python3
"""Does the screen's ranking and its active-site filter actually separate
enzymes that do the query's chemistry from same-fold enzymes that do not?

The Atlas carries no functional labels, so the screen could never answer this.
The same searches also hit AlphaFold DB, and about 1,700 of those hits carry a
UniProt EC number, which gives a labelled retrieval problem of exactly the kind
the pipeline claims to solve: among proteins that share the query's fold, find
the ones with the query's activity.

Positives are hits whose EC matches the query's EC (at the query's own level of
specificity, since several queries are annotated 1.7.1.- rather than to four
digits). Negatives are hits with some other EC, most of them enzymes from the
same family acting on a different sugar or substrate, which is the hard case.
Hits with no EC are not scored either way.

Usage: evaluate_retrieval.py <scored_tsv> <uniprot_labels_tsv> <query_labels_tsv> <out_dir>
"""
import os
import sys

import numpy as np
import pandas as pd

RMSD_OK, MAD_OK = 2.0, 1.5


def auroc(pos, neg):
    """Mann-Whitney U as the area under the ROC curve, ties counted as half."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan, len(pos), len(neg)
    allv = np.concatenate([pos, neg])
    order = allv.argsort()
    ranks = np.empty(len(allv), float)
    ranks[order] = np.arange(1, len(allv) + 1)
    # average ranks over ties
    s = pd.Series(allv)
    ranks = s.rank(method="average").to_numpy()
    u = ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2.0
    return u / (len(pos) * len(neg)), len(pos), len(neg)


def ec_match(query_ec, hit_ec):
    """A hit counts as a positive when it carries an EC that matches the
    query's at every level the query itself specifies."""
    if not isinstance(hit_ec, str) or not isinstance(query_ec, str):
        return False
    for q in query_ec.split(";"):
        q = q.strip()
        if not q:
            continue
        qp = q.split(".")
        for h in hit_ec.split(";"):
            hp = h.strip().split(".")
            if len(hp) != 4:
                continue
            ok = True
            for a, b in zip(qp, hp):
                if a == "-":
                    break
                if a != b:
                    ok = False
                    break
            if ok:
                return True
    return False


def main(scored_tsv, labels_tsv, qlabels_tsv, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    s = pd.read_csv(scored_tsv, sep="\t")
    lab = pd.read_csv(labels_tsv, sep="\t")
    ql = pd.read_csv(qlabels_tsv, sep="\t").rename(
        columns={"Entry": "acc", "EC number": "query_ec", "Pfam": "query_pfam"})

    s = s[s.db == "afdb50"].copy()
    s["uniprot"] = s.target.str.extract(r"AF-([A-Z0-9]+)-F")
    d = s.merge(lab[["uniprot", "ec", "pfam", "protein_name", "reviewed"]],
                on="uniprot", how="left").merge(
                    ql[["acc", "query_ec", "query_pfam"]], on="acc", how="left")

    # one row per (query enzyme, target): keep the window that aligned best
    d = d.sort_values("tm_score_qnorm", ascending=False).drop_duplicates(["acc", "uniprot"])
    d = d[d.ec.notna()].copy()
    d["label"] = [ec_match(q, h) for q, h in zip(d.query_ec, d.ec)]

    # the screen's own operating rule, written out
    d["catalytic_kept"] = np.where(
        d.n_catalytic > 0, d.n_catalytic_identical == d.n_catalytic, np.nan)
    d["geometry_ok"] = (d.site_rmsd <= RMSD_OK) & (
        d.site_dist_mad.isna() | (d.site_dist_mad <= MAD_OK))
    d["screen_call"] = (d.catalytic_kept.fillna(False).astype(bool)
                        & d.geometry_ok & (d.tm_score_qnorm >= 0.6))
    # Pfam baseline: does the hit carry any Pfam domain of the query?
    def shares_pfam(row):
        if not isinstance(row.pfam, str) or not isinstance(row.query_pfam, str):
            return False
        a = {x for x in row.pfam.split(";") if x}
        b = {x for x in row.query_pfam.split(";") if x}
        return bool(a & b)
    d["pfam_same_family"] = d.apply(shares_pfam, axis=1)

    rows = []
    for acc, sub in d.groupby("acc"):
        pos, neg = sub[sub.label], sub[~sub.label]
        rec = dict(query=acc, query_ec=sub.query_ec.iloc[0],
                   n_labelled=len(sub), n_positive=len(pos), n_negative=len(neg))
        for name, col in [("tm_score", "tm_score_qnorm"),
                          ("seq_identity", "seq_id_pct"),
                          ("site_retention", "frac_identical")]:
            a, np_, nn = auroc(pos[col], neg[col])
            rec["auroc_" + name] = round(a, 3) if a == a else np.nan
        for name, col in [("screen_rule", "screen_call"),
                          ("pfam_same_family", "pfam_same_family")]:
            tp = int((sub[col] & sub.label).sum())
            fp = int((sub[col] & ~sub.label).sum())
            fn = int((~sub[col] & sub.label).sum())
            prec = tp / (tp + fp) if tp + fp else np.nan
            recl = tp / (tp + fn) if tp + fn else np.nan
            rec[name + "_tp"] = tp
            rec[name + "_fp"] = fp
            rec[name + "_precision"] = round(prec, 3) if prec == prec else np.nan
            rec[name + "_recall"] = round(recl, 3) if recl == recl else np.nan
        rows.append(rec)

    out = pd.DataFrame(rows).sort_values("n_positive", ascending=False)
    out.to_csv(os.path.join(out_dir, "retrieval_by_query.tsv"), sep="\t", index=False)
    d.to_csv(os.path.join(out_dir, "retrieval_labelled_hits.tsv"), sep="\t", index=False)
    cols = ["query", "n_positive", "n_negative", "auroc_tm_score", "auroc_seq_identity",
            "auroc_site_retention", "screen_rule_tp", "screen_rule_fp",
            "screen_rule_precision", "screen_rule_recall",
            "pfam_same_family_precision", "pfam_same_family_recall"]
    print(out[cols].to_string(index=False))


if __name__ == "__main__":
    main(*sys.argv[1:5])
