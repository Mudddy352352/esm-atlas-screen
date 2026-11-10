#!/usr/bin/env python3
"""Numbers the referees asked for: intervals, a test, an attrition count.

Covers reviewer 1's points 3 and 4 and two of the minor points. Everything is
recomputed from the screen's own result files; nothing is transcribed.

Usage: stats_for_revision.py
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
_ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

V1 = _ROOT + "/screen_v1/results/"
V2 = _ROOT + "/screen_v2/results/"
BENCH = _ROOT + "/screen_v3/benchmark/results/"
OUT = _ROOT + "/screen_v3/revision/results/"


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (100 * (c - h) / d, 100 * (c + h) / d)


def section_3_4():
    """83.9% against 42.3%, with the interval and the test it lacked.

    Reproduced from the joined per-hit table the screen itself wrote, using the
    same definitions rank_candidates.py used: structure-only means not found by
    phmmer at 1e-3 and below 30% identity, and the comparison is restricted to
    the query windows whose UniProt record annotates a catalytic residue."""
    j = pd.read_csv(V2 + "active_site_all_hits_joined.tsv.gz", sep="\t")
    j = j[j.n_catalytic > 0].copy()
    j["kept"] = j.catalytic_class == "catalytic_kept"
    j["structure_only"] = (~j.seq_detectable) & (j.seq_id_pct < 30.0)
    seq = j[j.seq_detectable]
    stronly = j[j.structure_only]
    ks, ns = int(seq.kept.sum()), len(seq)
    ko, no = int(stronly.kept.sum()), len(stronly)
    odds, p = fisher_exact([[ks, ns - ks], [ko, no - ko]])
    print("Section 3.4, recomputed from active_site_all_hits_joined.tsv.gz")
    print("  sequence-detectable %d/%d = %.1f%%  95%% CI %.1f to %.1f"
          % (ks, ns, 100 * ks / ns, *wilson(ks, ns)))
    print("  structure-only      %d/%d = %.1f%%  95%% CI %.1f to %.1f"
          % (ko, no, 100 * ko / no, *wilson(ko, no)))
    print("  Fisher exact two-sided p = %.3g, odds ratio %.1f" % (p, odds))
    print("  query windows contributing: %s" % ", ".join(sorted(j["query"].unique())))
    print("  per window (kept/n, sequence-detectable then structure-only):")
    rows = []
    for q, g in j.groupby("query"):
        s_, o_ = g[g.seq_detectable], g[g.structure_only]
        rows.append(dict(query=q, sd_kept=int(s_.kept.sum()), sd_n=len(s_),
                         so_kept=int(o_.kept.sum()), so_n=len(o_)))
        print("    %-26s %4d/%-4d  %3d/%-3d"
              % (q, rows[-1]["sd_kept"], rows[-1]["sd_n"],
                 rows[-1]["so_kept"], rows[-1]["so_n"]))
    pd.DataFrame(rows).to_csv(OUT + "section_3_4_by_query.tsv", sep="\t", index=False)
    return dict(seq_k=ks, seq_n=ns, so_k=ko, so_n=no, p=p, odds=odds)


def attrition():
    """How 30 catalytic-kept candidates become a shortlist of 13."""
    c = pd.read_csv(V2 + "active_site_candidates.tsv", sep="\t")
    t1 = pd.read_csv(V2 + "tier1_candidates.tsv", sep="\t")
    kept = c[(c.n_catalytic > 0) & (c.n_catalytic_identical == c.n_catalytic)]
    geo = kept[(kept.site_rmsd <= 2.0)
               & (kept.site_dist_mad.isna() | (kept.site_dist_mad <= 1.5))]
    tm = geo[geo.tm_score_qnorm >= 0.6]
    print("\nAttrition from the catalytic-kept set to the shortlist")
    print("  all catalytic residues kept            %d" % len(kept))
    print("  and within the geometric thresholds    %d" % len(geo))
    print("  and TM-score >= 0.6                    %d" % len(tm))
    print("  and phmmer-blind at any threshold      %d  (tier 1)" % len(t1))
    by_q = kept.groupby("query").size().sort_values(ascending=False)
    print("  catalytic-kept by query: %s"
          % ", ".join("%s %d" % (q.split("_")[0], n) for q, n in by_q.items()))
    return dict(kept=len(kept), geo=len(geo), tm=len(tm), tier1=len(t1))


def shortlist_plddt():
    """Model confidence for the thirteen, which the paper gives only as a sample."""
    t1 = pd.read_csv(V2 + "tier1_candidates.tsv", sep="\t")
    pl = pd.read_csv(V1 + "candidate_plddt.tsv", sep="\t")
    conf = pd.read_csv(V2 + "docking/pose_scan_with_confidence.tsv", sep="\t")
    conf = conf[conf.target.notna()][["target", "plddt_mean", "plddt_site"]]
    m = t1[["target"]].merge(conf, on="target", how="left").merge(
        pl[["target", "plddt_mean"]].rename(columns={"plddt_mean": "plddt_sample"}),
        on="target", how="left")
    m["plddt_mean"] = m.plddt_mean.fillna(m.plddt_sample)
    m = m.drop(columns=["plddt_sample"])
    m.to_csv(OUT + "shortlist_plddt.tsv", sep="\t", index=False)
    print("\nModel confidence for the thirteen shortlisted proteins")
    print(m.to_string(index=False))
    ok = m.plddt_mean.dropna()
    if len(ok):
        print("  whole-model pLDDT: median %.1f, range %.1f to %.1f, n = %d"
              % (ok.median(), ok.min(), ok.max(), len(ok)))
    site = m.plddt_site.dropna()
    if len(site):
        print("  local pLDDT at the site: median %.2f, range %.2f to %.2f, n = %d"
              % (site.median(), site.min(), site.max(), len(site)))
    return m


def benchmark_intervals():
    """Table 16 as it now stands, with the intervals the table itself carries.

    baseline_comparison.py computes the Wilson interval and the number of calls
    for every precision, so this only reads them back rather than reconstructing
    them from a rounded precision, which an earlier version of this script did.
    """
    tab = pd.read_csv(BENCH + "baseline_comparison.tsv", sep="\t")
    tab = tab[tab.precision.notna()]
    print("\nTable 16 with 95% intervals")
    for (ls, band), g in tab.groupby(["labelled_set", "band"], sort=False):
        print("  -- %s, identity %s (%d positives, %d negatives)"
              % (ls, band, g.n_pos.iloc[0], g.n_neg.iloc[0]))
        for _, r in g.iterrows():
            print("     %-36s %.3f  (%.3f to %.3f) over %d calls"
                  % (r.method, r.precision, r.ci_low, r.ci_high, r.n_called))
    low = tab[(tab.band == "0-30") & (tab.labelled_set == "name-labelled")]
    rule = low[low.method == "screen rule"].iloc[0]
    tm = low[low.method.str.startswith("TM-score")].iloc[0]
    print("\n  below 30% identity the screen rule and the TM-score baseline have"
          " overlapping intervals:")
    print("     screen rule %.3f (%.3f to %.3f), TM-score %.3f (%.3f to %.3f)"
          % (rule.precision, rule.ci_low, rule.ci_high,
             tm.precision, tm.ci_low, tm.ci_high))
    exp_rule = rule.precision * 10
    exp_tm = tm.precision * 10
    print("  expected true glucuronidases among ten shortlisted GUS-fold"
          " proteins: %.1f on the screen rule, %.1f on the TM-score baseline"
          % (exp_rule, exp_tm))
    sens = pd.read_csv(BENCH + "baseline_sensitivity.tsv", sep="\t")
    s030 = sens[(sens.band == "0-30") & (sens.method == "screen rule")]
    print("  screen rule below 30%% identity under four exclusion rules: %s"
          % ", ".join("%s %.3f" % (r.exclusion, r.precision)
                      for _, r in s030.iterrows()))
    return tab


if __name__ == "__main__":
    import os
    os.makedirs(OUT, exist_ok=True)
    section_3_4()
    attrition()
    shortlist_plddt()
    benchmark_intervals()
