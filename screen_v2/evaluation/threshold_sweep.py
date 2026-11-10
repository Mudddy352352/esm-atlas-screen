#!/usr/bin/env python3
"""Threshold sweep for the tier-1 shortlist (manuscript section 3.10, Supplementary Note S4).

How much of the shortlist is the data and how much is the cut-off? The three
geometric thresholds were set once, by eye, on the query set, so the length of
the shortlist is partly a choice. This sweeps them and reports how many
proteins survive each setting, and which survive all of them.

Runs offline from two deposited files:
  screen_v1/results/candidates.tsv          TM-score and phmmer E-value per pair
  screen_v2/results/active_site_all_hits.tsv  active-site geometry per pair

and checks itself against screen_v2/results/tier1_candidates.tsv.

Usage:  python3 screen_v2/evaluation/threshold_sweep.py [--out DIR]
"""
import argparse, collections, csv, io, itertools, os, sys

ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# the thresholds the manuscript reports, and the ranges swept around them
TM, RMSD, MAD = 0.6, 2.0, 1.5
TM_RANGE, RMSD_RANGE, MAD_RANGE = (0.5, 0.7), (1.5, 3.0), (1.0, 2.5)
STEPS = 5  # 5 x 5 x 5 = the 125 combinations reported in section 3.10


def read_tsv(rel):
    with io.open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def norm(target):
    """Target names carry a .pdb/.pdb.gz suffix in some tables and not others."""
    return target.replace(".pdb.gz", "").replace(".pdb", "")


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def linspace(lo, hi, n):
    return [lo + (hi - lo) * i / (n - 1) for i in range(n)]


def load():
    cand = {(r["query"], norm(r["target"])): r for r in read_tsv("screen_v1/results/candidates.tsv")}
    site = {(r["query"], norm(r["target"])): r for r in read_tsv("screen_v2/results/active_site_all_hits.tsv")}
    tier1 = {(r["query"], norm(r["target"])) for r in read_tsv("screen_v2/results/tier1_candidates.tsv")}
    return cand, site, tier1


def passes(key, cand, site, tm_min=None, rmsd_max=None, mad_max=None,
           need_phmmer_blind=False, need_catalytic=True):
    """Apply the tier-1 rule, with any subset of its criteria."""
    c, s = cand[key], site.get(key)
    if s is None:
        return False
    if need_catalytic:
        # every curated catalytic residue present; queries with none annotated cannot qualify
        n_cat, n_identical = num(s["n_catalytic"]), num(s["n_catalytic_identical"])
        if not n_cat or n_cat < 1 or n_identical is None or n_identical < n_cat:
            return False
    rmsd, mad, tm = num(s["site_rmsd"]), num(s["site_dist_mad"]), num(c["tm_score_qnorm"])
    if rmsd_max is not None and (rmsd is None or rmsd > rmsd_max):
        return False
    if mad_max is not None and (mad is None or mad > mad_max):
        return False
    if tm_min is not None and (tm is None or tm < tm_min):
        return False
    if need_phmmer_blind and c["phmmer_evalue"].strip() not in ("", "NA", "nan", "None"):
        # phmmer must return no hit at all from the query, not even at E = 10
        return False
    return True


def shortlist(cand, site, tm_min, rmsd_max, mad_max):
    return {k for k in cand
            if passes(k, cand, site, tm_min, rmsd_max, mad_max, need_phmmer_blind=True)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "screen_v2", "evaluation"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    cand, site, tier1 = load()

    # attrition, as reported in section 3.10
    step1 = [k for k in cand if passes(k, cand, site)]
    step2 = [k for k in step1 if passes(k, cand, site, None, RMSD, MAD)]
    step3 = [k for k in step2 if passes(k, cand, site, TM, RMSD, MAD)]
    step4 = [k for k in step3 if passes(k, cand, site, TM, RMSD, MAD, need_phmmer_blind=True)]
    attrition = [
        ("candidate pairs", len(cand), 1643),
        ("keep every annotated catalytic residue", len(step1), 30),
        ("inside the geometric thresholds", len(step2), 17),
        ("also TM-score >= %.1f" % TM, len(step3), 15),
        ("also phmmer-undetected (the shortlist)", len(step4), 13),
    ]
    print("Attrition to the tier-1 shortlist")
    for label, got, expect in attrition:
        print("  %-40s %5d   (reported: %d)%s"
              % (label, got, expect, "" if got == expect else "   <-- MISMATCH"))
    agrees = set(step4) == tier1
    print("  shortlist identical to tier1_candidates.tsv: %s" % agrees)

    # the sweep itself
    grid = list(itertools.product(linspace(*TM_RANGE, STEPS),
                                  linspace(*RMSD_RANGE, STEPS),
                                  linspace(*MAD_RANGE, STEPS)))
    counts, sizes, rows = collections.Counter(), [], []
    for tm_min, rmsd_max, mad_max in grid:
        sl = shortlist(cand, site, tm_min, rmsd_max, mad_max)
        sizes.append(len(sl))
        rows.append((tm_min, rmsd_max, mad_max, len(sl)))
        for k in sl:
            counts[k] += 1
    n = len(grid)
    robust = sorted(t for q, t in tier1 if counts[(q, t)] == n)
    print("\nSweep over %d threshold combinations" % n)
    print("  shortlist size: %d to %d" % (min(sizes), max(sizes)))
    print("  of the %d shortlisted proteins, in every setting: %d" % (len(tier1), len(robust)))
    print("  %s" % ", ".join(robust))

    # which criteria actually bind at the reported setting
    print("\nWhich criteria bind")
    for mad_max in linspace(*MAD_RANGE, STEPS):
        print("  MAD  <= %5.3f (TM %.2f, RMSD %.2f) -> %2d"
              % (mad_max, TM, RMSD, len(shortlist(cand, site, TM, RMSD, mad_max))))
    for rmsd_max in linspace(*RMSD_RANGE, STEPS):
        print("  RMSD <= %5.3f (TM %.2f, MAD %.2f)  -> %2d"
              % (rmsd_max, TM, MAD, len(shortlist(cand, site, TM, rmsd_max, MAD))))
    for tm_min in linspace(*TM_RANGE, STEPS):
        print("  TM   >= %5.3f (RMSD %.2f, MAD %.2f) -> %2d"
              % (tm_min, RMSD, MAD, len(shortlist(cand, site, tm_min, RMSD, MAD))))

    grid_path = os.path.join(args.out, "threshold_sweep_grid.tsv")
    with io.open(grid_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["tm_min", "site_rmsd_max", "site_mad_max", "n_shortlisted"])
        for tm_min, rmsd_max, mad_max, k in rows:
            w.writerow(["%.4f" % tm_min, "%.4f" % rmsd_max, "%.4f" % mad_max, k])

    per_path = os.path.join(args.out, "threshold_sweep_per_protein.tsv")
    with io.open(per_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["target", "query", "n_settings_shortlisted", "n_settings", "fraction"])
        for q, t in sorted(tier1, key=lambda k: (-counts[k], k[1])):
            w.writerow([t, q, counts[(q, t)], n, "%.3f" % (counts[(q, t)] / n)])

    print("\nWrote %s\n      %s" % (grid_path, per_path))
    return 0 if agrees and all(g == e for _, g, e in attrition) else 1


if __name__ == "__main__":
    sys.exit(main())
