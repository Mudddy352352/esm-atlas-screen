#!/usr/bin/env python3
"""What a real nr search says about the thirteen shortlisted proteins.

Section 3.8 answers this with a BLAST-class aligner run on local sequence sets
and E-values scaled to the size of nr, because no sequence service was
reachable when the screen ran. These are the searches themselves. For each
protein the question is the one the projection stood in for: does a pairwise
search against nr link it to the query's specific activity, or only to the
family?

A hit counts as naming the query activity if its description names the enzyme
the query is (beta-glucuronidase, cytidine deaminase). A hit counts as
family-level if it names the fold or family without an activity, which for the
glucuronidase queries means glycoside hydrolase family 2 and its relatives.

Usage: analyse_nr.py <work_dir>
"""
import glob
import os
import re
import sys
import xml.etree.ElementTree as ET

import pandas as pd

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
_ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

ACTIVITY = {
    "beta-glucuronidase": r"glucuronidase",
    "cytidine deaminase": r"cytidine deaminase|cytosine nucleoside deaminase",
}
GENERIC = r"hypothetical|uncharacteri|DUF\d|domain-containing|glycoside hydrolase family|glycosyl hydrolase|MULTISPECIES: *$"
QUERY_ACTIVITY = {
    "P05804": "beta-glucuronidase", "uidA": "beta-glucuronidase",
    "gus": "beta-glucuronidase", "glucuronidase": "beta-glucuronidase",
    "P0ABF6": "cytidine deaminase", "P19079": "cytidine deaminase",
    "cdd": "cytidine deaminase", "deaminase": "cytidine deaminase",
}


def activity_of(query):
    q = str(query).lower()
    for key, act in QUERY_ACTIVITY.items():
        if key.lower() in q:
            return act
    return None


def hits(path):
    root = ET.parse(path).getroot()
    out = []
    for h in root.findall(".//Hit"):
        hsp = h.find(".//Hsp")
        if hsp is None:
            continue
        al = int(hsp.findtext("Hsp_align-len") or 0) or 1
        out.append(dict(
            desc=h.findtext("Hit_def", ""),
            evalue=float(hsp.findtext("Hsp_evalue")),
            pct_id=100.0 * int(hsp.findtext("Hsp_identity") or 0) / al,
            align_len=al))
    return out


def main(work):
    res = os.path.join(work, "results")
    tier1 = pd.read_csv(_ROOT + "/screen_v2/results/tier1_candidates.tsv", sep="\t")
    qmap = dict(zip(tier1.target, tier1["query"]))

    rows = []
    for path in sorted(glob.glob(os.path.join(res, "MGYP*.xml"))):
        acc = os.path.basename(path).split(".")[0]
        hs = hits(path)
        if not hs:
            rows.append(dict(protein=acc, n_hits=0))
            continue
        act = activity_of(qmap.get(acc, ""))
        pat = ACTIVITY.get(act, r"$^")
        named = [h for h in hs if re.search(pat, h["desc"], re.I)]
        best = min(hs, key=lambda h: h["evalue"])
        best_named = min(named, key=lambda h: h["evalue"]) if named else None
        rows.append(dict(
            protein=acc, query=qmap.get(acc, ""), query_activity=act, n_hits=len(hs),
            best_evalue=best["evalue"], best_pct_id=round(best["pct_id"], 1),
            best_desc=best["desc"][:80],
            best_is_generic=bool(re.search(GENERIC, best["desc"], re.I)),
            n_naming_activity=len(named),
            named_best_evalue=best_named["evalue"] if best_named else None,
            named_best_pct_id=round(best_named["pct_id"], 1) if best_named else None,
            named_best_desc=best_named["desc"][:80] if best_named else ""))

    df = pd.DataFrame(rows).sort_values("protein")
    df.to_csv(os.path.join(res, "nr_shortlist_summary.tsv"), sep="\t", index=False)

    print("Thirteen shortlisted proteins against nr\n")
    print("%-20s %-9s %7s %8s  %s" % ("protein", "best E", "best id", "activity", "best hit"))
    for _, r in df.iterrows():
        print("%-20s %-9.1e %6.1f%% %8s  %s"
              % (r.protein, r.best_evalue, r.best_pct_id,
                 "yes" if r.n_naming_activity else "no", r.best_desc[:56]))
    print()
    named = df[df.n_naming_activity > 0]
    print("proteins with any hit naming the query activity: %d of %d" % (len(named), len(df)))
    if len(named):
        print("%-20s %-10s %8s  %s" % ("protein", "E", "identity", "closest hit naming it"))
        for _, r in named.iterrows():
            print("%-20s %-10.1e %7.1f%%  %s"
                  % (r.protein, r.named_best_evalue, r.named_best_pct_id,
                     r.named_best_desc[:52]))
    print()
    print("proteins whose single best hit is a generic or unannotated name: %d of %d"
          % (int(df.best_is_generic.sum()), len(df)))
    print("median identity to the best nr hit: %.1f%%" % df.best_pct_id.median())
    print("proteins with a nr hit at 95%% identity or better: %d"
          % int((df.best_pct_id >= 95).sum()))


if __name__ == "__main__":
    main(*sys.argv[1:2])
