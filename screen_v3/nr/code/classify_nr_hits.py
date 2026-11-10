#!/usr/bin/env python3
"""Classify what each shortlisted protein's nr hits actually say.

A hit description from nr is several titles joined by ">", because identical
sequences are merged, so "names the activity" has to be asked title by title.
Three things are counted separately, since they mean different things:

  unambiguous   a title naming the query's activity and nothing else
                ("beta-glucuronidase [Organism]")
  ambiguous     a dual annotation that names the ambiguity rather than
                resolving it ("beta-galactosidase/beta-glucuronidase")
  query enzyme  the characterised enzyme the screen started from

Usage: classify_nr_hits.py <work_dir>
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

ACT = {"beta-glucuronidase": r"beta[- ]?glucuronidase",
       "cytidine deaminase": r"cytidine deaminase"}
DUAL = r"galactosidase\s*/\s*beta[- ]?glucuronidase|beta[- ]?glucuronidase\s*/\s*.*galactosidase"
QUERY_ORG = {"beta-glucuronidase": r"Escherichia coli",
             "cytidine deaminase": r"Escherichia coli"}


def activity_of(query):
    q = str(query).lower()
    if "glucuronid" in q or "gus" in q or "uida" in q:
        return "beta-glucuronidase"
    if "deaminase" in q or "cdd" in q:
        return "cytidine deaminase"
    return None


def main(work):
    res = os.path.join(work, "results")
    tier1 = pd.read_csv(_ROOT + "/screen_v2/results/tier1_candidates.tsv", sep="\t")
    qmap = dict(zip(tier1.target, tier1["query"]))
    rows = []
    for path in sorted(glob.glob(os.path.join(res, "deep_MGYP*.xml"))):
        acc = os.path.basename(path)[5:-4]
        act = activity_of(qmap.get(acc, ""))
        pat, dual_pat = ACT.get(act, r"$^"), DUAL
        root = ET.parse(path).getroot()
        n = 0
        best_u = best_d = best_q = None
        for h in root.findall(".//Hit"):
            hsp = h.find(".//Hsp")
            if hsp is None:
                continue
            n += 1
            desc = h.findtext("Hit_def", "")
            if not re.search(pat, desc, re.I):
                continue
            al = int(hsp.findtext("Hsp_align-len") or 1)
            ev = float(hsp.findtext("Hsp_evalue"))
            pid = round(100.0 * int(hsp.findtext("Hsp_identity")) / al, 1)
            titles = [t.strip() for t in desc.split(">")]
            for t in titles:
                if not re.search(pat, t, re.I):
                    continue
                if re.search(dual_pat, t, re.I):
                    if best_d is None or ev < best_d[0]:
                        best_d = (ev, pid, t[:90])
                else:
                    if best_u is None or ev < best_u[0]:
                        best_u = (ev, pid, t[:90])
                    if re.search(QUERY_ORG.get(act, r"$^"), t) and (
                            best_q is None or ev < best_q[0]):
                        best_q = (ev, pid, t[:90])
        rows.append(dict(
            protein=acc, activity=act, hits_examined=n,
            unambiguous_evalue=best_u[0] if best_u else None,
            unambiguous_pct_id=best_u[1] if best_u else None,
            unambiguous_title=best_u[2] if best_u else "",
            ambiguous_evalue=best_d[0] if best_d else None,
            ambiguous_pct_id=best_d[1] if best_d else None,
            query_organism_hit=best_q[2] if best_q else "",
            query_organism_evalue=best_q[0] if best_q else None,
            query_organism_pct_id=best_q[1] if best_q else None))
    df = pd.DataFrame(rows).sort_values("protein")
    df.to_csv(os.path.join(res, "nr_activity_classified.tsv"), sep="\t", index=False)

    u = df[df.unambiguous_evalue.notna() & (df.unambiguous_evalue <= 1e-3)]
    print("%-20s %-11s %6s  %s" % ("protein", "best E", "id", "closest unambiguous title"))
    for _, r in df.iterrows():
        print("%-20s %-11s %5s%%  %s" % (
            r.protein,
            "%.1e" % r.unambiguous_evalue if r.unambiguous_evalue is not None else "none",
            r.unambiguous_pct_id if r.unambiguous_pct_id is not None else "-",
            (r.unambiguous_title or "only a dual galactosidase/glucuronidase name")[:58]))
    print("\nlinked to an unambiguously named enzyme of the query activity at E <= 1e-3: %d of %d"
          % (len(u), len(df)))
    print("median identity to that hit: %.1f%%" % u.unambiguous_pct_id.median())
    q = df[df.query_organism_evalue.notna()]
    print("hits to the characterised query enzyme's own organism: %d, best E %s"
          % (len(q), "%.2g" % q.query_organism_evalue.min() if len(q) else "none"))
    print("of those, significant at E <= 1e-3: %d"
          % int((q.query_organism_evalue <= 1e-3).sum()))


if __name__ == "__main__":
    main(*sys.argv[1:2])
