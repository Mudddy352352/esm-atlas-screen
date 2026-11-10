#!/usr/bin/env python3
"""Does a pairwise search link each shortlisted protein to the query activity?

The plain nr searches cannot answer this cleanly. A shortlisted protein has
thousands of same-family hits at E = 0, so a weak alignment to a protein
actually annotated with the query's activity falls outside any sane hit list.
This asks the question directly: the same blastp against nr, restricted by an
Entrez query to proteins whose name carries the query's activity, with the
E-values still computed against nr. That is the measurement section 3.8
projected.

Submits, then collects; both are safe to rerun, and ids are kept in
probe_state.json so nothing is submitted twice.

Usage: nr_activity_probe.py <work_dir> submit|collect
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_nr_blast import BLAST, TOOL, sequences  # noqa: E402

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
_ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

TIER1 = _ROOT + "/screen_v2/results/tier1_candidates.tsv"
ENTREZ = {
    "beta-glucuronidase": "beta-glucuronidase[Protein Name]",
    "cytidine deaminase": "cytidine deaminase[Protein Name]",
}


def activity_of(query):
    q = str(query).lower()
    if "glucuronid" in q or "gus" in q or "uida" in q or "p05804" in q:
        return "beta-glucuronidase"
    if "deaminase" in q or "cdd" in q:
        return "cytidine deaminase"
    return None


def call(params, timeout=300):
    req = urllib.request.Request(BLAST, data=urllib.parse.urlencode(params).encode())
    return urllib.request.urlopen(req, timeout=timeout).read().decode("utf8", "replace")


def main(work, action):
    res = os.path.join(work, "results")
    path = os.path.join(res, "probe_state.json")
    state = json.load(open(path)) if os.path.exists(path) else {}
    tier1 = pd.read_csv(TIER1, sep="\t")
    seqs = sequences()

    if action == "submit":
        for _, row in tier1.iterrows():
            acc = row["target"]
            act = activity_of(row["query"])
            if not act or state.get(acc, {}).get("rid") or acc not in seqs:
                continue
            txt = call({"CMD": "Put", "PROGRAM": "blastp", "DATABASE": "nr",
                        "QUERY": seqs[acc], "HITLIST_SIZE": "50", "EXPECT": "100",
                        "ENTREZ_QUERY": ENTREZ[act], "TOOL": TOOL}, timeout=180)
            m = re.search(r"RID = (\S+)", txt)
            state[acc] = dict(rid=m.group(1) if m else None, activity=act, done=False)
            json.dump(state, open(path, "w"), indent=1)
            print("[%s] %s -> %s" % (acc, act, state[acc]["rid"]), flush=True)
            time.sleep(12)
        return

    rows, waiting = [], []
    for acc, rec in sorted(state.items()):
        rid = rec.get("rid")
        if not rid:
            continue
        info = call({"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid,
                     "TOOL": TOOL}, timeout=120)
        if "Status=WAITING" in info:
            waiting.append(acc)
            continue
        hits_present = "ThereAreHits=yes" in info
        best = None
        if hits_present:
            xml = call({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML",
                        "ALIGNMENTS": "50", "DESCRIPTIONS": "50", "TOOL": TOOL})
            open(os.path.join(res, "probe_%s.xml" % acc), "w").write(xml)
            import xml.etree.ElementTree as ET
            root = ET.fromstring(xml)
            for h in root.findall(".//Hit"):
                hsp = h.find(".//Hsp")
                if hsp is None:
                    continue
                al = int(hsp.findtext("Hsp_align-len") or 1)
                cand = dict(desc=h.findtext("Hit_def", "")[:90],
                            evalue=float(hsp.findtext("Hsp_evalue")),
                            pct_id=round(100.0 * int(hsp.findtext("Hsp_identity")) / al, 1))
                if best is None or cand["evalue"] < best["evalue"]:
                    best = cand
        rec.update(done=True)
        state[acc] = rec
        rows.append(dict(protein=acc, activity=rec["activity"],
                         best_evalue=best["evalue"] if best else None,
                         pct_id=best["pct_id"] if best else None,
                         significant=bool(best and best["evalue"] <= 1e-3),
                         best_hit=best["desc"] if best else "no hit"))
        time.sleep(3)
    json.dump(state, open(path, "w"), indent=1)
    if rows:
        df = pd.DataFrame(rows).sort_values("protein")
        df.to_csv(os.path.join(res, "nr_activity_probe.tsv"), sep="\t", index=False)
        print("%-20s %-20s %-11s %7s  %s"
              % ("protein", "activity searched", "best E", "id", "best hit"))
        for _, r in df.iterrows():
            print("%-20s %-20s %-11s %6s%%  %s"
                  % (r.protein, r.activity,
                     "%.1e" % r.best_evalue if r.best_evalue is not None else "none",
                     r.pct_id if r.pct_id is not None else "-", str(r.best_hit)[:48]))
        print("\nlinked at E <= 1e-3: %d of %d" % (int(df.significant.sum()), len(df)))
    print("still running: %d" % len(waiting))


if __name__ == "__main__":
    main(*sys.argv[1:3])
