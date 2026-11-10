#!/usr/bin/env python3
"""The same nr searches with a deep hit list, which is how the activity
question can actually be asked.

A shortlisted protein has thousands of same-family hits at E = 0, so a weak
alignment to a protein annotated with the query's own activity falls far
outside a hundred-hit list. Restricting the search with an Entrez query looked
like the fix and is not: on this endpoint every restricted search returns
nothing, including E. coli beta-glucuronidase searched against nr restricted to
beta-glucuronidase, which is the control that settles it. So the hit list is
simply taken deep instead, with the statistics still those of nr.

Usage: nr_deep.py <work_dir> submit|collect
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_nr_blast import BLAST, TOOL, sequences  # noqa: E402

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
_ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

TIER1 = _ROOT + "/screen_v2/results/tier1_candidates.tsv"
DEPTH = "5000"
PATTERN = {"beta-glucuronidase": r"glucuronidase",
           "cytidine deaminase": r"cytidine deaminase"}


def activity_of(query):
    q = str(query).lower()
    if "glucuronid" in q or "gus" in q or "uida" in q:
        return "beta-glucuronidase"
    if "deaminase" in q or "cdd" in q:
        return "cytidine deaminase"
    return None


def call(d, t=600):
    return urllib.request.urlopen(
        urllib.request.Request(BLAST, data=urllib.parse.urlencode(d).encode()),
        timeout=t).read().decode("utf8", "replace")


def main(work, action):
    res = os.path.join(work, "results")
    path = os.path.join(res, "deep_state.json")
    state = json.load(open(path)) if os.path.exists(path) else {}
    tier1 = pd.read_csv(TIER1, sep="\t")
    seqs = sequences()

    if action == "submit":
        for _, row in tier1.iterrows():
            acc = row["target"]
            if state.get(acc, {}).get("rid") or acc not in seqs:
                continue
            txt = call({"CMD": "Put", "PROGRAM": "blastp", "DATABASE": "nr",
                        "QUERY": seqs[acc], "HITLIST_SIZE": DEPTH,
                        "DESCRIPTIONS": DEPTH, "ALIGNMENTS": "0",
                        "EXPECT": "10", "TOOL": TOOL}, t=180)
            m = re.search(r"RID = (\S+)", txt)
            state[acc] = dict(rid=m.group(1) if m else None,
                              activity=activity_of(row["query"]), done=False)
            json.dump(state, open(path, "w"), indent=1)
            print("[%s] %s" % (acc, state[acc]["rid"]), flush=True)
            time.sleep(12)
        return

    rows, waiting = [], []
    for acc, rec in sorted(state.items()):
        if rec.get("done"):
            continue
        info = call({"CMD": "Get", "FORMAT_OBJECT": "SearchInfo",
                     "RID": rec["rid"], "TOOL": TOOL}, t=120)
        if "Status=WAITING" in info:
            waiting.append(acc)
            continue
        xml = call({"CMD": "Get", "RID": rec["rid"], "FORMAT_TYPE": "XML",
                    "DESCRIPTIONS": DEPTH, "ALIGNMENTS": DEPTH, "TOOL": TOOL})
        open(os.path.join(res, "deep_%s.xml" % acc), "w").write(xml)
        rec.update(done=True)
        json.dump(state, open(path, "w"), indent=1)
        print("[%s] collected" % acc, flush=True)
        time.sleep(3)

    for acc, rec in sorted(state.items()):
        f = os.path.join(res, "deep_%s.xml" % acc)
        if not os.path.exists(f):
            continue
        root = ET.parse(f).getroot()
        pat = PATTERN.get(rec.get("activity") or "", r"$^")
        n, best = 0, None
        for h in root.findall(".//Hit"):
            hsp = h.find(".//Hsp")
            if hsp is None:
                continue
            n += 1
            d = h.findtext("Hit_def", "")
            if re.search(pat, d, re.I):
                al = int(hsp.findtext("Hsp_align-len") or 1)
                c = dict(desc=d[:100], evalue=float(hsp.findtext("Hsp_evalue")),
                         pct_id=round(100.0 * int(hsp.findtext("Hsp_identity")) / al, 1),
                         ambiguous=bool(re.search(r"galactosidase/|/beta-glucuronidase|beta-galactosidase/beta-glucuronidase", d, re.I)))
                if best is None or c["evalue"] < best["evalue"]:
                    best = c
        rows.append(dict(protein=acc, activity=rec.get("activity"), n_hits=n,
                         best_named_evalue=best["evalue"] if best else None,
                         best_named_pct_id=best["pct_id"] if best else None,
                         ambiguous_name=best["ambiguous"] if best else None,
                         best_named_desc=best["desc"] if best else "none in %d hits" % n))
    if rows:
        df = pd.DataFrame(rows).sort_values("protein")
        df.to_csv(os.path.join(res, "nr_deep_activity.tsv"), sep="\t", index=False)
        print("\n%-20s %6s %-11s %7s  %s"
              % ("protein", "hits", "best E", "id", "closest hit naming the activity"))
        for _, r in df.iterrows():
            print("%-20s %6d %-11s %6s  %s"
                  % (r.protein, r.n_hits,
                     "%.1e" % r.best_named_evalue if r.best_named_evalue is not None else "none",
                     r.best_named_pct_id if r.best_named_pct_id is not None else "-",
                     str(r.best_named_desc)[:60]))
        sig = df[(df.best_named_evalue.notna()) & (df.best_named_evalue <= 1e-3)]
        print("\nlinked to the query activity at E <= 1e-3: %d of %d" % (len(sig), len(df)))
    print("still running: %d" % len(waiting))


if __name__ == "__main__":
    main(*sys.argv[1:3])
